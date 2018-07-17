#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_min_max.h>

#define CHUNK_SIZE 5

/*********define cluster************/
typedef struct cluster cluster;
struct cluster{
  int ndim;
  int nstrt,nsize,ndpth,nnson,nmbr;
  int ndscd;                         //number of descendants                                                                    
  double *bmin;                      //bounding box                                                                             
  double *bmax;
  double zwdth;
  cluster **pc_sons;
};

/*********define leafMtx***********/ //information of H-matrix
typedef struct leafmtx leafmtx;
struct leafmtx{
  int ltmtx;                         //kind of the matrix; 1:rk 2:full
  int kt;                            //rank of the partition
  int nstrtl,ndl;                    //the coordination of the first element of partition
  int nstrtt,ndt;                    //the length & width of the partition
  double *a1,*a2;                    //the elements of the partition
};

/*********define leafmtxp*********/  //whole H-matrix
typedef struct leafmtxp leafmtxp;
struct leafmtxp{
  int nlf;                           //number of partitions
  int nlfkt;                         //number ot partitions approximated
};

void supermatrix_construction_cog_leafmtrx(leafmtxp *st_leafmtxp,double (*gmid)[3],double param[],int *lod,int *lnmtx,int nofc,int nffc,int ndim);

void random(int *lod,int nd);
int contains(int *lod, int n, int x);

cluster * create_ctree_ssgeom(cluster *st_clt,double (*zgmid)[3],double param[],int *lod,int ndpth,int ndscd,int nsrt,int nd,int md,int ndim,int nclst);
void mm(double (*zgmid)[3],int *lod,int nd,int md, int ndim);
void mm1(double (*zgmid)[3],int *lod,int nd,int md, int ndim);
void mm2(double (*zgmid)[3],int *lod,int nd,int md, int ndim);
void swap1(double (*zgmid)[3],int *lod,int nd,int md, int ndim);
void swap2(double (*zgmid)[3],int *lod,int *templod,int nd,int md, int ndim);
double get_wall_time();
double get_cpu_time();

void checkClusterTree(FILE *f,cluster *st_clt);

int ncall[500];
int countlist[500];
double ntime[2][36];
int cs[3]={0,0,0};

int main(int argc, char **argv){
  if(argc >= 2){
    __cilkrts_set_param("nworkers",argv[1]);
  }
  int nworkers = __cilkrts_get_nworkers();
  printf("number of workers:%d\n",nworkers);
  /******** read file *********/
  FILE *file;
  int countOfNode=0;
  int count = 0;
  int i;
  double (*coordOfNode)[3];
  double (*coordOfFace)[3];
  file = fopen("bem_data/input_50ms.txt","r");
  if(file == NULL){
    printf("Error: Unable to input file 'input_50ms.txt'!\n");
  }else{
    char line[100];
    fgets(line,sizeof(line),file);
    sscanf(line, "%d", &countOfNode);
    coordOfNode = (double(*)[3])malloc(countOfNode*3*sizeof(double));
    for(i=0;i<countOfNode;i++){
      fgets(line,sizeof(line),file);
      double x,y,z;
      sscanf(line,"%lf %lf %lf",&x,&y,&z);
      coordOfNode[i][0] = x;
      coordOfNode[i][1] = y;
      coordOfNode[i][2] = z;
    }
    fgets(line,sizeof(line),file);
    sscanf(line,"%d",&count);
    coordOfFace = (double(*)[3])malloc(count*3*sizeof(double));

    fgets(line,sizeof(line),file);
    fgets(line,sizeof(line),file);
    fgets(line,sizeof(line),file);
    
    printf("count:%d\n",count);
    for(i=0;i<count;i++){
      fgets(line,sizeof(line),file);
      int X,Y,Z;
      sscanf(line,"%d %d %d",&X,&Y,&Z);
      coordOfFace[i][0] = (coordOfNode[X][0] + coordOfNode[Y][0] + coordOfNode[Z][0])/3;
      coordOfFace[i][1] = (coordOfNode[X][1] + coordOfNode[Y][1] + coordOfNode[Z][1])/3;
      coordOfFace[i][2] = (coordOfNode[X][2] + coordOfNode[Y][2] + coordOfNode[Z][2])/3;
    }
  }
  fclose(file);
  free(coordOfNode);
  double param[100];
  for(i=0;i<100;i++){
    param[i] = 0;
  }
  param[21] = 10.0;
  param[31] = 1.1;
  param[41] = 15.0;
  param[51] = 2.0;

  /***********clustering***********/
  leafmtxp *st_leafmtxp;
  int *lod;
  int *lnmtx;
  int nofc = count;
  int nffc = 1;
  int ndim = 3;
  lnmtx = (int *)malloc(3 * sizeof(int));
  for(i=0;i<3;i++){
    lnmtx[i] = 0;
  }
  st_leafmtxp = (leafmtxp *)malloc(sizeof(leafmtxp));
  lod = (int *)malloc(nofc*sizeof(int));
  for(i=0;i<nofc;i++){
    lod[i] = 0;
  }

  /*********construction of array of leafmatrices******/
  supermatrix_construction_cog_leafmtrx(st_leafmtxp,coordOfFace,param,lod,lnmtx,nofc,nffc,ndim);

  printf("Execution complete successfully!\n");
  return 0;
}

void supermatrix_construction_cog_leafmtrx(leafmtxp *st_leafmtxp,   //the H-matrix
					   double (*gmid)[3],       //coordination of objects
					   double param[],int *lod,
					   int *lnmtx,              //1:k-rank 2:dense 3:H-matrix
					   int nofc,int nffc,       //number of elements in same coordination
					   int ndim){
  /*****init****/
  cluster *st_clt = (cluster *)malloc(sizeof(cluster));
  int i,j,k,nfl,nflkt,ip,il,ig;
  int nd = nofc * nffc;
  int *lodfc,*templod;
  int lel;
  leafmtx **st_leafmtx;
  lodfc = (int *)malloc(nofc*sizeof(int));
  templod = (int *)malloc(nofc*sizeof(int));
  for(il=0;il<nofc;il++){
    lodfc[il] = il;
  }
  int nsrt = 0;
  int ndf = nofc;
  int nclst = 0;
  int ndpth = 0;
  int ndscd = 0;

  double start,end,spent;
  random(lodfc,ndf);

  start = get_wall_time();
  //swap1(gmid,lodfc,ndf,nofc,ndim);
  swap2(gmid,lodfc,templod,ndf,nofc,ndim);
  end = get_wall_time();
  spent = end - start;
  //printf("sequtial time spent:%lf\n",spent);
  start = get_wall_time();
  /*mm1(gmid,lodfc,ndf,nofc, ndim);
  end = get_wall_time();
  spent = end - start;
  printf("jump access parallelization time spent:%.10f\n",spent);
  start = get_wall_time();
  mm2(gmid,lodfc,ndf,nofc, ndim);
  end = get_wall_time();
  spent = end - start;
  printf("succession access parallelization time spent:%.10f\n",spent);*/
}


void mm1(double (*zgmid)[3],int *lod,int nd,int md, int ndim){
  int id,il;
  double zlmin[ndim],zlmax[ndim];
  CILK_C_REDUCER_MAX(max,double,-10.1);
  CILK_C_REDUCER_MIN(min,double,10.1);
  CILK_C_REGISTER_REDUCER(max);
  CILK_C_REGISTER_REDUCER(min);
  _Cilk_for(id=0;id<ndim;id++){
    REDUCER_VIEW(min) = 10.1;
    REDUCER_VIEW(max) = -10.1;
    _Cilk_for(il=0;il<nd;il++){
      CILK_C_REDUCER_MIN_CALC(min,zgmid[lod[il]][id]);
      CILK_C_REDUCER_MAX_CALC(max,zgmid[lod[il]][id]);
    }
    zlmin[id] = REDUCER_VIEW(min);
    zlmax[id] = REDUCER_VIEW(max);
  }
  CILK_C_UNREGISTER_REDUCER(max);
  CILK_C_UNREGISTER_REDUCER(min);
  for(id=0;id<ndim;id++){
    printf("%lf %lf\n",zlmin[id],zlmax[id]);
  }
}

void mm(double (*zgmid)[3],int *lod,int nd,int md, int ndim){
  int id,il;
  double zlmin[ndim],zlmax[ndim];
  for(id=0;id<ndim;id++){
    zlmin[id] = zgmid[lod[0]][id];
    zlmax[id] = zlmin[id];
    for(il=1;il<nd;il++){
      double zg = zgmid[lod[il]][id];
      if(zg < zlmin[id]){
        zlmin[id] = zg;
      }else if(zlmax[id] < zg){
        zlmax[id] = zg;
      }
    }
  }
}

void swap1(double (*zgmid)[3],int *lod,int nd,int md, int ndim){
  int id,il;
  double zlmin[ndim],zlmax[ndim];
  for(id=0;id<ndim;id++){
    zlmin[id] = zgmid[lod[0]][id];
    zlmax[id] = zlmin[id];
    for(il=1;il<nd;il++){
      double zg = zgmid[lod[il]][id];
      if(zg < zlmin[id]){
        zlmin[id] = zg;
      }else if(zlmax[id] < zg){
        zlmax[id] = zg;
      }
    }
  }
  double zcoef = 1.1;
  double zdiff = zlmax[0] - zlmin[0];
  int ncut = 0;
  for(id=0;id<ndim;id++){
    double zidiff = zlmax[id]-zlmin[id];
    if(zidiff > zcoef * zdiff){
      zdiff = zidiff;
      ncut = id;
    }
  }
  double zlmid = 0.5 * (zlmax[ncut] + zlmin[ncut]);
  int nl = 0;
  int nr = nd-1;
  double start = get_wall_time();
  while(nl < nr){
    while(nl < nd && zgmid[lod[nl]][ncut] <= zlmid){
      nl = nl + 1;
    }
    while(nr >= 0 && zgmid[lod[nr]][ncut] > zlmid){
      nr = nr - 1;
    }
    if(nl < nr){
      int nh = lod[nl];
      lod[nl] = lod[nr];
      lod[nr] = nh;
    }
  }
  double end = get_wall_time();
  double spent = end - start;
  printf("sequtial time spent:%lf\n",spent);
}

void swap2(double (*zgmid)[3],int *lod,int *templod,int nd,int md, int ndim){
  int id,il;
  double zlmin[ndim],zlmax[ndim];
  CILK_C_REDUCER_MAX(max_x,double,-10.1);
  CILK_C_REDUCER_MIN(min_x,double,10.1);
  CILK_C_REDUCER_MAX(max_y,double,-10.1);
  CILK_C_REDUCER_MIN(min_y,double,10.1);
  CILK_C_REDUCER_MAX(max_z,double,-10.1);
  CILK_C_REDUCER_MIN(min_z,double,10.1);
  CILK_C_REGISTER_REDUCER(max_x);
  CILK_C_REGISTER_REDUCER(min_x);
  CILK_C_REGISTER_REDUCER(max_y);
  CILK_C_REGISTER_REDUCER(min_y);
  CILK_C_REGISTER_REDUCER(max_z);
  CILK_C_REGISTER_REDUCER(min_z);
#pragma cilk grainsize = 128
  _Cilk_for(il=0;il<nd;il++){
    CILK_C_REDUCER_MIN_CALC(min_x,zgmid[lod[il]][0]);
    CILK_C_REDUCER_MAX_CALC(max_x,zgmid[lod[il]][0]);
    CILK_C_REDUCER_MIN_CALC(min_y,zgmid[lod[il]][1]);
    CILK_C_REDUCER_MAX_CALC(max_y,zgmid[lod[il]][1]);
    CILK_C_REDUCER_MIN_CALC(min_z,zgmid[lod[il]][2]);
    CILK_C_REDUCER_MAX_CALC(max_z,zgmid[lod[il]][2]);
  }
  zlmin[0] = REDUCER_VIEW(min_x);
  zlmax[0] = REDUCER_VIEW(max_x);
  zlmin[1] = REDUCER_VIEW(min_y);
  zlmax[1] = REDUCER_VIEW(max_y);
  zlmin[2] = REDUCER_VIEW(min_z);
  zlmax[2] = REDUCER_VIEW(max_z);
  CILK_C_UNREGISTER_REDUCER(max_x);
  CILK_C_UNREGISTER_REDUCER(min_x);
  CILK_C_UNREGISTER_REDUCER(max_y);
  CILK_C_UNREGISTER_REDUCER(min_y);
  CILK_C_UNREGISTER_REDUCER(max_z);
  CILK_C_UNREGISTER_REDUCER(min_z);
  double zcoef = 1.1;
  double zdiff = zlmax[0] - zlmin[0];
  int ncut = 0;
  for(id=0;id<ndim;id++){
    double zidiff = zlmax[id]-zlmin[id];
    if(zidiff > zcoef * zdiff){
      zdiff = zidiff;
      ncut = id;
    }
  }
  double zlmid = 0.5 * (zlmax[ncut] + zlmin[ncut]);
  int nl = 0;
  int nr = nd-1;
  double start = get_wall_time();
  int gn = (nd/CHUNK_SIZE)+1;
  int *lessNum = (int *)malloc(gn*sizeof(int));
  int *moreNum = (int *)malloc(gn*sizeof(int));
  int *lessStart = (int *)malloc(gn*sizeof(int));
  int *moreStart = (int *)malloc(gn*sizeof(int));
  for(id=0;id<gn;id++){
    lessNum[id] = 0;
    moreNum[id] = 0;
    lessStart[id] = 0;
    moreStart[id] = 0;
  }

  _Cilk_for(id=0;id<gn;id++){
    long long start = (long long)id*nd/gn;
    long long end = (long long)(id+1)*nd/gn;
    long long im;
    for(im=start;im<end;im++){
      if(zgmid[lod[im]][ncut]<=zlmid){
	lessNum[id]++;
      }else{
	moreNum[id]++;
      }
    }
  }

  for(id=0;id<gn;id++){
    nl += lessNum[id];
  }

  moreStart[0] = nl;
  if(nl != 0 && nl != nd){
    for(id=0;id<gn-1;id++){
      lessStart[id+1] = lessStart[id] + lessNum[id];
      moreStart[id+1] = moreStart[id] + moreNum[id];
    }
    free(moreNum);
    free(lessNum);
    
    _Cilk_for(id=0;id<gn;id++){
      long long start = (long long)id*nd/gn;
      long long end = (long long)(id+1)*nd/gn;
      long long im;
      for(im=start;im<end;im++){
	if(zgmid[lod[im]][ncut]<=zlmid){
	  templod[lessStart[id]] = lod[im];
	  lessStart[id]++;
	}else{
	  templod[moreStart[id]] = lod[im];
	  moreStart[id]++;
	}
      }
    }
    free(lessStart);
    free(moreStart);
  }
  double end = get_wall_time();
  double spent = end - start;
  printf("succession access parallelization time spent:%.10f\n",spent);
}

void random(int *lod,int nd){
  int i, j;
  int a,b;

  srand(time(NULL));

  for (i = 0; i < nd; i++) {
      a = rand() % nd;
      b = rand() % nd;
      int temp = lod[a];
      lod[a] = lod[b];
      lod[b] = temp;
  }
}

int contains(int *lod, int n, int x)
{
  while (n--) {
    if (lod[n] == x) return 1;
  }
  return 0;
}


void mm2(double (*zgmid)[3],int *lod,int nd,int md, int ndim){
  int il,id;
  double zlmin[ndim],zlmax[ndim];
  CILK_C_REDUCER_MAX(max_x,double,-10.1);
  CILK_C_REDUCER_MIN(min_x,double,10.1);
  CILK_C_REDUCER_MAX(max_y,double,-10.1);
  CILK_C_REDUCER_MIN(min_y,double,10.1);
  CILK_C_REDUCER_MAX(max_z,double,-10.1);
  CILK_C_REDUCER_MIN(min_z,double,10.1);
  CILK_C_REGISTER_REDUCER(max_x);
  CILK_C_REGISTER_REDUCER(min_x);
  CILK_C_REGISTER_REDUCER(max_y);
  CILK_C_REGISTER_REDUCER(min_y);
  CILK_C_REGISTER_REDUCER(max_z);
  CILK_C_REGISTER_REDUCER(min_z);
#pragma cilk grainsize = 128
  _Cilk_for(il=0;il<nd;il++){
    CILK_C_REDUCER_MIN_CALC(min_x,zgmid[lod[il]][0]);
    CILK_C_REDUCER_MAX_CALC(max_x,zgmid[lod[il]][0]);
    CILK_C_REDUCER_MIN_CALC(min_y,zgmid[lod[il]][1]);
    CILK_C_REDUCER_MAX_CALC(max_y,zgmid[lod[il]][1]);
    CILK_C_REDUCER_MIN_CALC(min_z,zgmid[lod[il]][2]);
    CILK_C_REDUCER_MAX_CALC(max_z,zgmid[lod[il]][2]);
  }
  zlmin[0] = REDUCER_VIEW(min_x);
  zlmax[0] = REDUCER_VIEW(max_x);
  zlmin[1] = REDUCER_VIEW(min_y);
  zlmax[1] = REDUCER_VIEW(max_y);
  zlmin[2] = REDUCER_VIEW(min_z);
  zlmax[2] = REDUCER_VIEW(max_z);
  CILK_C_UNREGISTER_REDUCER(max_x);
  CILK_C_UNREGISTER_REDUCER(min_x);
  CILK_C_UNREGISTER_REDUCER(max_y);
  CILK_C_UNREGISTER_REDUCER(min_y);
  CILK_C_UNREGISTER_REDUCER(max_z);
  CILK_C_UNREGISTER_REDUCER(min_z);
  for(id=0;id<ndim;id++){
    printf("%lf %lf\n",zlmin[id],zlmax[id]);
  }
}

/*****create a cluster*********/
cluster * create_cluster(int nmbr,int ndpth,int nstrt,int nsize,int ndim,int nson){
  cluster *st_clt;
  st_clt = (cluster *)malloc(sizeof(cluster));
  nmbr = nmbr + 1;
  st_clt->nstrt = nstrt;
  st_clt->nsize = nsize;
  st_clt->ndim = ndim;
  st_clt->nnson = nson;
  st_clt->nmbr = nmbr;
  st_clt->ndpth = ndpth;
  st_clt->pc_sons = (cluster **)malloc(nson * sizeof(cluster*));

  return st_clt;
}

/****create cluster tree******/
cluster * create_ctree_ssgeom(cluster *st_clt,   //the current node
                              double (*zgmid)[3],     //coordination of objects
                              double param[],    //param[21] is minGroup, param[31] is ? 
                              int *lod,          //hash table of pre-index and pro-index
                              int ndpth,         //depth of the tree
                              int ndscd,
                              int nsrt,          //the start index of list
                              int nd,            //the length of list
                              int md,            //number of data
                              int ndim,          //number of dimension (default is 3)
                              int nclst){        //number of cluster
  int id,il,nson,ic;
  double minsz = param[21];
  double zcoef = param[31];
  double zlmin[ndim],zlmax[ndim];
  
  ndpth = ndpth + 1;
  if(nd <= minsz){
    nson = 0;
    st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
  }else{
    if(ndpth < 8){
      CILK_C_REDUCER_MAX(max,double,-10.1);
      CILK_C_REDUCER_MIN(min,double,10.1);
      CILK_C_REGISTER_REDUCER(max);
      CILK_C_REGISTER_REDUCER(min);
      _Cilk_for(id=0;id<ndim;id++){
        REDUCER_VIEW(min) = 10.1;
        REDUCER_VIEW(max) = -10.1;
        _Cilk_for(il=0;il<nd;il++){
          CILK_C_REDUCER_MIN_CALC(min,zgmid[lod[il]][id]);
          CILK_C_REDUCER_MAX_CALC(max,zgmid[lod[il]][id]);
        }
        zlmin[id] = REDUCER_VIEW(min);
        zlmax[id] = REDUCER_VIEW(max);
      }
      CILK_C_UNREGISTER_REDUCER(max);
      CILK_C_UNREGISTER_REDUCER(min);
    }else{
      for(id=0;id<ndim;id++){
        zlmin[id] = zgmid[lod[0]][id];
        zlmax[id] = zlmin[id];
        for(il=1;il<nd;il++){
          double zg = zgmid[lod[il]][id];
          if(zg < zlmin[id]){
            zlmin[id] = zg;
          }else if(zlmax[id] < zg){
            zlmax[id] = zg;
          }
        }
      }
    }
    double zdiff = zlmax[0] - zlmin[0];
    int ncut = 0;
    for(id=0;id<ndim;id++){
      double zidiff = zlmax[id]-zlmin[id];
      if(zidiff > zcoef * zdiff){
        zdiff = zidiff;
        ncut = id;
      }
    }
    
    double zlmid = 0.5 * (zlmax[ncut] + zlmin[ncut]);
    int nl = 0;
    int nr = nd-1;
    if(ndpth < 8){
      int gs = 5000000;
      int gn = nd/gs + 1;
      int lessNum[gn];
      int moreNum[gn];
      int lessStart[gn];
      int moreStart[gn];
      for(id=0;id<gn;id++){
        lessNum[id] = 0;
        moreNum[id] = 0;
        lessStart[id] = 0;
        moreStart[id] = 0;
      }

      _Cilk_for(id=0;id<gn;id++){
        int im;
        for(im=nd/gn*id;im<nd/gn*(id+1);im++){
          if(zgmid[lod[im]][ncut]<=zlmid){
            lessNum[id]++;
          }else{
            moreNum[id]++;
          }
        }
      }

      for(id=0;id<gn;id++){
        nl += lessNum[id];
      }

      moreStart[0] = nl;
      if(nl != 0 && nl != nd){
        for(id=0;id<gn-1;id++){
          lessStart[id+1] = lessStart[id] + lessNum[id];
          moreStart[id+1] = moreStart[id] + moreNum[id];
        }

        int *templod;
        templod = (int *)malloc(nd*sizeof(int));        

        _Cilk_for(id=0;id<gn;id++){
          int im;
          for(im=nd/gn*id;im<nd/gn*(id+1);im++){
            if(zgmid[lod[im]][ncut]<=zlmid){
              templod[lessStart[id]] = lod[im];
              lessStart[id]++;
            }else{
              templod[moreStart[id]] = lod[im];
              moreStart[id]++;
            }
          }
        }
        cilk_for(id=0;id<nd;id++){
          lod[id] = templod[id];
        }
        free(templod);
      }
    }else{
      while(nl < nr){
        while(nl < nd && zgmid[lod[nl]][ncut] <= zlmid){
          nl = nl + 1;
        }
        while(nr >= 0 && zgmid[lod[nr]][ncut] > zlmid){
          nr = nr - 1;
        }
        if(nl < nr){
          int nh = lod[nl];
          lod[nl] = lod[nr];
          lod[nr] = nh;
        }
      }
    }
    if(nl == nd || nl == 0){
      nson = 1;
      st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
      st_clt->pc_sons[0] = create_ctree_ssgeom(st_clt->pc_sons[0],zgmid,param,lod,
                                               ndpth,ndscd,nsrt,nd,md,ndim,nclst);
    }else{
      nson = 2;
      st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
      int nsrt1 = nsrt;
      int nd1 = nl;

      if(ndpth<8){
        st_clt->pc_sons[0] = _Cilk_spawn create_ctree_ssgeom(st_clt->pc_sons[0],zgmid,param,lod,
                                                             ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
      }else{
        st_clt->pc_sons[0] = create_ctree_ssgeom(st_clt->pc_sons[0],zgmid,param,lod,
                                                 ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
      }

      nsrt1 = nsrt + nl;
      nd1 = nd - nl;
      st_clt->pc_sons[1] = create_ctree_ssgeom(st_clt->pc_sons[1],zgmid,param,&lod[nl],
                                               ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
      if(ndpth<8){
        _Cilk_sync;
      }
    }
  }
  st_clt->ndscd = nd;
  return st_clt;
}


double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){ 
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
  return (double)clock() / CLOCKS_PER_SEC;
}

