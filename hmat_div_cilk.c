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

#define INPUT_DEFAULT "bem_data/input_50ms.txt"
#define PARA_LEVEL 11
#define MM_LEVEL 8
#define SPAWN_LEVEL 11
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
void qsort_col_leafmtx(leafmtx **st_leafmtx,int first,int last);
void qsort_row_leafmtx(leafmtx **st_leafmtx,int first,int last);
int med3(int nl,int nr,int nlr2);
void create_leafmtx(cluster *st_cltl,cluster *st_cltt);
double dist_2cluster(cluster *st_cltl,cluster *stcltt);
void count_lntmx2(cluster *st_cltl,cluster *st_cltt,double param[],int **ln,int nffc);
void count_lntmx(cluster *st_cltl,cluster *st_cltt,double param[],int *lntmx,int nffc);
void cal_bndbox_cog(cluster *st_clt,double (*zgmid)[3],int nofc);
void set_bndbox_cog(cluster *st_clt,double (*zgmid)[3],int nofc);
cluster * create_cluster(int nmbr,int ndpth,int nstrt,int nsize,int ndim,int nson);
void free_st_clt(cluster *st_clt);
cluster * create_ctree_ssgeom(cluster *st_clt,double (*zgmid)[3],double (*tempzgmid)[3],double param[],int ndpth,int ndscd,int nsrt,int nd,int md,int ndim,int nclst);
double get_wall_time();
double get_cpu_time();
int depth_max(cluster *st_clt);

void checkClusterTree(FILE *f,cluster *st_clt);

int countlist[600];
int depth_m = 0;
leafmtx **temp_leafmtx;
int **conn;
double param[100];
int nffc;

int main(int argc, char **argv){
  if(argc >= 2){
    __cilkrts_set_param("nworkers",argv[1]);
  }
  int nworkers = __cilkrts_get_nworkers();
  printf("number of workers:%d\n",nworkers);
  /******** read file *********/
  char *fname;
  FILE *file;
  int countOfNode=0;
  int count = 0;
  int i;
  double (*coordOfNode)[3];
  double (*coordOfFace)[3];
  fname = (argc >= 3)?argv[2]:INPUT_DEFAULT;
  file = fopen(fname,"r");
  if(file == NULL){
    printf("Error: Unable to input file '%s'!\n", fname);
    exit (99);
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
  nffc = 1;
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
  int *lodfc;
  int lel;
  leafmtx *st_leafmtx;
  lodfc = (int *)malloc(nofc*sizeof(int));
  for(il=0;il<nofc;il++){
    lodfc[il] = il;
  }
  double (*tempgmid)[3];
  tempgmid = (double(*)[3])malloc(nofc*3*sizeof(double));
  int nsrt = 0;
  int ndf = nofc;
  int nclst = 0;
  int ndpth = 0;
  int ndscd = 0;
  int nworkers = __cilkrts_get_nworkers();
  double start,end,spent,start1,end1,spent1,start2,end2,spent2;

  /*******create cluster tree*******/
  start = get_wall_time();
  st_clt = create_ctree_ssgeom(st_clt,gmid,tempgmid,param,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst);
  end = get_wall_time();
  spent = end - start;
  printf("cluster tree time spent:%.10f\n",spent);

  /*****add info of boundary box for cluster*****/
  set_bndbox_cog(st_clt,gmid,nofc);

  ndpth = 0;
  free(tempgmid);
  free(gmid);
  /*****count number of leafmatrices*****/
  start1 = get_wall_time();
  count_lntmx(st_clt,st_clt,param,lnmtx,nffc);
  end1 = get_wall_time();
  spent1 = end1 - start1;
  printf("count time:%.10f\n",spent1);
  int nlf;
  
  /*****create block cluster tree*****/
  st_leafmtxp->nlfkt = lnmtx[0];
  nlf = lnmtx[0] + lnmtx[1];
  st_leafmtx = (leafmtx *)malloc(nlf*sizeof(leafmtx));
  st_leafmtxp->nlf = nlf;
  printf("nlf:%d\n",nlf);
  
  //int *countlist;
  //countlist = (int *)malloc(nworkers*sizeof(int));
  for(i=0;i<600;i++){
    countlist[i] = 0;
  }

  //lel = nlf/nworkers*3;
  temp_leafmtx = (leafmtx **)malloc(nworkers * sizeof(leafmtx*));
  conn = (int **)malloc(nworkers * sizeof(int*));
  for(i=0;i<nworkers;i++){
    temp_leafmtx[i] = (leafmtx *)malloc(nlf * sizeof(leafmtx));
    conn[i] = (int *)malloc(sizeof(int));
    conn[i][0] = 0;
  }
    
  int sum = 0;
  start2 = get_wall_time();
  create_leafmtx(st_clt,st_clt);
  end2 = get_wall_time();
  spent2 = end2 - start2;
  printf("block cluster tree time spent:%.10f\n",spent2);  
  //long sumk=0;
  /*for(i=0;i<nworkers;i++){
    long st = (long)i * (long)lel;
    for(j=0;j<countlist[i];j++){
      int a = sum + j;
      long b = st + (long)j;
      st_leafmtx[a] = &temp_leafmtx[b];
    }
    sum += countlist[i];
    }*/
  for(i=0;i<nworkers;i++){
    for(j=0;j<conn[i][0];j++){
      st_leafmtx[sum+j] = temp_leafmtx[i][j];
    }
    sum += conn[i][0];
  }

  //printf("block cluster tree time spent:%.10f\n",spent2);
  printf("sum:%d\n",sum);

  for(i=0;i<nworkers;i++){
    printf("%ld,",countlist[i]/*conn[i][0]*/);
  }
  printf("\n");
  /*printf("\n");
  for(i=0;i<nworkers;i++){
    printf("%d,",ncall[i]);
  }
  printf("\n");
  printf("\n");
  for(i=0;i<nworkers;i++){
    printf("judge time:%f,node time:%f\n",ntime[0][i],ntime[1][i]);
    }*/
  //free(temp_leafmtx);

  /*****sort the array of leafmatrices*****/

  /*qsort_row_leafmtx(st_leafmtx,0,nlf-1);
  int ilp = 0;
  int ips = 0;
  for(ip=0;ip<nlf;ip++){
    il = st_leafmtx[ip]->nstrtl;
    if(il < ilp){
      printf("Error!: supermatrix_construction_cog/leafmtx row_sort\n");
    }else if(il > ilp){
      qsort_col_leafmtx(st_leafmtx,ips,ip-1);
      ilp = il;
      ips = ip;
    }
    if(ip == nlf-1){
      qsort_col_leafmtx(st_leafmtx,ips,nlf-1);
    }
    }*/
  /*FILE *f = fopen("result.txt","w");
  for(i=0;i<st_leafmtxp->nlf;i++){
    //printf("nstrtl:%d ndl:%d nstrtt:%d ndt:%d\n",st_leafmtx[i]->nstrtl,st_leafmtx[i]->ndl,st_leafmtx[i]->nstrtt,st_leafmtx[i]->ndt);
    fprintf(f,"%d %d %d %d\n",st_leafmtx[i]->nstrtl,st_leafmtx[i]->ndl,st_leafmtx[i]->nstrtt,st_leafmtx[i]->ndt);
  }
  fclose(f);
  free_st_clt(st_clt);*/

  /*for(il=0;il<nofc;il++){
    for(ig=0;ig<nffc;ig++){
      int is = ig + il * nffc;
      lod[is] = lodfc[il];
    }
    }*/
}

void qsort_row_leafmtx(leafmtx **st_leafmtx,int first,int last){
  int i,j,pivot;
  leafmtx *st_www;
  if(first<last){
    pivot=(first+last)/2;
    i=first;
    j=last;
    st_www = st_leafmtx[pivot];
    st_leafmtx[pivot] = st_leafmtx[1];
    st_leafmtx[1] = st_www;
    pivot = first;
    while(i<j){
      while(st_leafmtx[i]->nstrtl <= st_leafmtx[pivot]->nstrtl && i<last)
        i++;
      while(st_leafmtx[j]->nstrtl>st_leafmtx[pivot]->nstrtl)
        j--;
      if(i<j){
        st_www = st_leafmtx[i];
        st_leafmtx[i]=st_leafmtx[j];
        st_leafmtx[j]=st_www;
      }
    }
    st_www = st_leafmtx[pivot];
    st_leafmtx[pivot] = st_leafmtx[j];
    st_leafmtx[j] = st_www;
    _Cilk_spawn qsort_row_leafmtx(st_leafmtx,first,j-1);
    qsort_row_leafmtx(st_leafmtx,j+1,last);
    _Cilk_sync;
  }
}

void qsort_col_leafmtx(leafmtx **st_leafmtx,int first,int last){
  int i,j,pivot;
  leafmtx *st_www;
  if(first<last){
    pivot=first;
    i=first;
    j=last;
    while(i<j){
      while(st_leafmtx[i]->nstrtt <= st_leafmtx[pivot]->nstrtt && i<last)
	i++;
      while(st_leafmtx[j]->nstrtt>st_leafmtx[pivot]->nstrtt)
	j--;
      if(i<j){
	st_www = st_leafmtx[i];
	st_leafmtx[i]=st_leafmtx[j];
	st_leafmtx[j]=st_www;
      }
    }
    st_www = st_leafmtx[pivot];
    st_leafmtx[pivot] = st_leafmtx[j];
    st_leafmtx[j] = st_www;
    _Cilk_spawn qsort_col_leafmtx(st_leafmtx,first,j-1);
    qsort_col_leafmtx(st_leafmtx,j+1,last);
    _Cilk_sync;
  }
}

int med3(int nl,int nr,int nlr2){
  int med3;
  if(nl < nr){
    if(nr < nlr2){
      med3 = nr;
    }else if(nlr2 < nl){
      med3 = nl;
    }else{
      med3 = nlr2;
    }
  }else{
    if(nlr2 < nr){
      med3 = nr;
    }else if(nl < nlr2){
      med3 = nl;
    }else{
      med3 = nlr2;
    }
  }
  return med3;
}

void create_leafmtx(cluster *st_cltl,cluster *st_cltt){
  //int ndl = st_cltl->nsize * nffc;
  //int ndt = st_cltt->nsize * nffc;
  //int nstrtl = st_cltl->nstrt;
  //int nstrtt = st_cltt->nstrt;
  //int nnsonl = st_cltl->nnson;
  //int nnsont = st_cltt->nnson;
  //int st_cltl_depth = st_cltl->ndpth;
  //int st_cltt_depth = st_cltt->ndpth;
  //int il;
  //int my_num = __cilkrts_get_worker_number();

  double nleaf = param[41];
  double zeta = param[51];
  double zdistlt = dist_2cluster(st_cltl,st_cltt);

  if((st_cltl->zwdth * zeta <= zdistlt || st_cltt->zwdth * zeta <= zdistlt) && (st_cltl->nsize * nffc >= nleaf && st_cltt->nsize * nffc >= nleaf)){
    int my_num = __cilkrts_get_worker_number();    
    int lnlf = countlist[my_num];//conn[my_num][0];
    //leafmtx (* restrict a)[nlf] = (leafmtx(*)[nlf])&temp_leafmtx[0];
    temp_leafmtx[my_num][lnlf].nstrtl = st_cltl->nstrt;
    temp_leafmtx[my_num][lnlf].ndl = st_cltl->nsize * nffc;
    temp_leafmtx[my_num][lnlf].nstrtt = st_cltt->nstrt;
    temp_leafmtx[my_num][lnlf].ndt = st_cltt->nsize * nffc;
    temp_leafmtx[my_num][lnlf].kt = 0;
    temp_leafmtx[my_num][lnlf].ltmtx = 1;

    //conn[my_num][0]++;
    countlist[my_num]++;
  }else if(st_cltl->nnson == 0 || st_cltt->nnson == 0 || st_cltl->nsize * nffc <= nleaf || st_cltt->nsize * nffc <= nleaf){
    int my_num = __cilkrts_get_worker_number();
    int lnlf = countlist[my_num];//conn[my_num][0];
    //leafmtx (* restrict a)[nlf] = (leafmtx(*)[nlf])&temp_leafmtx[0];
    temp_leafmtx[my_num][lnlf].nstrtl = st_cltl->nstrt;
    temp_leafmtx[my_num][lnlf].ndl = st_cltl->nsize * nffc;
    temp_leafmtx[my_num][lnlf].nstrtt = st_cltt->nstrt;
    temp_leafmtx[my_num][lnlf].ndt = st_cltt->nsize * nffc;
    temp_leafmtx[my_num][lnlf].kt = 0;
    temp_leafmtx[my_num][lnlf].ltmtx = 2;

    //conn[my_num][0]++;
    countlist[my_num];
  }else{
    if(st_cltl->ndpth < 16 && st_cltt->ndpth < 16){
      int il;
      cilk_for(il=0;il<st_cltl->nnson;il++){
	int it;
        for(it=0;it<st_cltt->nnson;it++){
          create_leafmtx(st_cltl->pc_sons[il],st_cltt->pc_sons[it]);
        }
      }
    }else{
      int il;
      for(il=0;il<st_cltl->nnson;il++){
	int it;
        for(it=0;it<st_cltt->nnson;it++){
          create_leafmtx(st_cltl->pc_sons[il],st_cltt->pc_sons[it]);
        }
      }
    }
  }
}

double dist_2cluster(cluster *st_cltl,cluster *st_cltt){
  double zs = 0.0;
  int id;
  for(id=0;id<st_cltl->ndim;id++){
    if(st_cltl->bmax[id] < st_cltt->bmin[id]){
      zs = zs + (st_cltt->bmin[id] - st_cltl->bmax[id]) * (st_cltt->bmin[id] - st_cltl->bmax[id]);
    }else if(st_cltt->bmax[id]< st_cltl->bmin[id]){
      zs = zs + (st_cltl->bmin[id] - st_cltt->bmax[id]) * (st_cltl->bmin[id] - st_cltt->bmax[id]);
    }
  }
  return sqrt(zs);
}

void count_lntmx2(cluster *st_cltl,cluster *st_cltt,double param[],int **ln,int nffc){
  int il,it;
  int ndl = st_cltl->nsize * nffc;
  int ndt = st_cltt->nsize * nffc;
  int nstrtl = st_cltl->nstrt;
  int nstrtt = st_cltt->nstrt;
  int nnsonl = st_cltl->nnson;
  int nnsont = st_cltt->nnson;
  int my_num = __cilkrts_get_worker_number();

  double nleaf = param[41];
  double zeta = param[51];
  double zdistlt = dist_2cluster(st_cltl,st_cltt);
  if ((st_cltl->zwdth * zeta <= zdistlt || st_cltt->zwdth * zeta <= zdistlt) && (ndl >= nleaf && ndt >= nleaf)){
    ln[my_num][0] = ln[my_num][0] + 1; 
  }else{
    if(nnsonl == 0 || nnsont == 0 || ndl <= nleaf || ndt <= nleaf){
      ln[my_num][1] = ln[my_num][1] + 1;
    }else{
      ln[my_num][2] = ln[my_num][2] + 1;
      _Cilk_for(il=0;il<nnsonl;il++){
	_Cilk_for(it=0;it<nnsont;it++){
	  count_lntmx2(st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,ln,nffc);
	}
      }
    }
  }
}

void count_lntmx(cluster *st_cltl,cluster *st_cltt,double param[],int *lntmx,int nffc){
  int il,it;
  int ndl = st_cltl->nsize * nffc;
  int ndt = st_cltt->nsize * nffc;
  int nstrtl = st_cltl->nstrt;
  int nstrtt = st_cltt->nstrt;
  int nnsonl = st_cltl->nnson;
  int nnsont = st_cltt->nnson;

  double nleaf = param[41];
  double zeta = param[51];
  double zdistlt = dist_2cluster(st_cltl,st_cltt);
  if ((st_cltl->zwdth * zeta <= zdistlt || st_cltt->zwdth * zeta <= zdistlt) && (ndl >= nleaf && ndt >= nleaf)){
    lntmx[0] = lntmx[0] + 1;
  }else{
    if(nnsonl == 0 || nnsont == 0 || ndl <= nleaf || ndt <= nleaf){
      lntmx[1] = lntmx[1] + 1;
    }else{
      lntmx[2] = lntmx[2] + 1;
      for(il=0;il<nnsonl;il++){
        for(it=0;it<nnsont;it++){
          count_lntmx(st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lntmx,nffc);
        }
      }
    }
  }
}

void cal_bndbox_cog(cluster *st_clt,double (*zgmid)[3],int nofc){
  int ndim = st_clt->ndim;
  int id,il;
  st_clt->bmin = (double *)malloc(3*sizeof(double));
  st_clt->bmax = (double *)malloc(3*sizeof(double));
  double zeps = 1.0e-5;
  if(st_clt->nnson > 0){
    for(id=0;id<ndim;id++){
      st_clt->bmin[id] = st_clt->pc_sons[0]->bmin[id];
      st_clt->bmax[id] = st_clt->pc_sons[0]->bmax[id];
    }
    for(il=1;il<st_clt->nnson;il++){
      for(id=0;id<ndim;id++){
	if(st_clt->pc_sons[il]->bmin[id] < st_clt->bmin[id]){
	  st_clt->bmin[id] = st_clt->pc_sons[il]->bmax[id];
	}
	if(st_clt->bmax[id] < st_clt->pc_sons[il]->bmax[id]){
	  st_clt->bmax[id] = st_clt->pc_sons[il]->bmax[id];
	}
      }
    }
  }else{
    for(id=0;id<ndim;id++){
      st_clt->bmin[id] = zgmid[0][id];
      st_clt->bmax[id] = zgmid[0][id];
    }
    for(il=1;il<st_clt->nsize;il++){
      for(id=0;id<ndim;id++){
	if(zgmid[il][id] < st_clt->bmin[id]){
	  st_clt->bmin[id] = zgmid[il][id];
	}
	if(st_clt->bmax[id] < zgmid[il][id]){
	  st_clt->bmax[id] = zgmid[il][id];
	}
      }
    }
  }
  double zwdth = (st_clt->bmax[0] - st_clt->bmin[0]) * (st_clt->bmax[0] - st_clt->bmin[0]);
  for(id=1;id<ndim;id++){
    zwdth = zwdth + (st_clt->bmax[id] - st_clt->bmin[id]) * (st_clt->bmax[id] - st_clt->bmin[id]);
  }
  zwdth = sqrt(zwdth);
  for(id=0;id<ndim;id++){
    double bdiff = st_clt->bmax[id] - st_clt->bmin[id];
    if(bdiff < zeps * zwdth){
      st_clt->bmax[id] = st_clt->bmax[id] + 0.5 * (zeps * zwdth - bdiff);
      st_clt->bmin[id] = st_clt->bmin[id] - 0.5 * (zeps * zwdth - bdiff);
    }
  }
  zwdth = (st_clt->bmax[0] - st_clt->bmin[0]) * (st_clt->bmax[0] - st_clt->bmin[0]);
  for(id=1;id<ndim;id++){
    zwdth = (st_clt->bmax[id] - st_clt->bmin[id]) * (st_clt->bmax[id] - st_clt->bmin[id]);
  }
  st_clt->zwdth = sqrt(zwdth);
}

void set_bndbox_cog(cluster *st_clt, double (*zgmid)[3], int nofc){
  int ic;

  //  int l[st_clt->nnson];
  int l;
  /*for(ic=0;ic<st_clt->nnson;ic++){
    if(ic == 0){
      l[ic] = 0;
    }else{
      l[ic] = l[ic-1] + st_clt->pc_sons[ic-1]->nsize;
    }
    }*/
  for(ic=0;ic<st_clt->nnson;ic++){
    if(ic == 0){
      l = 0;
    }else{
      l = l + st_clt->pc_sons[ic-1]->nsize;
    }
    //set_bndbox_cog(st_clt->pc_sons[ic],&zgmid[l[ic]],nofc);
    set_bndbox_cog(st_clt->pc_sons[ic],&zgmid[l],nofc);
  }
  cal_bndbox_cog(st_clt,zgmid,nofc);
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

void checkClusterTree(FILE *f,cluster *st_clt){
  if(st_clt->ndpth<PARA_LEVEL){
    fprintf(f,"%d %d %d %d %lf\n",st_clt->nstrt,st_clt->nsize,st_clt->ndpth,st_clt->nnson,st_clt->zwdth);
  }
  if(st_clt->nnson==0){
    return;
  }else if(st_clt->nnson==1){
    checkClusterTree(f,st_clt->pc_sons[0]);
  }else{
    checkClusterTree(f,st_clt->pc_sons[0]);
    checkClusterTree(f,st_clt->pc_sons[1]);
  }
}

int depth_max(cluster *st_clt){
  int i;
  if(st_clt->ndpth > depth_m){
    depth_m = st_clt->ndpth;
  }
  for(i=0;i<st_clt->nnson;i++){
    depth_max(st_clt->pc_sons[i]);
  }
}

void free_st_clt(cluster *st_clt){
  int ic;
  int nnson = st_clt->nnson;
  for(ic=0;ic<nnson;ic++){
    free_st_clt(st_clt->pc_sons[ic]);
  }
  free(st_clt->bmin);
  free(st_clt->bmax);
  free(st_clt->pc_sons);
}


/****create cluster tree******/
cluster * create_ctree_ssgeom(cluster *st_clt,   //the current node
			      double (*zgmid)[3],     //coordination of objects
			      double (*tempzgmid)[3],
			      double param[],    //param[21] is minGroup, param[31] is ? 
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
    if(ndpth<PARA_LEVEL && ndpth%2==0){
#pragma simd
      _Cilk_for(id=0;id<nd;id++){
	for(il=0;il<ndim;il++){
	  zgmid[id][il] = tempzgmid[id][il];
	}
      }
    }
    st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
  }else{
    if(ndpth < MM_LEVEL){
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
	CILK_C_REDUCER_MIN_CALC(min_x,zgmid[il][0]);
	CILK_C_REDUCER_MAX_CALC(max_x,zgmid[il][0]);
	CILK_C_REDUCER_MIN_CALC(min_y,zgmid[il][1]);
	CILK_C_REDUCER_MAX_CALC(max_y,zgmid[il][1]);
	CILK_C_REDUCER_MIN_CALC(min_z,zgmid[il][2]);
	CILK_C_REDUCER_MAX_CALC(max_z,zgmid[il][2]);
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
    }else{
      for(id=0;id<ndim;id++){
	zlmin[id] = zgmid[0][id];
	zlmax[id] = zlmin[id];
	for(il=1;il<nd;il++){
          double zg = zgmid[il][id];
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
    if(ndpth < PARA_LEVEL){
      int gn = nd/CHUNK_SIZE;
      int *lessNum = (int *)malloc(gn*sizeof(int));
      int *moreNum = (int *)malloc(gn*sizeof(int));
      int *lessStart = (int *)malloc(gn*sizeof(int));
      int *moreStart = (int *)malloc(gn*sizeof(int));
      //#pragma cilk grainsize = 10000
      _Cilk_for(id=0;id<gn;id++){
	lessNum[id] = 0;
	moreNum[id] = 0;
	lessStart[id] = 0;
	moreStart[id] = 0;
      }

      _Cilk_for(id=0;id<gn;id++){
	long start = (long)id*nd/gn;
	long end = (long)(id+1)*nd/gn;
	long im;
	for(im=start;im<end;im++){
	  if(zgmid[im][ncut]<=zlmid){
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
	  long start = (long)id*nd/gn;
	  long end = (long)(id+1)*nd/gn;
	  long im;
	  for(im=start;im<end;im++){
	    if(zgmid[im][ncut]<=zlmid){
	      tempzgmid[lessStart[id]][0] = zgmid[im][0];
	      tempzgmid[lessStart[id]][1] = zgmid[im][1];
	      tempzgmid[lessStart[id]][2] = zgmid[im][2];
	      lessStart[id]++;
	    }else{
	      tempzgmid[moreStart[id]][0] = zgmid[im][0];
	      tempzgmid[moreStart[id]][1] = zgmid[im][1];
	      tempzgmid[moreStart[id]][2] = zgmid[im][2];
	      moreStart[id]++;
	    }
	  }
	}
	free(lessStart);
	free(moreStart);
      }
    }else{
      while(nl < nr){
	while(nl < nd && zgmid[nl][ncut] <= zlmid){
	  nl = nl + 1;
	}
	while(nr >= 0 && zgmid[nr][ncut] > zlmid){
	  nr = nr - 1;
	}
	if(nl < nr){
	  for(id=0;id<ndim;id++){
	    double nh = zgmid[nl][id];
	    zgmid[nl][id] = zgmid[nr][id];
	    zgmid[nr][id] = nh;
	  }
	}
      }
    }
    if(nl == nd || nl == 0){
      nson = 1;
      st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
      if(ndpth<PARA_LEVEL){
	st_clt->pc_sons[0] = create_ctree_ssgeom(st_clt->pc_sons[0],tempzgmid,zgmid,param,
					       ndpth,ndscd,nsrt,nd,md,ndim,nclst);
      }else{
	st_clt->pc_sons[0] = create_ctree_ssgeom(st_clt->pc_sons[0],zgmid,tempzgmid,param,
						 ndpth,ndscd,nsrt,nd,md,ndim,nclst);
      }
    }else{
      if(ndpth<PARA_LEVEL){
	nson = 2;
	st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
	int nsrt1 = nsrt;
	int nd1 = nl;
	st_clt->pc_sons[0] = _Cilk_spawn create_ctree_ssgeom(st_clt->pc_sons[0],tempzgmid,zgmid,param,
							     ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);

	nsrt1 = nsrt + nl;
	nd1 = nd - nl;
	st_clt->pc_sons[1] = create_ctree_ssgeom(st_clt->pc_sons[1],&tempzgmid[nl],&zgmid[nl],param,
						 ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
	_Cilk_sync;
      }else if(ndpth<SPAWN_LEVEL){
	nson = 2;
        st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
        int nsrt1 = nsrt;
        int nd1 = nl;
	st_clt->pc_sons[0] = _Cilk_spawn create_ctree_ssgeom(st_clt->pc_sons[0],zgmid,tempzgmid,param,
						 ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
	nsrt1 = nsrt + nl;
	nd1 = nd - nl;
	st_clt->pc_sons[1] = create_ctree_ssgeom(st_clt->pc_sons[1],&zgmid[nl],&tempzgmid[nl],param,
						 ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);	
	_Cilk_sync;
      }else{
	nson = 2;
        st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
        int nsrt1 = nsrt;
        int nd1 = nl;
        st_clt->pc_sons[0] = create_ctree_ssgeom(st_clt->pc_sons[0],zgmid,tempzgmid,param,
							     ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
        nsrt1 = nsrt + nl;
        nd1 = nd - nl;
        st_clt->pc_sons[1] = create_ctree_ssgeom(st_clt->pc_sons[1],&zgmid[nl],&tempzgmid[nl],param,
                                                 ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
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
