#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <pthread.h>

#define INPUT_DEFAULT "bem_data/input_50ms.txt"

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
int NLF;
void create_leafmtx(leafmtx **st_leafmtx,cluster *stcltl,cluster *st_cltt,double param[],int *lnmtx,int nffc);
double dist_2cluster(cluster *st_cltl,cluster *stcltt);
void count_lntmx(cluster *st_cltl,cluster *st_cltt,double param[],int *lnmtx,int nffc);
void cal_bndbox_cog(cluster *st_clt,double (*zgmid)[3],int *lod,int nofc);
void set_bndbox_cog(cluster *st_clt,double (*zgmid)[3],int *lod,int nofc);
cluster * create_cluster(int nmbr,int ndpth,int nstrt,int nsize,int ndim,int nson);
void free_st_clt(cluster *st_clt);
cluster * create_ctree_ssgeom(cluster *st_clt,double (*zgmid)[3],double param[],int *lod,int ndpth,int ndscd,int nsrt,int nd,int md,int ndim,int nclst);
double get_wall_time();
double get_cpu_time();

/****defination of cilk_lock****/
pthread_mutex_t mut;

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
  int i,nfl,nflkt,ip,il,ig;
  int nd = nofc * nffc;
  int *lodfc;
  leafmtx **st_leafmtx;
  lodfc = (int *)malloc(nofc*sizeof(int));
  for(il=0;il<nofc;il++){
    lodfc[il] = il;
  }
  int nsrt = 0;
  int ndf = nofc;
  int nclst = 0;
  int ndpth = 0;
  int ndscd = 0;

  double start,end,spent;

  /*******create cluster tree*******/
  start = get_wall_time();
  st_clt = create_ctree_ssgeom(st_clt,gmid,param,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst);
  end = get_wall_time();
  spent = end - start;
  printf("cluster tree time spent:%.10f\n",spent);

  /*****add info of boundary box for cluster*****/
  set_bndbox_cog(st_clt,gmid,lodfc,nofc);

  ndpth = 0;

  /*****count number of leafmatrices*****/
  count_lntmx(st_clt,st_clt,param,lnmtx,nffc);

  /*****create block cluster tree*****/
  st_leafmtxp->nlfkt = lnmtx[0];
  int nlf = lnmtx[0] + lnmtx[1];
  st_leafmtx = (leafmtx **)malloc(nlf*sizeof(leafmtx*));
  st_leafmtxp->nlf = nlf;
  printf("nlf:%d\n",nlf);
  
  NLF = 0;
  pthread_mutex_init(&mut,NULL);
  start = get_wall_time();
  create_leafmtx(st_leafmtx,st_clt,st_clt,param,lnmtx,nffc);
  end = get_wall_time();
  spent = end - start;
  printf("block cluster tree time spent:%.10f\n",spent/5);
  printf("NLF:%d\n",NLF);

  /*****sort the array of leafmatrices*****/
  qsort_row_leafmtx(st_leafmtx,0,nlf-1);
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
  }

  //for(i=0;i<st_leafmtxp->nlf;i++){
  //  printf("nstrtl:%d ndl:%d nstrtt:%d ndt:%d\n",st_leafmtx[i]->nstrtl,st_leafmtx[i]->ndl,st_leafmtx[i]->nstrtt,st_leafmtx[i]->ndt);
  //}

  free_st_clt(st_clt);

  for(il=0;il<nofc;il++){
    for(ig=0;ig<nffc;ig++){
      int is = ig + il * nffc;
      lod[is] = lodfc[il];
    }
  }
}

void qsort_row_leafmtx(leafmtx **st_leafmtx,int first,int last){
  int i,j,pivot;
  leafmtx *st_www;
  if(first<last){
    pivot=first;
    i=first;
    j=last;
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
    qsort_row_leafmtx(st_leafmtx,first,j-1);
    qsort_row_leafmtx(st_leafmtx,j+1,last);
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
    qsort_col_leafmtx(st_leafmtx,first,j-1);
    qsort_col_leafmtx(st_leafmtx,j+1,last);
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

void create_leafmtx(leafmtx **st_leafmtx,cluster *st_cltl,cluster *st_cltt,
                    double param[],int *lnmtx,int nffc){
  int ndl = st_cltl->nsize * nffc;
  int ndt = st_cltt->nsize * nffc;
  int nstrtl = st_cltl->nstrt;
  int nstrtt = st_cltt->nstrt;
  int nnsonl = st_cltl->nnson;
  int nnsont = st_cltt->nnson;
  int il,it;

  double nleaf = param[41];
  double zeta = param[51];
  double zdistlt = dist_2cluster(st_cltl,st_cltt);
  int condition1 = (st_cltl->zwdth * zeta <= zdistlt || st_cltt->zwdth * zeta <= zdistlt) && (ndl >= nleaf && ndt >= nleaf);
  int condition2 = nnsonl == 0 || nnsont == 0 || ndl <= nleaf || ndt <= nleaf;

  if(condition1 || condition2){
    int my_nlf;

    pthread_mutex_lock(&mut);
    my_nlf = NLF++;
    pthread_mutex_unlock(&mut);

    st_leafmtx[my_nlf] = (leafmtx *)malloc(sizeof(leafmtx));
    st_leafmtx[my_nlf]->nstrtl = nstrtl;
    st_leafmtx[my_nlf]->ndl = ndl;
    st_leafmtx[my_nlf]->nstrtt = nstrtt;
    st_leafmtx[my_nlf]->ndt = ndt;
    if(condition2){
      st_leafmtx[my_nlf]->ltmtx = 2;
    }else{
      st_leafmtx[my_nlf]->kt = 0;
      st_leafmtx[my_nlf]->ltmtx = 1;
    }
  }else{
    for(il=0;il<nnsonl;il++){
      for(it=0;it<nnsont;it++){
	if(il == nnsonl-1 && it == nnsont-1){
	  create_leafmtx(st_leafmtx,st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lnmtx,nffc);
	}else{
	  cilk_spawn create_leafmtx(st_leafmtx,st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lnmtx,nffc);
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
      zs = zs+ (st_cltl->bmin[id] - st_cltt->bmax[id]) * (st_cltl->bmin[id] - st_cltt->bmax[id]);
    }
  }
  return sqrt(zs);
}

void count_lntmx(cluster *st_cltl,cluster *st_cltt,double param[],int *lnmtx,int nffc){
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
  if ((st_cltl->zwdth <= zeta * zdistlt || st_cltt->zwdth <= zeta * zdistlt) && (ndl >= nleaf && ndt >= nleaf)){
    lnmtx[0] = lnmtx[0] + 1;
  }else{
    if(nnsonl == 0 || nnsont == 0 || ndl <= nleaf || ndt <= nleaf){
      lnmtx[1] = lnmtx[1] + 1;
    }else{
      lnmtx[2] = lnmtx[2] + 1;
      for(il=0;il<nnsonl;il++){
	for(it=0;it<nnsont;it++){
	  count_lntmx(st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lnmtx,nffc);
	}
      }
    }
  }
}

void cal_bndbox_cog(cluster *st_clt,double (*zgmid)[3],int *lod,int nofc){
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
      st_clt->bmin[id] = zgmid[lod[0]][id];
      st_clt->bmax[id] = zgmid[lod[0]][id];
    }
    for(id=0;id<ndim;id++){
      for(il=1;il<st_clt->nsize;il++){
	if(zgmid[lod[il]][id] < st_clt->bmin[id]){
	  st_clt->bmin[id] = zgmid[lod[il]][id];
	}
	if(st_clt->bmax[id] < zgmid[lod[il]][id]){
	  st_clt->bmax[id] = zgmid[lod[il]][id];
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

void set_bndbox_cog(cluster *st_clt, double (*zgmid)[3], int *lod, int nofc){
  int ic,l;

  for(ic=0;ic<st_clt->nnson;ic++){
    if(ic == 0){
      l = 0;
    }else{
      l = l + st_clt->pc_sons[ic-1]->nsize;
    }
    set_bndbox_cog(st_clt->pc_sons[ic],zgmid,&lod[l],nofc);
  }
  cal_bndbox_cog(st_clt,zgmid,lod,nofc);
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
			      double param[],    //param[21] is minGroup, param[31] is ? 
			      int *lod,          //hash table of pre-index and pro-index
			      int ndpth,         //depth of the tree
			      int ndscd,
			      int nsrt,          //the start index of list
			      int nd,            //the length of list
			      int md,            //number of data
			      int ndim,          //number of dimension (default is 3)
			      int nclst){        //number of cluster
  int id,il,nson;
  double minsz = param[21];
  double zcoef = param[31];
  double zlmin[ndim],zlmax[ndim];
  ndpth = ndpth + 1;
  if(nd <= minsz){
    nson = 0;
    st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
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
      st_clt->pc_sons[0] = cilk_spawn create_ctree_ssgeom(st_clt->pc_sons[0],zgmid,param,lod,
					       ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);

      nsrt1 = nsrt + nl;
      nd1 = nd - nl;
      st_clt->pc_sons[1] = create_ctree_ssgeom(st_clt->pc_sons[1],zgmid,param,&lod[nl],
					       ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
      cilk_sync;
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
