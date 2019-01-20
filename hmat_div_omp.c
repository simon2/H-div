#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include "data/bem_file.h"

#include <omp.h>

#define INPUT_DEFAULT "bem_data/input_50ms.txt"
#define OMPS 9
#define OMPL 512
#define PN 10000

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

typedef struct ompnode ompnode;
struct ompnode{
  cluster *st_clt;
  double (*zgmid)[3];
  double (*tempgmid)[3];
  int ndpth;
  int ndscd;
  int nsrt;
  int nd;
  int md;
  int ndim;
  int nclst;
  int nl;
};

void supermatrix_construction_cog_leafmtrx(leafmtxp *st_leafmtxp,double (*gmid)[3],double param[],int *lod,int *lnmtx,int nofc,int nffc,int ndim);
void qsort_col_leafmtx(leafmtx *st_leafmtx,int first,int last);
void qsort_row_leafmtx(leafmtx *st_leafmtx,int first,int last);
int med3(int nl,int nr,int nlr2);
void create_leafmtx(leafmtx *st_leafmtx,cluster *st_cltl,cluster *st_cltt,double param[],int *lnmtx,int nffc,int *nlf);
double dist_2cluster(cluster *st_cltl,cluster *st_cltt);
void count_lntmx(cluster *st_cltl,cluster *st_cltt,double param[],int *lnmtx,int nffc);
void cal_bndbox_cog(cluster *st_clt,double (*zgmid)[3],int *lod,int nofc);
void set_bndbox_cog(cluster *st_clt,double (*zgmid)[3],int *lod,int nofc);
cluster * create_cluster(int nmbr,int ndpth,int nstrt,int nsize,int ndim,int nson);
void free_st_clt(cluster *st_clt);
ompnode createOMPNode(cluster *st_clt,double (*zgmid)[3],double (*tempgmid)[3],int ndpth,int ndscd,int nsrt,int nd,int md,int ndim,int nclst,int nl);
cluster * create_ctree_ssgeom_t(cluster *st_clt,double (*zgmid)[3],double (*tempgmid)[3],double param[],int ndpth,int ndscd,int nsrt,int nd,int md,int ndim,int nclst);
void create_ctree_ssgeom_b(double param[]);
double get_wall_time();
double get_cpu_time();

void checkClusterTree(cluster *st_clt);

int depth_max;
int count_node;
int countOMP;
int nodecount = 0;
ompnode nodelist[OMPL];

int main(int argc, char **argv){
  if(argc >= 2){
    omp_set_num_threads(atoi(argv[1]));
  }
  printf("num of threads:%ld\n",omp_get_max_threads());
  /******** read file *********/
  char *fname;
  FILE *file;
  int countOfNode=0;
  int count = 0;
  int i;
  double (*coordOfNode)[3];
  double (*coordOfFace)[3];
  struct bem_input bi;
  fname = (argc >= 3)?argv[2]:INPUT_DEFAULT;
  file = fopen(fname,"r");
  if(file == NULL){
    printf("Error: Unable to input file '%s'!\n", fname);
    exit (99);
  }else{
    if (read_bem_input (file, &bi, BI_AUTO) == -1) {
      fprintf (stderr, "Bem input file read error!\n");
      exit (99);
    }
    char line[100];
    countOfNode = bi.nNode;
    coordOfNode = bi.coordOfNode;
    count = bi.nFace;
    coordOfFace = bi.coordOfFace;
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

  //double start,end;
  //start = get_wall_time();
  supermatrix_construction_cog_leafmtrx(st_leafmtxp,coordOfFace,param,lod,lnmtx,nofc,nffc,ndim);  // construction of Leaf-matrix 
  //end = get_wall_time();
  //double spent = end - start;
  //printf("total time spent:%.13f\n",spent);

  return 0;
}

#ifdef DEBUG
int *lod0;
#endif

void supermatrix_construction_cog_leafmtrx(leafmtxp *st_leafmtxp,    //the H-matrix
					   double (*gmid)[3],            //coordination of objects
					   double param[],int *lod,
					   int *lnmtx,              //1:k-rank 2:dense 3:H-matrix
					   int nofc,int nffc,       //number of elements in same coordination
					   int ndim){
  cluster *st_clt = (cluster *)malloc(sizeof(cluster));
  int i,nfl,nflkt,ip,il,ig;
  int nd = nofc * nffc;
  int *lodfc;
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

  depth_max = 0;
  count_node = 0;

  double start,end,spent;
  start = get_wall_time();
#ifdef DEBUG
  lod0 = lodfc;
#endif
  st_clt = create_ctree_ssgeom_t(st_clt,gmid,tempgmid,param,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst);
  create_ctree_ssgeom_b(param);
  end = get_wall_time();
  spent = end - start;
  printf("cluster tree time spent:%.10f\n",spent);
  printf("nodecount:%ld\n",nodecount);
  //printf("%ld\n",*(nodelist[0].st_clt)->ndpth);
  //checkClusterTree(st_clt);
  set_bndbox_cog(st_clt,gmid,lodfc,nofc);
  
  /*FILE *f = fopen("8_sequential.txt", "w");
  checkClusterTree(f,st_clt);
  fclose(f);*/

  ndpth = 0;
  start = get_wall_time();
  count_lntmx(st_clt,st_clt,param,lnmtx,nffc);
  end = get_wall_time();
  spent = end - start;
  printf("count time:%.10f\n",spent);

  st_leafmtxp->nlfkt = lnmtx[0];
  int nlf = lnmtx[0] + lnmtx[1];
  st_leafmtx = (leafmtx *)malloc(nlf*sizeof(leafmtx));
  st_leafmtxp->nlf = nlf;
  printf("nlf:%d\n",nlf);
  
  nlf = 0;
  start = get_wall_time();
  create_leafmtx(st_leafmtx,st_clt,st_clt,param,lnmtx,nffc,&nlf);
  end = get_wall_time();
  spent = end - start;
  printf("nlf:%d\n",nlf);
  printf("block cluster tree time spent:%.10f\n",spent);
  printf("depth_max:%d  count_node:%d\n",depth_max,count_node);
  /*qsort_row_leafmtx(st_leafmtx,0,nlf-1);
  int ilp = 0;
  int ips = 0;
  for(ip=0;ip<nlf;ip++){
    il = st_leafmtx[ip].nstrtl;
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

  FILE * fp;
  fp = fopen ("result_plot.txt","w");  
  for(i=0;i<st_leafmtxp->nlf;i++){
    fprintf(fp,"%ld, %ld, %ld, %ld, %ld, %ld\n",0,st_leafmtx[i].nstrtl,st_leafmtx[i].nstrtt,
	    st_leafmtx[i].nstrtl+st_leafmtx[i].ndl-1,st_leafmtx[i].nstrtt+st_leafmtx[i].ndt-1,st_leafmtx[i].ltmtx);
  }
  fclose (fp);
  free_st_clt(st_clt);

  for(il=0;il<nofc;il++){
    for(ig=0;ig<nffc;ig++){
      int is = ig + il * nffc;
      lod[is] = lodfc[il];
    }
    }*/
}

void qsort_row_leafmtx(leafmtx *st_leafmtx,int first,int last){
  int i,j,pivot;
  leafmtx st_www;
  if(first<last){
    pivot=first;
    i=first;
    j=last;
    while(i<j){
      while(st_leafmtx[i].nstrtl <= st_leafmtx[pivot].nstrtl && i<last)
        i++;
      while(st_leafmtx[j].nstrtl>st_leafmtx[pivot].nstrtl)
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

void qsort_col_leafmtx(leafmtx *st_leafmtx,int first,int last){
  int i,j,pivot;
  leafmtx st_www;
  if(first<last){
    pivot=first;
    i=first;
    j=last;
    while(i<j){
      while(st_leafmtx[i].nstrtt <= st_leafmtx[pivot].nstrtt && i<last)
	i++;
      while(st_leafmtx[j].nstrtt>st_leafmtx[pivot].nstrtt)
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

void create_leafmtx(leafmtx *st_leafmtx,cluster *st_cltl,cluster *st_cltt,
                    double param[],int *lnmtx,int nffc,int *nlf){
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
  if((st_cltl->zwdth * zeta <= zdistlt || st_cltt->zwdth * zeta <= zdistlt) && (ndl >= nleaf && ndt >= nleaf)){
    //st_leafmtx[*nlf] = (leafmtx *)malloc(sizeof(leafmtx));
    st_leafmtx[*nlf].nstrtl = nstrtl;
    st_leafmtx[*nlf].ndl = ndl;
    st_leafmtx[*nlf].nstrtt = nstrtt;
    st_leafmtx[*nlf].ndt = ndt;
    st_leafmtx[*nlf].kt = 0;
    st_leafmtx[*nlf].ltmtx = 1;
    *nlf = *nlf + 1;
  }else{
    if(nnsonl == 0 || nnsont == 0 || ndl <= nleaf || ndt <= nleaf){
      //st_leafmtx[*nlf]= (leafmtx *)malloc(sizeof(leafmtx));
      st_leafmtx[*nlf].nstrtl = nstrtl;
      st_leafmtx[*nlf].ndl = ndl;
      st_leafmtx[*nlf].nstrtt = nstrtt;
      st_leafmtx[*nlf].ndt = ndt;
      st_leafmtx[*nlf].ltmtx = 2;
      *nlf = *nlf + 1;
    }else{
      for(il=0;il<nnsonl;il++){
        for(it=0;it<nnsont;it++){
          create_leafmtx(st_leafmtx,st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lnmtx,nffc,nlf);
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
  if ((st_cltl->zwdth * zeta <= zdistlt || st_cltt->zwdth * zeta <= zdistlt) && (ndl >= nleaf && ndt >= nleaf)){
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

void checkClusterTree(cluster *st_clt){
  //if(st_clt->ndpth<11){
    printf("%d %d %d %d %lf\n",st_clt->nstrt,st_clt->nsize,st_clt->ndpth,st_clt->nnson,st_clt->zwdth);
    //}
  if(st_clt->nnson==0){
    //fprintf(f,"%lf %lf %lf %lf %lf %lf ",st_clt->bmin[0],st_clt->bmin[1],st_clt->bmin[2],st_clt->bmax[0],st_clt->bmax[1],st_clt->bmax[2]);
    //fprintf(f,"%d %d %d %d %lf\n",st_clt->nstrt,st_clt->nsize,st_clt->ndpth,st_clt->nnson,st_clt->zwdth);
    return;
  }else if(st_clt->nnson==1){
    checkClusterTree(st_clt->pc_sons[0]);
  }else{
    checkClusterTree(st_clt->pc_sons[0]);
    checkClusterTree(st_clt->pc_sons[1]);
  }
}

ompnode createOMPNode(cluster *st_clt,double (*zgmid)[3],double (*tempgmid)[3],int ndpth,int ndscd,int nsrt,int nd,int md,int ndim,int nclst,int nl){
  ompnode newnode;
  newnode.st_clt = st_clt;
  newnode.zgmid = zgmid;
  newnode.tempgmid = tempgmid;
  newnode.ndpth = ndpth;
  newnode.ndscd = ndscd;
  newnode.nsrt = nsrt;
  newnode.nd = nd;
  newnode.md = md;
  newnode.ndim = ndim;
  newnode.nclst = nclst;
  newnode.nl = nl;
  return newnode;
}

/****create cluster tree******/
cluster * create_ctree_ssgeom_t(cluster *st_clt,   //the current node
				double (*zgmid)[3],     //coordination of objects
				double (*tempgmid)[3],
				double param[],    //param[21] is minGroup, param[31] is ? 
				int ndpth,         //depth of the tree
				int ndscd,
				int nsrt,          //the start index of list
				int nd,            //the length of list
				int md,            //number of data
				int ndim,          //number of dimension (default is 3)
				int nclst){        //number of cluster
  int id,il,nson,ic,i,im;
  double minsz = param[21];
  double zcoef = param[31];
  double zlmin[ndim],zlmax[ndim];
  ndpth = ndpth + 1;
  
  if(nd <= minsz){
    nson = 0;
    st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
    if(ndpth > depth_max){
      depth_max = ndpth;
    }
    count_node++;
  }else{
    if(nd>PN && ndpth <= OMPS+1){
      double max_x = -DBL_MAX;
      double min_x = DBL_MAX;
      double max_y = -DBL_MAX;
      double min_y = DBL_MAX;
      double max_z = -DBL_MAX;
      double min_z = DBL_MAX;
#pragma omp parallel for reduction(max:max_x,max_y,max_z) reduction(min:min_x,min_y,min_z)
      for(il=0;il<nd;il++){
	if(zgmid[il][0] > max_x){
	  max_x = zgmid[il][0];
	}
	if(zgmid[il][0] < min_x){
	  min_x = zgmid[il][0];
	}
	if(zgmid[il][1] > max_y){
	  max_y = zgmid[il][1];
	}
	if(zgmid[il][1] < min_y){
	  min_y = zgmid[il][1];
	}
	if(zgmid[il][2] > max_z){
	  max_z = zgmid[il][2];
	}
	if(zgmid[il][2] < min_z){
	  min_z = zgmid[il][2];
	}
      }
      zlmin[0] = min_x;
      zlmax[0] = max_x;
      zlmin[1] = min_y;
      zlmax[1] = max_y;
      zlmin[2] = min_z;
      zlmax[2] = max_z;
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

    if(nd > PN && ndpth <= OMPS+1){
      int gn = omp_get_max_threads();
      int chunk_size = nd / gn;
      int *lessNum = (int *)malloc(gn*sizeof(int));
      int *moreNum = (int *)malloc(gn*sizeof(int));
      int *lessStart = (int *)malloc(gn*sizeof(int));
      int *moreStart = (int *)malloc(gn*sizeof(int));
      int ssum = 0;
#pragma omp parallel
      {
	int id = omp_get_thread_num();
	int ln = 0, mn = 0;
	int start = id * chunk_size;
	int end = id == gn-1 ? nd : start+chunk_size;
	for(im=start;im<end;im++){
	  if(zgmid[im][ncut]<=zlmid){
	    ln++;
	  }else{
	    mn++;
	  }
	}
	lessNum[id] = ln;
	moreNum[id] = mn;
      }
      for(il=0;il<gn;il++){
	ssum += lessNum[il];
      }
      nl = ssum;

      lessStart[0] = 0;
      moreStart[0] = nl;
      if(nl != 0 && nl != nd){
	int tl = 0, tm = nl;
	for(il=0;il<gn-1;il++){
	  tl = tl + lessNum[il];
	  tm = tm + moreNum[il];
	  lessStart[il+1] = tl;
	  moreStart[il+1] = tm;
	}

#pragma omp parallel
	{
	  int id = omp_get_thread_num();
	  int ls = lessStart[id];
	  int ms = moreStart[id];
	  int start = id * chunk_size;
	  int end = id == gn-1 ? nd : start+chunk_size;
	  for(im=start;im<end;im++){
	    if(zgmid[im][ncut]<=zlmid){
	      tempgmid[ls][0] = zgmid[im][0];
	      tempgmid[ls][1] = zgmid[im][1];
	      tempgmid[ls][2] = zgmid[im][2];
	      ls++;
	    }else{
	      tempgmid[ms][0] = zgmid[im][0];
	      tempgmid[ms][1] = zgmid[im][1];
	      tempgmid[ms][2] = zgmid[im][2];
	      ms++;
	    }
	  }
	}
      }
      free(lessNum);
      free(moreNum);
      free(lessStart);
      free(moreStart);
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
      nson = 0;
      st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
      if(ndpth > depth_max){
        depth_max = ndpth;
      }
      count_node++;
    }else{
      nson = 2;
      st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
      if(ndpth > depth_max){
	depth_max = ndpth;
      }
      count_node++;

      if(ndpth == OMPS+1){
	if((ndpth & 1) && nd > PN){
	  for(i=0;i<nd;i++){
	    zgmid[i][0] = tempgmid[i][0];
	    zgmid[i][1] = tempgmid[i][1];
	    zgmid[i][2] = tempgmid[i][2];
	  }
	  nodelist[nodecount] = createOMPNode(st_clt,zgmid,tempgmid,ndpth,ndscd,nsrt,nd,md,ndim,nclst,nl);
	}else{
	  nodelist[nodecount] = createOMPNode(st_clt,tempgmid,zgmid,ndpth,ndscd,nsrt,nd,md,ndim,nclst,nl);
	}
	nodecount++;
      }else{
	if(ndpth > OMPS+1 || nd <= PN){
	  int nsrt1 = nsrt;
	  int nd1 = nl;
	  st_clt->pc_sons[0] = create_ctree_ssgeom_t(st_clt->pc_sons[0],zgmid,tempgmid,param,
						     ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
	  nsrt1 = nsrt + nl;
	  nd1 = nd - nl;
	  st_clt->pc_sons[1] = create_ctree_ssgeom_t(st_clt->pc_sons[1],&zgmid[nl],&tempgmid[nl],param,
						     ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
	}else{
	  int nsrt1 = nsrt;
          int nd1 = nl;
	  if((ndpth & 1) && nd1 <= PN){
	    for(i=0;i<nl;i++){
	      zgmid[i][0] = tempgmid[i][0];
	      zgmid[i][1] = tempgmid[i][1];
	      zgmid[i][2] = tempgmid[i][2];
	    }
	    st_clt->pc_sons[0] = create_ctree_ssgeom_t(st_clt->pc_sons[0],zgmid,tempgmid,param,
						       ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
          }else{
	    st_clt->pc_sons[0] = create_ctree_ssgeom_t(st_clt->pc_sons[0],tempgmid,zgmid,param,
                                                       ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
	  }
	  nsrt1 = nsrt + nl;
          nd1 = nd - nl;
	  if((ndpth & 1) && nd1 <= PN){
	    for(i=nl;i<nd;i++){
	      zgmid[i][0] = tempgmid[i][0];
	      zgmid[i][1] = tempgmid[i][1];
	      zgmid[i][2] = tempgmid[i][2];
	    }
	    st_clt->pc_sons[1] = create_ctree_ssgeom_t(st_clt->pc_sons[1],&zgmid[nl],&tempgmid[nl],param,
                                                       ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
	  }else{
	    st_clt->pc_sons[1] = create_ctree_ssgeom_t(st_clt->pc_sons[1],&tempgmid[nl],&zgmid[nl],param,
						       ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
	  }
	}
      }
    }
  }
  st_clt->ndscd = nd;
  return st_clt;
}

void create_ctree_ssgeom_b(double param[]){
  int i;
  int nsrt1,nd1;
#pragma omp parallel for private(nsrt1,nd1) schedule(dynamic,1)
  for(i=0;i<nodecount;i++){
    nsrt1 = nodelist[i].nsrt;
    nd1 = nodelist[i].nl;
    nodelist[i].st_clt->pc_sons[0] = create_ctree_ssgeom_t(nodelist[i].st_clt->pc_sons[0],nodelist[i].zgmid,nodelist[i].tempgmid,
							   param,nodelist[i].ndpth,nodelist[i].ndscd,nsrt1,nd1,nodelist[i].md,nodelist[i].ndim,nodelist[i].nclst);
    nsrt1 = nodelist[i].nsrt + nodelist[i].nl;
    nd1 = nodelist[i].nd - nodelist[i].nl;
    nodelist[i].st_clt->pc_sons[1] = create_ctree_ssgeom_t(nodelist[i].st_clt->pc_sons[1],&nodelist[i].zgmid[nodelist[i].nl],&nodelist[i].tempgmid[nodelist[i].nl],
							   param,nodelist[i].ndpth,nodelist[i].ndscd,nsrt1,nd1,nodelist[i].md,nodelist[i].ndim,nodelist[i].nclst);
  }
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
