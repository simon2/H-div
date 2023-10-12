#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mkl.h>
#include "data/bem_file.h"

#ifdef DEBUG
#include <assert.h>
#endif

// #define INPUT_DEFAULT "bem_data/input_10ts.txt"             //small data for dbg
// #define INPUT_DEFAULT "bem_data/input_100ts.txt"          //Sphere
// #define INPUT_DEFAULT "bem_data/input_10ts_c30_3_3_4.txt"   //SphereCube
#define INPUT_DEFAULT "bem_data/input_10ts_p30_4.txt"       //SpherePyramid
// #define INPUT_DEFAULT "bem_data/input_216h_5x10.txt"        //Humans
// #define INPUT_DEFAULT "bem_data/input_84tp7_30_2p.txt"      //SpherePyramidPyramid


/*********define cluster************/
typedef struct cluster cluster;
struct cluster{
  int ndim;
  int nstrt,nsize,ndpth,nnson,nmbr;
  int ndscd;                         //number of descendants
  double bmin[3];
  double bmax[3];
  double zwdth;
  int offsets[2];
  int nnsons;
  long nnnd;
};

/*********define leafMtx***********/ //information of H-matrix
typedef struct leafmtx leafmtx;
struct leafmtx{
  int ltmtx;                         //kind of the matrix; 1:rk 2:full
  int kt;                            //rank of the partition
  int nstrtl,ndl;                    //the coordination of the first element of partition
  int nstrtt,ndt;                    //the length & width of the partition
  double *a1;
  double *a2;                    //the elements of the partition, a2 is row, a1 is col
  int ndpth;
};

/*********define leafmtxp*********/  //whole H-matrix
typedef struct leafmtxp leafmtxp;
struct leafmtxp{
  leafmtx* st_leafmtx;               //list of partitions (BCT leaf-nodes)
  int nlf;                           //number of partitions
  int nlfkt;                         //number ot partitions approximated
};

void supermatrix_construction_cog_leafmtrx(leafmtxp *st_leafmtxp,double param[],int *lod,int *lnmtx,int nofc,int nffc,int ndim);
int med3(int nl,int nr,int nlr2);
void qsort_col_leafmtx(leafmtx *st_leafmtx,int first,int last);
void qsort_row_leafmtx(leafmtx *st_leafmtx,int first,int last);
void create_leafmtx(leafmtx *st_leafmtx,int st_cltl,int st_cltt,double param[],int *lnmtx,int nffc,int *nlf);
double dist_2cluster(int st_cltl,int st_cltt);
void count_lntmx(int st_cltl,int st_cltt,double param[],int *lnmtx,int nffc);
int create_cluster(int nmbr,int ndpth,int nstrt,int nsize,int ndim,int nson);
int create_ctree_ssgeom(int st_clt,double (*zgmid)[3],int (*face2node)[3],double param[],int ndpth,int ndscd,int nsrt,int nd,int md,int ndim,int nclst);

int acaplus(double* zaa, double* zab, int ndl, int ndt, int nstrtl, int nstrtt, int kmax, double eps, double znrmmat, double pACA_EPS);
void fill_leafmtx(leafmtx *st_lf, double znrmmat, int *lnmtx, int nd, int nlf);
int minabsvalloc_d(double* za, int nd);
int maxabsvalloc_d(double* za, int nd);
double entry_ij(int i, int j);
double face_integral2(double xs[], double ys[], double zs[], double x, double y, double z);
void comp_col(double* zaa, double *zab, int ndl, int ndt, int k, int it, double* col, int nstrtl, int nstrtt, int* lrow_done);
void comp_row(double* zaa, double* zab, int ndl, int ndt, int k, int il, double* row, int nstrtl, int nstrtt, int* lrow_done);
void cross_product(double* u, double* v, double* w);
double dot_product(double* v, double* u, int n);
void adotsub_dsm(double* zr, double* zaa, double* zu, int it, int ndl, int ndt, int mdl, int mdt);
void adot_dsm(double* zau, double* zaa, double* zu, int im, int ndl, int ndt, int mdl, int mdt);
int max(int a, int b);
int min(int a, int b);

double get_wall_time();
double get_cpu_time();

void checkClusterTree(int st_clt);
void checkClusterTree2(int ndpth, int st_clt);
int depth_max;
int count_node;

double (*zgmid)[3];
double (*bgmid)[3];
int (*f2n)[3];

int LN = 30000000;
cluster* CTlist;
int countCT = 0; 

double middle;
int middle_flag = 0;

int bal[32];
int main(int argc, char **argv){
  /******** read file *********/
  char *fname;
  FILE *file;
  int countOfNode=0;
  int count = 0;
  int i;
  double (*coordOfNode)[3];
  double (*coordOfFace)[3];
  int (*face2node)[3];
  struct bem_input bi;
  fname = (argc >= 2)?argv[1]:INPUT_DEFAULT;
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
    bgmid = bi.coordOfNode;
    count = bi.nFace;
    zgmid = bi.coordOfFace;
    f2n = bi.face2node;
  }
  fclose(file);
  double param[100];
  for(i=0;i<100;i++){
    param[i] = 0;
  }
  param[21] = 50.0;
  param[31] = 1.1;
  param[41] = 40.0;
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
  double start = get_wall_time();
  supermatrix_construction_cog_leafmtrx(st_leafmtxp,param,lod,lnmtx,nofc,nffc,ndim);  // construction of Leaf-matrix 
  double end = get_wall_time();
  printf("MP time:%f\n", (end - start));
  return 0;
}

#ifdef DEBUG
int *lod0;
#endif

void supermatrix_construction_cog_leafmtrx(leafmtxp *st_leafmtxp,    //the H-matrix
					   double param[],int *lod,
					   int *lnmtx,              //1:k-rank 2:dense 3:H-matrix
					   int nofc,int nffc,       //number of elements in same coordination
					   int ndim){
  int st_clt = 0;
  int i,nfl,nflkt,ip,il,ig;
  int nd = nofc * nffc;
  int *lodfc;
  leafmtx *st_leafmtx;
  lodfc = (int *)malloc(nofc*sizeof(int));
  for(il=0;il<nofc;il++){
    lodfc[il] = il;
  }
  int nsrt = 0;
  int ndf = nofc;
  int nclst = 0;
  int ndpth = 0;
  int ndscd = 0;

  depth_max = 0;
  count_node = 0;

  double start,end,spent;
  double start2,end2,spent2;
  start = get_wall_time();
#ifdef DEBUG
  lod0 = lodfc;
#endif 

  CTlist = (cluster*)malloc(sizeof(cluster) * LN);
  st_clt = create_ctree_ssgeom(st_clt,zgmid,f2n,param,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst);
  end = get_wall_time();
  spent = end - start;
  printf("cluster tree time spent:%.10f\n",spent);
  printf("breakdown:%.10f\n",middle-start);

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
  start2 = get_wall_time();
  create_leafmtx(st_leafmtx,st_clt,st_clt,param,lnmtx,nffc,&nlf);
  end2 = get_wall_time();
  spent2 = end - start;
  printf("nlf:%d\n",nlf);
  printf("block cluster tree time spent:%.10f\n",spent2);
  printf("depth_max:%d  count_node:%d\n",depth_max,count_node);

  qsort_row_leafmtx(st_leafmtx,0,nlf-1);
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

  start = get_wall_time();
  fill_leafmtx(st_leafmtx, 0.0, lnmtx, ndf, nlf);
  end = get_wall_time();
  spent = end - start;
  spent2 = end - start2;
  printf("filling time:%.10f\n",spent);
  printf("BCT + filling time:%.10f\n",spent2);
  printf("nlf:%d\n",nlf);
  st_leafmtxp->st_leafmtx = st_leafmtx;
  st_leafmtxp->nlf = nlf;
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

void create_leafmtx(leafmtx *st_leafmtx,int st_cltl,int st_cltt,
                    double param[],int *lnmtx,int nffc,int *nlf){
  int ndl = CTlist[st_cltl].nsize * nffc;
  int ndt = CTlist[st_cltt].nsize * nffc;
  int nstrtl = CTlist[st_cltl].nstrt;
  int nstrtt = CTlist[st_cltt].nstrt;
  int nnsonl = CTlist[st_cltl].nnson;
  int nnsont = CTlist[st_cltt].nnson;
  int ndpth = CTlist[st_cltt].ndpth;
  int il,it;

  double nleaf = param[41];
  double zeta = param[51];
  double zdistlt = dist_2cluster(st_cltl,st_cltt);
  if((CTlist[st_cltl].zwdth <= zdistlt * zeta || CTlist[st_cltt].zwdth <= zdistlt * zeta) && (ndl >= nleaf && ndt >= nleaf)){
    st_leafmtx[*nlf].nstrtl = nstrtl;
    st_leafmtx[*nlf].ndl = ndl;
    st_leafmtx[*nlf].nstrtt = nstrtt;
    st_leafmtx[*nlf].ndt = ndt;
    st_leafmtx[*nlf].kt = 0;
    st_leafmtx[*nlf].ltmtx = 1;
    st_leafmtx[*nlf].ndpth = ndpth;
    *nlf = *nlf + 1;
  }else{
    if(nnsonl == 0 || nnsont == 0 || ndl <= nleaf || ndt <= nleaf){
      st_leafmtx[*nlf].nstrtl = nstrtl;
      st_leafmtx[*nlf].ndl = ndl;
      st_leafmtx[*nlf].nstrtt = nstrtt;
      st_leafmtx[*nlf].ndt = ndt;
      st_leafmtx[*nlf].ltmtx = 2;
      st_leafmtx[*nlf].ndpth = ndpth;
      *nlf = *nlf + 1;
    }else{
      for(il=0;il<nnsonl;il++){
        for(it=0;it<nnsont;it++){
          create_leafmtx(st_leafmtx,CTlist[st_cltl].offsets[il],CTlist[st_cltt].offsets[it],param,lnmtx,nffc,nlf);
        }
      }
    }
  }
}

double dist_2cluster(int st_cltl,int st_cltt){
  double zs = 0.0;
  int id;
  for(id=0;id<CTlist[st_cltl].ndim;id++){
    if(CTlist[st_cltl].bmax[id] < CTlist[st_cltt].bmin[id]){
      zs = zs + (CTlist[st_cltt].bmin[id] - CTlist[st_cltl].bmax[id]) * (CTlist[st_cltt].bmin[id] - CTlist[st_cltl].bmax[id]);
    }else if(CTlist[st_cltt].bmax[id]< CTlist[st_cltl].bmin[id]){
      zs = zs+ (CTlist[st_cltl].bmin[id] - CTlist[st_cltt].bmax[id]) * (CTlist[st_cltl].bmin[id] - CTlist[st_cltt].bmax[id]);
    }
  }
  return sqrt(zs);
}

void count_lntmx(int st_cltl,int st_cltt,double param[],int *lnmtx,int nffc){
  int il,it;
  int ndl = CTlist[st_cltl].nsize * nffc;
  int ndt = CTlist[st_cltt].nsize * nffc;
  int nstrtl = CTlist[st_cltl].nstrt;
  int nstrtt = CTlist[st_cltt].nstrt;
  int nnsonl = CTlist[st_cltl].nnson;
  int nnsont = CTlist[st_cltt].nnson;

  double nleaf = param[41];
  double zeta = param[51];
  double zdistlt = dist_2cluster(st_cltl,st_cltt);
  if ((CTlist[st_cltl].zwdth <= zdistlt * zeta || CTlist[st_cltt].zwdth <= zdistlt * zeta) && (ndl >= nleaf && ndt >= nleaf)){
    lnmtx[0] = lnmtx[0] + 1;
  }else{
    if(nnsonl == 0 || nnsont == 0 || ndl <= nleaf || ndt <= nleaf){
      lnmtx[1] = lnmtx[1] + 1;
    }else{
      lnmtx[2] = lnmtx[2] + 1;
      for(il=0;il<nnsonl;il++){
	      for(it=0;it<nnsont;it++){
	        count_lntmx(CTlist[st_cltl].offsets[il],CTlist[st_cltt].offsets[it],param,lnmtx,nffc);
	      }
      }
    }
  }
}

/*****create a cluster*********/
int create_cluster(int nmbr,int ndpth,int nstrt,int nsize,int ndim,int nson){
  int st_clt;
  st_clt = countCT;
  nmbr = nmbr + 1;
  countCT++;
  CTlist[st_clt].nstrt = nstrt;
  CTlist[st_clt].nsize = nsize;
  CTlist[st_clt].ndim = ndim;
  CTlist[st_clt].nnson = nson;
  CTlist[st_clt].nmbr = nmbr;
  CTlist[st_clt].ndpth = ndpth;

  return st_clt;
}

/****create cluster tree******/
int create_ctree_ssgeom(int st_clt,   //the current node
			      double (*zgmid)[3],     //coordination of objects
            int (*face2node)[3],
			      double param[],    //param[21] is minGroup, param[31] is ? 
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
#ifdef DEBUG
  printf ("ndepth = %ld, nd = %ld, range = %ld-%ld ", ndpth, nd, lod-lod0, lod-lod0+nd);
#endif
  
  if(nd <= minsz){
#ifdef DEBUG
    fprintf (stdout, "... leaf\n");
#endif
    nson = 0;
    st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
    if(ndpth > depth_max){
      depth_max = ndpth;
    }
    if(nd > bal[ndpth]){
      bal[ndpth] = nd;
    }
    count_node++;
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
        for(id=0;id<ndim;id++){
          int mh = face2node[nl][id];
          face2node[nl][id] = face2node[nr][id];
          face2node[nr][id] = mh;
        }
      }
    }

#ifdef DEBUG
    fprintf (stdout, "nl = %ld, nr = %ld\n", nl, nr);
    assert ( nl-nr == 1 );
    // assert ( nl < nd );
#endif
    
    if(nl == nd || nl == 0){
      fprintf (stdout, "nl = %ld, nr = %ld\n", nl, nr);
      nson = 0;
      st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
      if(ndpth > depth_max){
        depth_max = ndpth;
      }
      if(nd > bal[ndpth]){
	bal[ndpth] = nd;
      }
      count_node++;
    }else{
      if(middle_flag == 0){
        middle = get_wall_time();
        middle_flag++;
      }
      nson = 2;
      st_clt = create_cluster(nclst,ndpth,nsrt,nd,ndim,nson);
      if(ndpth > depth_max){
	depth_max = ndpth;
      }
      if(nd > bal[ndpth]){
	bal[ndpth] = nd;
      }
      count_node++;
      int nsrt1 = nsrt;
      int nd1 = nl;
      CTlist[st_clt].offsets[0] = create_ctree_ssgeom(CTlist[st_clt].offsets[0],zgmid,face2node,param,
					       ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);

      nsrt1 = nsrt + nl;
      nd1 = nd - nl;
      CTlist[st_clt].offsets[1] = create_ctree_ssgeom(CTlist[st_clt].offsets[1],&zgmid[nl],&face2node[nl],param,
					       ndpth,ndscd,nsrt1,nd1,md,ndim,nclst);
    }
  }
  CTlist[st_clt].ndscd = nd;
  //bounding box
  double zeps = 1.0e-5;
  if(CTlist[st_clt].nnson > 0){
    for(id=0;id<ndim;id++){
      CTlist[st_clt].bmin[id] = CTlist[CTlist[st_clt].offsets[0]].bmin[id];
      CTlist[st_clt].bmax[id] = CTlist[CTlist[st_clt].offsets[0]].bmax[id];
    }
    for(il=1;il<CTlist[st_clt].nnson;il++){
      for(id=0;id<ndim;id++){
	      if(CTlist[CTlist[st_clt].offsets[il]].bmin[id] < CTlist[st_clt].bmin[id]){
	        CTlist[st_clt].bmin[id] = CTlist[CTlist[st_clt].offsets[il]].bmin[id];
	      }
	      if(CTlist[st_clt].bmax[id] < CTlist[CTlist[st_clt].offsets[il]].bmax[id]){
	        CTlist[st_clt].bmax[id] = CTlist[CTlist[st_clt].offsets[il]].bmax[id];
	      }
      }
    }
  }else{
    for(id=0;id<ndim;id++){
      CTlist[st_clt].bmin[id] = zgmid[0][id];
      CTlist[st_clt].bmax[id] = zgmid[0][id];
    }
    for(id=0;id<ndim;id++){
      for(il=1;il<CTlist[st_clt].nsize;il++){
	      if(zgmid[il][id] < CTlist[st_clt].bmin[id]){
	        CTlist[st_clt].bmin[id] = zgmid[il][id];
	      }
	      if(CTlist[st_clt].bmax[id] < zgmid[il][id]){
	        CTlist[st_clt].bmax[id] = zgmid[il][id];
	      }
      }
    }
  }
  double zwdth = (CTlist[st_clt].bmax[0] - CTlist[st_clt].bmin[0]) * (CTlist[st_clt].bmax[0] - CTlist[st_clt].bmin[0]);
  for(id=1;id<ndim;id++){
    zwdth = zwdth + (CTlist[st_clt].bmax[id] - CTlist[st_clt].bmin[id]) * (CTlist[st_clt].bmax[id] - CTlist[st_clt].bmin[id]);
  }
  zwdth = sqrt(zwdth);
  for(id=0;id<ndim;id++){
    double bdiff = CTlist[st_clt].bmax[id] - CTlist[st_clt].bmin[id];
    if(bdiff < zeps * zwdth){
      CTlist[st_clt].bmax[id] = CTlist[st_clt].bmax[id] + 0.5 * (zeps * zwdth - bdiff);
      CTlist[st_clt].bmin[id] = CTlist[st_clt].bmin[id] - 0.5 * (zeps * zwdth - bdiff);
    }
  }
  zwdth = (CTlist[st_clt].bmax[0] - CTlist[st_clt].bmin[0]) * (CTlist[st_clt].bmax[0] - CTlist[st_clt].bmin[0]);
  for(id=1;id<ndim;id++){
    zwdth += (CTlist[st_clt].bmax[id] - CTlist[st_clt].bmin[id]) * (CTlist[st_clt].bmax[id] - CTlist[st_clt].bmin[id]);
  }
  CTlist[st_clt].zwdth = sqrt(zwdth);
  //end of bounding box
  return st_clt;
}

int acaplus(double* zaa, double* zab, int ndl, int ndt, int nstrtl, int nstrtt, int kmax, double eps, double znrmmat, double pACA_EPS){
  double *prow, *pcol, *pb_ref, *pa_ref;
  int *lrow_done, *lcol_done;

  double (*zaa2)[ndl];
  double (*zab2)[ndt];
  zaa2 = (double(*)[ndl])zaa;
  zab2 = (double(*)[ndt])zab;
  
  int il,it,ib;

  int INCY = 1;
  double za_ACA_EPS = 1.0e-10;

  double znrm = znrmmat * sqrt((double)ndl * (double)ndt);
  double ACA_EPS = pACA_EPS;

  int ntries = max(ndl, ndt) + 1;
  int ntries_row = 6;
  int ntries_col = 6;

  lrow_done = (int *)calloc(ndl, sizeof(int));
  lcol_done = (int *)calloc(ndt, sizeof(int));
  int k = 0;

  int j_ref = 0;
  pa_ref = (double *)malloc(ndl * sizeof(double));

  comp_col(zaa, zab, ndl, ndt, k, j_ref, pa_ref, nstrtl, nstrtt, lrow_done);

  double colnorm = cblas_dnrm2(ndl, pa_ref, INCY);

  int i_ref = minabsvalloc_d(pa_ref, ndl);
  double rownorm = fabs(pa_ref[i_ref]);
  pb_ref = (double *)malloc(ndt * sizeof(double));
  comp_row(zaa, zab, ndl, ndt, k, i_ref, pb_ref, nstrtl, nstrtt, lcol_done);

  rownorm = cblas_dnrm2(ndt, pb_ref, INCY);

  double apxnorm = 0.0;
  int lstop_aca = 0;
  
  double col_maxval, row_maxval;
  while(k < kmax && (ntries_row > 0 || ntries_col > 0) && ntries > 0){
    ntries--;
    pcol = &zaa2[k][0];
    prow = &zab2[k][0];
    col_maxval = 0.0;
    int i = maxabsvalloc_d(pa_ref, ndl);
    col_maxval = fabs(pa_ref[i]);
    row_maxval = 0.0;
    int j = maxabsvalloc_d(pb_ref, ndt);
    row_maxval = fabs(pb_ref[j]);

    double zinvmax;
    if(row_maxval > col_maxval){
      if(j != j_ref){
        comp_col(zaa, zab, ndl, ndt, k, j, pcol, nstrtl, nstrtt, lrow_done);
      }else{
        for(il=0;il<ndl;il++){
          pcol[il] = pa_ref[il];
        }
      }
      i = maxabsvalloc_d(pcol, ndl);
      col_maxval = fabs(pcol[i]);

      if(col_maxval < ACA_EPS && k >= 1){
        lstop_aca = 1;
      }else{
        comp_row(zaa, zab, ndl, ndt, k, i, prow, nstrtl, nstrtt, lcol_done);
        if(fabs(pcol[i]) > 1.0e-20){
          zinvmax = 1.0 / pcol[i];
        }else{
          k = max(k-1, 0);
          break;
        }
        for(il=0;il<ndl;il++){
          pcol[il] *= zinvmax;
        } 
      }
    }else{
      if(i != i_ref){
        comp_row(zaa, zab, ndl, ndt, k, i, prow, nstrtl, nstrtt, lcol_done);
      }else{
        for(il=0;il<ndt;il++){
          prow[il] = pb_ref[il];
        }
      }
      j = maxabsvalloc_d(prow, ndt);
      row_maxval = fabs(prow[j]);

      if(row_maxval < ACA_EPS && k >= 1){
        lstop_aca = 1;
      }else{
        comp_col(zaa, zab, ndl, ndt, k, j, pcol, nstrtl, nstrtt, lrow_done);
        if(fabs(prow[j]) > 1.0e-20){
          zinvmax = 1.0 / prow[j];
        }else{
          k = max(k-1, 0);
          break;
        }
        for(il=0;il<ndt;il++){
          prow[il] *= zinvmax;
        } 
      }
    }
    lrow_done[i] = 1;
    lcol_done[j] = 1;

    if(i != i_ref){
      zinvmax = -pcol[i_ref];
      for(il=0;il<ndt;il++){
        pb_ref[il] += prow[il] * zinvmax;
      }
      rownorm = cblas_dnrm2(ndt, pb_ref, INCY);
    }
    if(i == i_ref || rownorm < ACA_EPS){
      if(i == i_ref){
        ntries_row++;
      }
      if(ntries_row > 0){
        rownorm = 0.0;
        i = i_ref;
        while(i != (i_ref + ndl - 1) % ndl && rownorm < za_ACA_EPS && ntries_row > 0){
          if(lrow_done[i] == 0){
            comp_row(zaa, zab, ndl, ndt, k+1, i, pb_ref, nstrtl, nstrtt, lcol_done);
            rownorm = cblas_dnrm2(ndt, pb_ref, INCY);
            if(rownorm < ACA_EPS){
              lrow_done[i] = 1;
            }
            ntries_row--;
          }else{
            rownorm = 0.0;
          }
          i = (i+1) % ndl;
        }
        i_ref = (i + ndl - 1) % ndl;
      }
    }

    if(j != j_ref){
      zinvmax = -prow[j_ref];
      for(il=0;il<ndl;il++){
        pa_ref[il] += pcol[il] * zinvmax;
      }
      colnorm = cblas_dnrm2(ndl, pa_ref, INCY);
    }
    if(j == j_ref || colnorm < ACA_EPS){
      if(j == j_ref){
        ntries_col++;
      }
      if(ntries_col > 0){
        colnorm = 0.0;
        j = j_ref;
        while(j != (j_ref + ndt - 1) % ndt && colnorm < za_ACA_EPS && ntries_col > 0){
          if(lcol_done[j] == 0){
            comp_col(zaa, zab, ndl, ndt, k+1, j, pa_ref, nstrtl, nstrtt, lrow_done);
            colnorm = cblas_dnrm2(ndl, pa_ref, INCY);
            if(colnorm < ACA_EPS){
              lcol_done[j] = 1;
            }
            ntries_col--;
          }else{
            colnorm = 0.0;
          }
          j = (j+1) % ndt;
        }
        j_ref = (j + ndt - 1) % ndt;
      }
    }

    if(colnorm < ACA_EPS && rownorm < ACA_EPS && k >= 1){
      lstop_aca = 1;
      k = k + 1;
    }

    if(lstop_aca == 0){
      double blknorm = cblas_dnrm2(ndl, pcol, INCY) * cblas_dnrm2(ndt, prow, INCY);
      if(k == 0){
        apxnorm = blknorm;
      }else{
        double compared = apxnorm * eps;
        if(blknorm < compared && rownorm < compared && colnorm < compared && k >= 1){
          lstop_aca = 1;
        }
      }
    }
    if(lstop_aca == 1 && k >= 1){
      break;
    }

    k++;
  }

  if(k < 1){
    printf("alert!\n");
    printf("colnorm=%f rownorm=%f ACA_EPS=%f\n", colnorm, rownorm, ACA_EPS);
    printf("col_maxval=%f row_maxval=%f\n", col_maxval, row_maxval);
    printf("ntries_row=%d ntries_col=%d ntries=%d\n", ntries_row, ntries_col, ntries);
    printf("k=%d\n", k);
  }

  return k;

}

void fill_leafmtx(leafmtx *st_lf, double znrmmat, int *lnmtx, int nd, int nlf){
  double start,end,spent;
  
  double eps = 1.0e-8;
  double ACA_EPS = 0.9 * eps;
  int kparam = 50;
  int ip,il,it;

  int k_max = -1;
  int k_min = kparam;
  long k_sum = 0;
  int k_count = 0;

  long sum_full_entries = 0;
  long sum_acap_entries = 0;
  double sum_time_full = 0.0;
  double sum_time_acap = 0.0;

  for(ip=0;ip<nlf;ip++){
    
    int ndl = st_lf[ip].ndl;
    int ndt = st_lf[ip].ndt;
    int ns = ndl * ndt;
    int nstrtl = st_lf[ip].nstrtl;
    int nstrtt = st_lf[ip].nstrtt;
    int ltmtx = st_lf[ip].ltmtx;
    

    if(ltmtx == 1){
      start = get_wall_time();

      st_lf[ip].a1 = (double*)malloc(sizeof(double) * ndt * kparam);
      st_lf[ip].a2 = (double*)malloc(sizeof(double) * ndl * kparam);
      if(!st_lf[ip].a1 && !st_lf[ip].a2){
        printf("allocate a1 or a2 failed!\n");
        exit(99);
      }

      int kt = acaplus(st_lf[ip].a2, st_lf[ip].a1, ndl, ndt, nstrtl, nstrtt, kparam, eps, znrmmat, ACA_EPS);
      st_lf[ip].kt = kt;
      // printf("DEBUGING: kt=%d, kparam=%d, nstrtl=%d, nstrtt=%d, ndl=%d, ndt=%d\n", kt, kparam, nstrtl, nstrtt, ndl, ndt);
      
      if(kt > k_max) k_max = kt;
      else if(kt < k_min) k_min = kt;
      k_sum += (long)kt;
      k_count++;

      if(kt > kparam){ 
        printf("WARNING: Insufficient k: kt=%d, kparam=%d, nstrtl=%d, nstrtt=%d, ndl=%d, ndt=%d\n", kt, kparam, nstrtl, nstrtt, ndl, ndt);
      }

      st_lf[ip].a1 = (double *)realloc(st_lf[ip].a1, kt * ndt * sizeof(double));
      st_lf[ip].a2 = (double *)realloc(st_lf[ip].a2, kt * ndl * sizeof(double));

      end = get_wall_time();
      spent = end - start;
      sum_time_acap += spent;
      int t_entries = (kt * (ndt + ndl));
      sum_acap_entries += t_entries;
      printf("low-time: %lf\n low-entires: %d,\n time_per_entries:%lf\n",spent, t_entries, (spent/(double)t_entries));
      // if(end-start > 1.0)
      //   printf("Low-rank filling time: k=%d nstrtl=%d, nstrtt=%d, ndl=%d, ndt=%d time=%E\n", kt, nstrtl, nstrtt, ndl, ndt, end-start);
    }else if(st_lf[ip].ltmtx == 2){
      start = get_wall_time();

      st_lf[ip].a1 = (double *)malloc(sizeof(double) * ns);
      if(!st_lf[ip].a1){
        printf("allocate a1 failed!\n");
        exit(99);
      }

      double (*tempa1)[ndt];
      tempa1 = (double(*)[ndt])st_lf[ip].a1;
  
      for(il=0;il<ndl;il++){
        int ill = il + nstrtl;
        for(it=0;it<ndt;it++){
          int itt = it + nstrtt;
          tempa1[il][it] = entry_ij(ill, itt);
        }
      }
      end = get_wall_time();
      spent = end - start;
      sum_time_full += spent;
      int t_entries = (ndt * ndl);
      sum_full_entries += t_entries;
      // printf("full\n");
      printf("full-time: %lf,\n full-entires: %d,\n time_per_entries:%lf\n",spent, t_entries, (spent/(double)t_entries));
      // if(end-start > 1.0)
      //   printf("Dense mtx filling time:  nstrtl=%d, nstrtt=%d, ndl=%d, ndt=%d time=%E\n", nstrtl, nstrtt, ndl, ndt, end-start);
    }
  }
  double k_average = (double)k_sum / k_count;
  printf("k_max:%d\nk_min:%d\nk_average:%f\n",k_max,k_min,k_average);
  printf("acap entries:%d\tacap time:%lf\tfull entries:%d\tfull time:%lf\n",sum_acap_entries,sum_time_acap,sum_full_entries,sum_time_full);
}

void comp_row(double* zaa, double* zab, int ndl, int ndt, int k, int il, double* row, int nstrtl, int nstrtt, int* lrow_done){
  int it;
  for(it=0;it<ndt;it++){
    if(lrow_done[it] == 0){
      int ill = il + nstrtl;
      int itt = it + nstrtt;
      row[it] = entry_ij(ill, itt);
    }
  }

  if(k == 0){
    return;
  }

  adotsub_dsm(row, zab, zaa, il, ndt, k, ndt, ndl);

  for(it=0;it<ndt;it++){
    if(lrow_done[it] != 0){
      row[it] = 0.0;
    }
  }
}

void comp_col(double* zaa, double* zab, int ndl, int ndt, int k, int it, double* col, int nstrtl, int nstrtt, int* lrow_done){
  int il;

  for(il=0;il<ndl;il++){
    if(lrow_done[il] == 0){
      int ill = il + nstrtl;
      int itt = it + nstrtt;
      col[il] = entry_ij(ill, itt);
    }
  }

  if(k == 0){
    return;
  }

  adotsub_dsm(col, zaa, zab, it, ndl, k, ndl, ndt);

  for(il=0;il<ndl;il++){
    if(lrow_done[il] != 0){
      col[il] = 0.0;
    }
  }
}

int minabsvalloc_d(double* za, int nd){
  int il = 0;
  double zz = fabs(za[0]);
  int it;
  for(it=1;it<nd;it++){
    if(fabs(za[it]) < zz){
      il = it;
      zz = fabs(za[it]);
    }
  }
  return il;
}

int maxabsvalloc_d(double* za, int nd){
  int il = 0;
  double zz = 0.0;
  int it;
  for(it=0;it<nd;it++){
    if(fabs(za[it]) > zz){
      il = it;
      zz = fabs(za[it]);
    }
  }
  return il;
}

double entry_ij(int i, int j){
  int il;
  int n[3];
  double xf[3], yf[3], zf[3];
  double xp, yp, zp;

  xp = zgmid[i][0];
  yp = zgmid[i][1];
  zp = zgmid[i][2];
  
  for(il=0;il<3;il++){
    n[il] = f2n[j][il];
  }
  for(il=0;il<3;il++){
    xf[il] = bgmid[n[il]][0];
    yf[il] = bgmid[n[il]][1];
    zf[il] = bgmid[n[il]][2];
  }

  double result = face_integral2(xf, yf, zf, xp, yp, zp);

  return result;
}

double face_integral2(double xs[], double ys[], double zs[], double x, double y, double z){
  int il;
  double PI = 3.1415927;
  double EPSILON_0 = 8.854188 * 1e-12;

  double r[3];
  double xi, xj, yi, dx, dy, t, l, m, d, ti, tj;
  double theta, omega, q, g, zp, zpabs;

  int i, j;
  double *u, *v, *w;
  double ox, oy, oz;

  for(il=0;il<3;il++){
    r[il] = sqrt( pow((xs[il] - x),2.0) + pow((ys[il] - y),2.0) + pow((zs[il] - z),2.0) );
  }

  u = (double*)malloc(sizeof(double)*3);
  v = (double*)malloc(sizeof(double)*3);
  w = (double*)malloc(sizeof(double)*3);

  u[0] = xs[1] - xs[0];
  u[1] = ys[1] - ys[0];
  u[2] = zs[1] - zs[0];

  v[0] = xs[2] - xs[1];
  v[1] = ys[2] - ys[1];
  v[2] = zs[2] - zs[1];
  
  cross_product(u, v, w);
  
  double dw = sqrt( dot_product(w,w,3));
  for(il=0;il<3;il++){
    w[il] = w[il] / dw;
  }
  u[0] = x - xs[0];
  u[1] = y - ys[0];
  u[2] = z - zs[0];
  zp = dot_product(u,w,3);
  
  ox = x - zp * w[0];
  oy = y - zp * w[1];
  oz = z - zp * w[2];
  zpabs = fabs(zp);

  double face_integral = 0.0;
  for(i=0;i<3;i++){
    j = (i + 1) % 3;
    u[0] = xs[j] - ox;
    u[1] = ys[j] - oy;
    u[2] = zs[j] - oz;
    xj = sqrt( dot_product(u,u,3) );

    for(il=0;il<3;il++){
      u[il] = u[il] / xj;
    }
    cross_product(w, u, v);
    xi = (xs[i] - ox) * u[0] + (ys[i] - oy) * u[1] + (zs[i] - oz) * u[2];
    yi = (xs[i] - ox) * v[0] + (ys[i] - oy) * v[1] + (zs[i] - oz) * v[2];

    dx = xj - xi;
    dy = - yi;
    t = sqrt ((dx*dx) + (dy*dy));
    l = dx / t;
    m = dy / t;
    d = (l * yi) - (m * xi);
    ti = (l * xi) + (m * yi);
    tj = l * xj;

    theta = atan2(yi, xi);
    omega = theta - atan2( r[i] * d, zpabs * ti ) + atan2( r[j] * d, zpabs * tj );
    q = log( (r[j] + tj) / (r[i] + ti) );
    g = d * q - zpabs * omega;
    face_integral = face_integral + g;
  }
  
  return fabs(face_integral) / (4.0 * PI * EPSILON_0);

}

void cross_product(double* u, double* v, double* w){
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - u[0] * v[2];
  w[2] = u[0] * v[1] - u[1] * v[0];
}

void adotsub_dsm(double* zr, double* zaa, double* zab, int it, int ndl, int ndt, int mdl, int mdt){
  int il;
  double* zau = (double*)calloc(ndl,sizeof(double));

  adot_dsm(zau,zaa,zab,it,ndl,ndt,mdl,mdt);
  for(il=0;il<ndl;il++){
    zr[il] = zr[il] - zau[il];
  }
  free(zau);
}

void adot_dsm(double* zau, double* zaa, double* zab, int im, int ndl, int ndt, int mdl, int mdt){
  int it,il;
  double (*zaa2)[mdl];
  double (*zab2)[mdt];
  zaa2 = (double(*)[mdl])zaa;
  zab2 = (double(*)[mdt])zab;
  for(it=0;it<ndt;it++){
    for(il=0;il<ndl;il++){
      zau[il] = zau[il] + zaa2[it][il] * zab2[it][im];
    }
  }
}

double dot_product(double* v, double* u, int n){
  double result = 0.0;
  int i;
  for (i = 0; i < n; i++){
    result += v[i] * u[i];
  }
  return result;
}

int max(int a, int b){
  if(a >= b){
    return a;
  }
  return b;
}

int min(int a, int b){
  if(a <= b){
    return a;
  }
  return b;
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
