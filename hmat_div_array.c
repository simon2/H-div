#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
// #include <string.h>
#include <mkl.h>
#include "data/bem_file.h"

#ifdef DEBUG
#include <assert.h>
#endif

#define INPUT_DEFAULT "bem_data/input_10ts.txt"

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
  int nlf;                           //number of partitions
  int nlfkt;                         //number ot partitions approximated
};

void supermatrix_construction_cog_leafmtrx(leafmtxp *st_leafmtxp,double (*gmid)[3],double (*bmid)[3],int (*face2node)[3],
                                            double param[],int *lod,int *lnmtx,int nofc,int nffc,int ndim);
int med3(int nl,int nr,int nlr2);
void create_leafmtx(leafmtx *st_leafmtx,int st_cltl,int st_cltt,
                    double param[],int *lnmtx,int nffc,int *nlf);
double dist_2cluster(int st_cltl,int st_cltt);
void count_lntmx(int st_cltl,int st_cltt,double param[],int *lnmtx,int nffc);
int create_cluster(int nmbr,int ndpth,int nstrt,int nsize,int ndim,int nson);
//void free_st_clt(int st_clt);
int create_ctree_ssgeom(int st_clt,double (*zgmid)[3],int (*face2node)[3],double param[],int ndpth,int ndscd,int nsrt,int nd,int md,int ndim,int nclst);

int acaplus(double* zaa, double* zab, int ndl, int ndt, int nstrtl, int nstrtt, double (*face)[3], double (*node)[3], int (*face2node)[3], int kmax, double eps, double znrmmat, double pACA_EPS);
void fill_leafmtx(leafmtx *st_lf, double (*face)[3], double (*node)[3], int (*face2node)[3], double znrmmat, int *lnmtx, int nd, int nlf);
int minabsvalloc_d(double* za, int nd);
int maxabsvalloc_d(double* za, int nd);
double entry_ij(int i, int j, double (*face)[3], double (*node)[3], int (*face2node)[3]);
double face_integral2(double xs[], double ys[], double zs[], double x, double y, double z);
void comp_col(double* zaa, double *zab, int ndl, int ndt, int k, int it, double* col, int nstrtl, int nstrtt, double (*face)[3], double (*node)[3], int (*face2node)[3], int* lrow_done);
void comp_row(double* zaa, double* zab, int ndl, int ndt, int k, int il, double* row, int nstrtl, int nstrtt, double (*face)[3], double (*node)[3], int (*face2node)[3], int* lrow_done);
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
// int sumn[100];
// int sums[100];
// int max[100];

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
    coordOfNode = bi.coordOfNode;
    count = bi.nFace;
    coordOfFace = bi.coordOfFace;
    face2node = bi.face2node;
  }
  fclose(file);
  // free(coordOfNode);
  double param[100];
  for(i=0;i<100;i++){
    param[i] = 0;
  }
  param[21] = 10.0;
  param[31] = 1.1;
  param[41] = 15.0;
  param[51] = 2.0;

  // for(i=0;i<count;i++){
  //   printf("face2node:%d %d %d\n",face2node[i][0],face2node[i][1],face2node[i][2]);
  // }

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

  supermatrix_construction_cog_leafmtrx(st_leafmtxp,coordOfFace,coordOfNode,face2node,param,lod,lnmtx,nofc,nffc,ndim);  // construction of Leaf-matrix 

//void fill_leafmtx(leafmtx *st_lf, double (*face)[3], double (*node)[3], int (*face2node)[3], double znrmmat, int *lpmd, int *lnmtx, int nd, int nlf/*, int lnps, int lnpe*/);
  return 0;
}

#ifdef DEBUG
int *lod0;
#endif

void supermatrix_construction_cog_leafmtrx(leafmtxp *st_leafmtxp,    //the H-matrix
					   double (*gmid)[3],            //coordination of objects
             double (*bmid)[3],
             int (*face2node)[3],
					   double param[],int *lod,
					   int *lnmtx,              //1:k-rank 2:dense 3:H-matrix
					   int nofc,int nffc,       //number of elements in same coordination
					   int ndim){
  //cluster *st_clt = (cluster *)malloc(sizeof(cluster));
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
  
  // for(il=0;il<100;il++){
  //   sumn[il] = 0;
  //   sums[il] = 0;
  //   max[il] = 0;
  // }

  double start,end,spent;
  start = get_wall_time();
#ifdef DEBUG
  lod0 = lodfc;
#endif

// int csf = 0;
//   for(il=0;il<ndf;il++){
//     for(ip=0;ip<ndim;ip++){
//       double n = 0.0;
//       for(ig=0;ig<ndim;ig++){
//         n = n + bmid[face2node[il][ig]][ip];
//       }
//       n = n / 3.0;
//       printf("%f ",n - gmid[il][ip] );
//     }
//   }
  

  CTlist = (cluster*)malloc(sizeof(cluster) * LN);
  //printf("cluster tree time spent:%.10f\n",spent);
  st_clt = create_ctree_ssgeom(st_clt,gmid,face2node,param,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst);
  end = get_wall_time();
  spent = end - start;
  printf("cluster tree time spent:%.10f\n",spent);
  printf("breakdown:%.10f\n",middle-start);

  // csf = 0;
  // for(il=0;il<ndf;il++){
  //   for(ip=0;ip<ndim;ip++){
  //     double n = 0.0;
  //     for(ig=0;ig<ndim;ig++){
  //       n = n + bmid[face2node[il][ig]][ip];
  //     }
  //     n = n / 3.0;
  //     printf("%f ",n - gmid[il][ip] );
  //   }
  // }


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

  start = get_wall_time();
  fill_leafmtx(st_leafmtx, gmid, bmid, face2node, 0.0, lnmtx, ndf, nlf);
  end = get_wall_time();
  spent = end - start;
  printf("filling time:%.10f\n",spent);
  // printf("nlf:%d\n",nlf);
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
  if((CTlist[st_cltl].zwdth * zeta <= zdistlt || CTlist[st_cltt].zwdth * zeta <= zdistlt) && (ndl >= nleaf && ndt >= nleaf)){
    //st_leafmtx[*nlf] = (leafmtx *)malloc(sizeof(leafmtx));
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
      //st_leafmtx[*nlf]= (leafmtx *)malloc(sizeof(leafmtx));
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
          //create_leafmtx(st_leafmtx,st_cltl->pc_sons[il],st_cltt->pc_sons[it],param,lnmtx,nffc,nlf);
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
  if ((CTlist[st_cltl].zwdth * zeta <= zdistlt || CTlist[st_cltt].zwdth * zeta <= zdistlt) && (ndl >= nleaf && ndt >= nleaf)){
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
          int nh = face2node[nl][id];
          face2node[nl][id] = face2node[nr][id];
          face2node[nr][id] = nh;
        }
      }
    }

#ifdef DEBUG
    fprintf (stdout, "nl = %ld, nr = %ld\n", nl, nr);
    assert ( nl-nr == 1 );
    // assert ( nl < nd );
#endif
    
    if(nl == nd || nl == 0){
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
    zwdth = (CTlist[st_clt].bmax[id] - CTlist[st_clt].bmin[id]) * (CTlist[st_clt].bmax[id] - CTlist[st_clt].bmin[id]);
  }
  CTlist[st_clt].zwdth = sqrt(zwdth);
  //end of bounding box
  return st_clt;
}

int acaplus(double* zaa, double* zab, int ndl, int ndt, int nstrtl, int nstrtt, double (*face)[3], double (*node)[3], int (*face2node)[3], int kmax, double eps, double znrmmat, double pACA_EPS){
  double *prow, *pcol, *pb_ref, *pa_ref;
  int *lrow_done, *lcol_done;

  // printf("ndl=%d ndt=%d\n", ndl, ndt); 

  double (*zaa2)[ndl];
  double (*zab2)[ndt];
  zaa2 = (double(*)[ndl])zaa;
  zab2 = (double(*)[ndt])zab;
  
  int il,it,ib;

  int INCY = 1;
  double za_ACA_EPS = 1.0e-10;

  double znrm = znrmmat * sqrt((double)ndl * (double)ndt);
  double ACA_EPS = pACA_EPS;
  printf("ACA_EPS=%f\n",ACA_EPS);

  int ntries = max(ndl, ndt) + 1; //why +1 here?
  int ntries_row = 6;
  int ntries_col = 6;

  lrow_done = (int *)calloc(ndl, sizeof(int));
  lcol_done = (int *)calloc(ndt, sizeof(int));
  int k = 0;

  int j_ref = 0;
  pa_ref = (double *)malloc(ndl * sizeof(double));

  comp_col(zaa, zab, ndl, ndt, k, j_ref, pa_ref, nstrtl, nstrtt, face, node, face2node, lrow_done);
  // for(ib=0;ib<ndl;ib++){
  //   printf("%f ",pa_ref[ib]);
  // }
  // printf("\n");
  double colnorm = cblas_dnrm2(ndl, pa_ref, INCY);

  int i_ref = minabsvalloc_d(pa_ref, ndl);
  double rownorm = fabs(pa_ref[i_ref]);
  pb_ref = (double *)malloc(ndt * sizeof(double));
  comp_row(zaa, zab, ndl, ndt, k, i_ref, pb_ref, nstrtl, nstrtt, face, node, face2node, lcol_done);
  // for(ib=0;ib<ndl;ib++){
  //   printf("%f ",pb_ref[ib]);
  // }
  // printf("\n");
  rownorm = cblas_dnrm2(ndt, pb_ref, INCY);

  // printf("colnorm=%f rownorm=%f\n", colnorm, rownorm);
  // printf("i_ref=%d j_ref=%d\n", i_ref, j_ref);

  double apxnorm = 0.0;
  int lstop_aca = 0;
  
  double col_maxval, row_maxval;
  while(k < kmax && (ntries_row > 0 || ntries_col > 0) && ntries > 0){
    printf("ntries:%d, ntries_row:%d, ntries_col:%d\n",ntries,ntries_row,ntries_col);
    ntries--;
    double *pcol = &zaa2[k][0]; //Fortran: pcol => zaa(:,k+1) zaa and zab are defined as target in Fortran
    double *prow = &zab2[k][0]; //Fortran: prow => zab(:,k+1)
    col_maxval = 0.0;
    int i = maxabsvalloc_d(pa_ref, ndl);
    col_maxval = fabs(pa_ref[i]);
    row_maxval = 0.0;
    int j = maxabsvalloc_d(pb_ref, ndt);
    row_maxval = fabs(pb_ref[j]);
    // printf("col_maxval=%f row_maxval=%f\n", col_maxval, row_maxval);
    // printf("i=%d j=%d\n", i, j);

    double zinvmax; //no definiatino in Fortran
    if(row_maxval > col_maxval){
      // printf("row_maxval > col_maxval\n");
      if(j != j_ref){
        // printf("ndl=%d, ndt=%d, k=%d, j=%d, nstrtl=%d, nstrtt=%d\n",ndl,ndt,k,j,nstrtl,nstrtt);
        comp_col(zaa, zab, ndl, ndt, k, j, pcol, nstrtl, nstrtt, face, node, face2node, lrow_done);
      }else{
        // printf("here2\n");
        for(il=0;il<ndl;il++){
          pcol[il] = pa_ref[il];
        }
      }
      // printf("here3\n");
      i = maxabsvalloc_d(pcol, ndl);
      col_maxval = fabs(pcol[i]);
      // printf("here4\n");
      // printf("col_maxval=%f row_maxval=%f\n", col_maxval, row_maxval);

      if(col_maxval < ACA_EPS && k >= 1){
        printf("here1\n");
        lstop_aca = 1;
      }else{
        // printf("here6\n");
        // printf("ndl=%d, ndt=%d, k=%d, i=%d, nstrtl=%d, nstrtt=%d\n",ndl,ndt,k,i,nstrtl,nstrtt);
        comp_row(zaa, zab, ndl, ndt, k, i, prow, nstrtl, nstrtt, face, node, face2node, lcol_done);
        if(fabs(pcol[i]) > 1.0e-20){
          // printf("here7\n");
          zinvmax = 1.0 / pcol[i];
        }else{
          printf("here2!\n");
          k = max(k-1, 0);
          break;
        }
        // printf("here9\n");
        for(il=0;il<ndl;il++){
          pcol[il] *= zinvmax;
        } 
      }
    }else{
      // printf("col_maxval > row_maxval\n");
      if(i != i_ref){
        comp_row(zaa, zab, ndl, ndt, k, i, prow, nstrtl, nstrtt, face, node, face2node, lcol_done);
      }else{
        for(il=0;il<ndt;il++){
          prow[il] = pb_ref[il];
        }
      }
      // printf("prow1:");
      // for(ib=0;ib<ndt;ib++){
      //   printf("%d:%f ",ib,prow[ib]);
      // }
      // printf("\n");
      j = maxabsvalloc_d(prow, ndt);
      row_maxval = fabs(prow[j]);
      // printf("prow2:");
      // for(ib=0;ib<ndt;ib++){
      //   printf("%d:%f ",ib,prow[ib]);
      // }
      // printf("\n");
      // printf("col_maxval=%f row_maxval=%f\n", col_maxval, row_maxval);

      if(row_maxval < ACA_EPS && k >= 1){
        printf("here3\n");
        lstop_aca = 1;
      }else{
        comp_col(zaa, zab, ndl, ndt, k, j, pcol, nstrtl, nstrtt, face, node, face2node, lrow_done);
        // printf("fabs(prow[%d])=%f\n",j,fabs(prow[j]));
        if(fabs(prow[j]) > 1.0e-20){
          zinvmax = 1.0 / prow[j];
        }else{
          printf("here4!\n");
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
    // printf("get here!\n");

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
            comp_row(zaa, zab, ndl, ndt, k+1, i, pb_ref, nstrtl, nstrtt, face, node, face2node, lcol_done);
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
            comp_col(zaa, zab, ndl, ndt, k+1, j, pa_ref, nstrtl, nstrtt, face, node, face2node, lrow_done);
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
    // printf("colnorm=%f rownorm=%f\n", colnorm, rownorm);
    // printf("i_ref=%d j_ref=%d\n", i_ref, j_ref);

    if(colnorm < ACA_EPS && rownorm < ACA_EPS && k >= 1){
      printf("here5\n");
      lstop_aca = 1;
      k = k + 1;
    }

    if(lstop_aca == 0){
      double blknorm = cblas_dnrm2(ndl, pcol, INCY) * cblas_dnrm2(ndt, prow, INCY);
      printf("blknorm:%f, apxnorm:%f, eps:%f, rownorm:%f, colnorm:%f, k:%d\n", blknorm, apxnorm, eps, rownorm, colnorm, k);
      if(k == 0){
        printf("here6\n");
        apxnorm = blknorm;
      }else{
        double compared = apxnorm * eps;
        if(blknorm < compared && rownorm < compared && colnorm < compared && k >= 1){
          // printf("blknorm:%f, apxnorm:%f, eps:%f, rownorm:%f, colnorm:%f, k:%d\n", blknorm, apxnorm, eps, rownorm, colnorm, k);
          lstop_aca = 1;
        }
      }
    }
    if(lstop_aca == 1 && k >= 1){
      printf("here7\n");
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

void fill_leafmtx(leafmtx *st_lf, double (*face)[3], double (*node)[3], int (*face2node)[3], double znrmmat, int *lnmtx, int nd, int nlf){
  double eps = 1.0e-4;
  double ACA_EPS = 0.0 * eps;
  int kparam = 200;
  int ip,il,it;
  int real_max_k = 0; 

  for(ip=0;ip<nlf;ip++){
    int ndl = st_lf[ip].ndl;
    int ndt = st_lf[ip].ndt;
    int ns = ndl * ndt;
    int nstrtl = st_lf[ip].nstrtl;
    int nstrtt = st_lf[ip].nstrtt;
    int ltmtx = st_lf[ip].ltmtx;
    

    if(ltmtx == 1){
      // printf("low-rank mtx!\n");
      double * zaa = (double*)malloc(sizeof(double) * ndt * kparam);
      double * zab = (double*)malloc(sizeof(double) * ndl * kparam);
      // st_lf[ip].a1 = (double*)malloc(sizeof(double) * ndt * kparam);
      // st_lf[ip].a2 = (double*)malloc(sizeof(double) * ndl * kparam);
      // if(!st_lf[ip].a1 && !st_lf[ip].a2){
      //   printf("allocate a1 or a2 failed!\n");
      //   exit(99);
      // }

      // int kt = acaplus(st_lf[ip].a2, st_lf[ip].a1, ndl, ndt, nstrtl, nstrtt, face, node, face2node, kparam, eps, znrmmat, ACA_EPS);
      int kt = acaplus(zab, zaa, ndl, ndt, nstrtl, nstrtt, face, node, face2node, kparam, eps, znrmmat, ACA_EPS);
      // printf("DEBUGING: kt=%d, kparam=%d, nstrtl=%d, nstrtt=%d, ndl=%d, ndt=%d\n", kt, kparam, nstrtl, nstrtt, ndl, ndt);

      if(kt > kparam){ //Fortran: kt > kparam-1. kt is the rank.
        printf("WARNING: Insufficient k: kt=%d, kparam=%d, nstrtl=%d, nstrtt=%d, ndl=%d, ndt=%d\n", kt, kparam, nstrtl, nstrtt, ndl, ndt);
      }

      if(kt > real_max_k) real_max_k = kt;

      // st_lf[ip].a1 = (double *)realloc(st_lf[ip].a1, kt * ndt * sizeof(double));//check if realloc is right
      // st_lf[ip].a2 = (double *)realloc(st_lf[ip].a2, kt * ndl * sizeof(double));

      st_lf[ip].a1 = (double*)malloc(sizeof(double) * ndt * kt);
      st_lf[ip].a2 = (double*)malloc(sizeof(double) * ndl * kt);

      for(il=0;il<kt;il++){
        for(it=0;it<ndt;it++){
          st_lf[ip].a1[il*ndt+it] = zaa[il*ndt+it];
        }
      }
      for(il=0;il<kt;il++){
        for(it=0;it<ndl;it++){
          st_lf[ip].a2[il*ndl+it] = zab[il*ndl+it];
        }
      }
      free(zaa);free(zab);
      break;

    }else if(st_lf[ip].ltmtx == 2){
      // printf("dense mtx!\n");
      // printf("DEBUGING:  nstrtl=%d, nstrtt=%d, ndl=%d, ndt=%d\n", nstrtl, nstrtt, ndl, ndt);
      st_lf[ip].a1 = (double *)malloc(sizeof(double) * ns);

      double (*tempa1)[ndt];
      tempa1 = (double(*)[ndt])st_lf[ip].a1;
  
      for(il=0;il<ndl;il++){
        int ill = il + nstrtl;
        for(it=0;it<ndt;it++){
          int itt = it + nstrtt;
          tempa1[il][it] = entry_ij(ill, itt, face, node, face2node);
        }
      }
    }
    // printf("ip=%d\n",ip);
  }
  printf("max_kt:%d\n",real_max_k);
}

void comp_row(double* zaa, double* zab, int ndl, int ndt, int k, int il, double* row, int nstrtl, int nstrtt, double (*face)[3], double (*node)[3], int (*face2node)[3], int* lrow_done){
  int it;
  for(it=0;it<ndt;it++){
    if(lrow_done[it] == 0){
      int ill = il + nstrtl;
      int itt = it + nstrtt;
      row[it] = entry_ij(ill, itt, face, node, face2node);
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

void comp_col(double* zaa, double *zab, int ndl, int ndt, int k, int it, double* col, int nstrtl, int nstrtt, double (*face)[3], double (*node)[3], int (*face2node)[3], int* lrow_done){
  int il;

  for(il=0;il<ndl;il++){
    if(lrow_done[il] == 0){
      int ill = il + nstrtl;
      int itt = it + nstrtt;
      col[il] = entry_ij(ill, itt, face, node, face2node);
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

double entry_ij(int i, int j, double (*face)[3], double (*node)[3], int (*face2node)[3]){
  int il;
  int n[3];
  double xf[3], yf[3], zf[3];
  double xp, yp, zp;

  xp = face[i][0];
  yp = face[i][1];
  zp = face[i][2];

  for(il=0;il<3;il++){
    n[il] = face2node[j][il];
  }
  for(il=0;il<3;il++){
    xf[il] = node[n[il]][0];
    yf[il] = node[n[il]][1];
    zf[il] = node[n[il]][2];
  }

  return face_integral2(xf, yf, zf, xp, yp, zp);
}

double face_integral2(double xs[], double ys[], double zs[], double x, double y, double z){
  int il;
  const double PI = 3.141592653589793238462643383279;
  const double EPSILON_0 = 8.854187818 * 1e-12l;

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
    w[il] = w[il] / dw; //it is strange here!
  }

  u[0] = x - xs[0];
  u[1] = y - ys[0];
  u[2] = z - zs[0];
  zp = dot_product(u,w,3);
  ox = x - zp * w[0];
  oy = y - zp * w[1];
  oz = z - zp * w[2];
  zpabs = fabs(zp); //Fortran did not use dabs here.

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
  double* zau = (double*)calloc(ndl,sizeof(double));
  adot_dsm(zau,zaa,zab,it,ndl,ndt,mdl,mdt);
  int il;
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
