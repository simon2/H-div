#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "data/bem_file.h"

#ifdef DEBUG
#include <assert.h>
#endif

#define INPUT_DEFAULT "bem_data/input_338ts.txt"

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
void qsort_col_leafmtx(leafmtx *st_leafmtx,int first,int last);
void qsort_row_leafmtx(leafmtx *st_leafmtx,int first,int last);
void create_leafmtx(leafmtx *st_leafmtx,int st_cltl,int st_cltt,
                    double param[],int *lnmtx,int nffc,int *nlf,double (*gmid)[3],double (*bmid)[3],int (*face2node)[3]);
double dist_2cluster(int st_cltl,int st_cltt);
void count_lntmx(int st_cltl,int st_cltt,double param[],int *lnmtx,int nffc);
int create_cluster(int nmbr,int ndpth,int nstrt,int nsize,int ndim,int nson);
int create_ctree_ssgeom(int st_clt,double (*zgmid)[3],int (*face2node)[3],double param[],int ndpth,int ndscd,int nsrt,int nd,int md,int ndim,int nclst);

double get_wall_time();
double get_cpu_time();

void checkClusterTree(int st_clt);
void checkClusterTree2(int ndpth, int st_clt);
int depth_max;
int count_node;

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
  param[21] = 50.0;
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
  double start = get_wall_time();
  supermatrix_construction_cog_leafmtrx(st_leafmtxp,coordOfFace,coordOfNode,face2node,param,lod,lnmtx,nofc,nffc,ndim);  // construction of Leaf-matrix 
  double end = get_wall_time();
  printf("MP time:%f\n", (end - start));
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
  start = get_wall_time();
#ifdef DEBUG
  lod0 = lodfc;
#endif 

  CTlist = (cluster*)malloc(sizeof(cluster) * LN);
  st_clt = create_ctree_ssgeom(st_clt,gmid,face2node,param,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst);
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
  start = get_wall_time();
  create_leafmtx(st_leafmtx,st_clt,st_clt,param,lnmtx,nffc,&nlf,gmid, bmid, face2node);
  end = get_wall_time();
  spent = end - start;
  printf("nlf:%d\n",nlf);
  printf("block cluster tree time spent:%.10f\n",spent);
  printf("depth_max:%d  count_node:%d\n",depth_max,count_node);
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
                    double param[],int *lnmtx,int nffc,int *nlf, double (*face)[3], double (*node)[3], int (*face2node)[3]){
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
          create_leafmtx(st_leafmtx,CTlist[st_cltl].offsets[il],CTlist[st_cltt].offsets[it],param,lnmtx,nffc,nlf,face, node, face2node);
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
