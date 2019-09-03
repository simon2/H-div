#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <time.h>
#include <sys/time.h>

#define LENGTH 25

void prefix_sum_seq1(int* input, int length);
void prefix_sum_par1(int* input, int length, int init);
int* prefix_sum_seq2(int* input, int length);
void prefix_sum_par2(int* input, int length);
double get_wall_time();

int main(int argc, char **argv){
  printf("begin!\n");
  if(argc >= 2){
    __cilkrts_set_param("nworkers",argv[1]);
  }
  int nworkers = __cilkrts_get_nworkers();
  printf("number of workers:%d\n",nworkers);
  int i;
  int *input = (int*)malloc(LENGTH * sizeof(int));
  double start,end,spent;
  srand(time(NULL));
  for(i=0;i<LENGTH;i++){
    input[i] = rand()%10;
  }
  for(i=0;i<LENGTH;i++){
    printf("%ld ",input[i]);
  }
  printf("\n");
  start = get_wall_time();
#ifdef A
  prefix_sum_par2(input, LENGTH);
#endif
#ifdef B
  prefix_sum_seq1(input, LENGTH);
#endif
#ifdef C
  prefix_sum_par1(input, LENGTH, 0);
#endif
#ifdef D
  prefix_sum_seq2(input, LENGTH);
#endif
  end = get_wall_time();
  spent = end - start;
  printf("time:%.10f\n",spent);

  for(i=0;i<LENGTH;i++){
    printf("%ld ",input[i]);
  }
  printf("\nend!\n");
}

void prefix_sum_seq1(int* input, int length){
  int i, tl = 0, tl2;
  for(i=0;i<length;i++){
    tl2 = input[i];
    input[i] = tl;
    tl = tl + tl2;
  }
}

int* prefix_sum_seq2(int* input, int length){
  int i;
  int* output = (int*)malloc(length*sizeof(int));
  output[0] = 0;
  for(i=1;i<length;i++){
    output[i] = input[i-1] + output[i-1];
  }
  return output;
}

void prefix_sum_par1(int* input, int length, int init){
  int i,j,t,m;
  int k = 0;
  int tl = length;
  while (tl >>= 1) k++;
  int l = 1 << k;
  //printf("k:%ld l:%ld ",k,l);
  for(i=0;i<k;i++){
    t = 1 << (i+1);
    m = 1 << i;
    cilk_for(j=0;j<l;j+=t){
      input[j+t-1] = input[j+m-1] + input[j+t-1];
    }
  }
  int sum = input[l-1] + init;
  //printf("sum:%ld\n",sum);
  input[l-1] = init;
  for(i=k-1;i>=0;i--){
    t = 1 << (i+1);
    m = 1 << i;
    cilk_for(j=0;j<l;j+=t){
      int temp = input[j+m-1];
      input[j+m-1] = input[j+t-1];
      input[j+t-1] = temp + input[j+t-1];
    }
  }

  if(length-l>1){
    prefix_sum_par1(&input[l],length-l,sum);
  }else if(length-l == 1){
    input[l] = sum;
  }
}

void prefix_sum_par2(int* input, int length){
  int i,j;
  int chunk = 2;
  int chunk_len = length/chunk;
  int* sl = (int*)malloc(chunk*sizeof(int));
  cilk_for(i=0;i<chunk;i++){
    int start = i * chunk_len;
    int len = i==chunk-1?length-start:chunk_len;
    int sum = input[start+len-1];
    int tl = 0, tl2;
    for(j=start;j<start+len;j++){
      tl2 = input[j];
      input[j] = tl;
      tl = tl + tl2;
    }
    sum += input[start+len-1];
    sl[i] = sum;
  }
  prefix_sum_seq1(sl,chunk);
  cilk_for(i=0;i<chunk;i++){
    int start = i * chunk_len;
    int end = i==chunk-1?length:start+chunk_len;
    for(j=start;j<end;j++){
      input[j] += sl[i];
    }
  }
}

double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
