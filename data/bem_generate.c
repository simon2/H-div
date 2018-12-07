// bem_generate:
// Generate a new BEM input file by putting copies of the object of a base input file.
// Programmed by Tasuku HIRAISHI <tasuku@media.kyoto-u.ac.jp> in 2018
// Original program is written by Yuki Noseda in 2013 as "make_input.c"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<stdint.h>
#include<inttypes.h>
#include<float.h>
#include<assert.h>
#include<sys/time.h>
#include "bem_file.h"

#define USAGE_STRING \
  "Usage: %s <command> <base_input_file> [<output_file>]\n" \
  "\n" \
  "<command> is one of:\n" \
  "\"p <distance> <n_layer>\"\n" \
  "-> Places 1+4+9+...+<n_layer>^2 copies of the object to form a pyramid\n" \
  "\"c <distance> <nx> <ny> <nz>\"\n" \
  "-> Places <nx>*<ny>*<nz> copies of the object to form a cuboid\n" \
  "\n" \
  "<distance> is specified as the ratio to the diameter of the object.\n" \
  "<base_input_file> is a BEM input file that describes the object to be copied.\n" \
  "When <output_file> is omitted, the output file name will be automatically generated based on <command>.\n"


void usage (int argc, char** argv) 
{
  fprintf (stderr, USAGE_STRING, argv[0]);
  exit (99);
}

// Malloc and returns string = ifile+"_"+argv[i1]+"_"+...+"_"+argv[i2-1]
char* make_ofile_name (char *ifile, char** argv, int i1, int i2)
{
  size_t len = strlen (ifile);
  for (int i=i1; i<i2; i++) {
    len += strlen (argv[i]);
  }
  size_t sz = len + (i2 - i1) + 1; // (Number of '_') + \0
  char *ofile = malloc (sizeof(char) * sz);
  char *p = ofile;
  p = stpcpy (p, ifile);
  for (int i=i1; i<i2; i++) {
    p = stpcpy (p, "_");
    p = stpcpy (p, argv[i]);
  }
  return ofile;
}

// Call usage() if argc < n
void check_argc (int n, int argc, char **argv)
{
  if (argc < n) {
    fprintf (stderr, "Error: insufficient number of arguments.\n\n");
    usage (argc, argv);
  }
}


/* Initialize new_bi based on bi */
/* ncopies: estimated copy of the bi object to be putted */   
void initialize_bem_input (struct bem_input *bi, struct bem_input *new_bi, int ncopies)
{
  int64_t est_nNode = ncopies * bi->nNode;
  int64_t est_nFace = ncopies * bi->nFace;
  fprintf (stderr, 
	   "# of object copies = %d\n"
	   "Estimated # of nodes = %"PRId64"\n" "Estimated # of faces = %"PRId64"\n",
	   ncopies, est_nNode, est_nFace);
  new_bi->nNode = 0;
  new_bi->coordOfNode = (double(*)[3]) malloc (est_nNode * 3 * sizeof(double));
  new_bi->nFace = 0;
  new_bi->nNodePerFace = bi->nNodePerFace;
  new_bi->nIFValue = bi->nIFValue;
  new_bi->nDFValue = bi->nDFValue;
  new_bi->idOfFace = (int64_t*) malloc (est_nFace * bi->nNodePerFace * sizeof(int64_t));
  new_bi->coordOfFace = NULL;   // Face's centers are not treaded in bem_generate
  new_bi->IFValue = (int64_t*) malloc (est_nFace * sizeof(int64_t)* bi->nIFValue);
  new_bi->DFValue = (double*) malloc (est_nFace * sizeof(double)* bi->nDFValue);
}


double diameter (struct bem_input *bi)
{
  double x_min, x_max;
  double y_min, y_max;
  double z_min, z_max;
  x_min = x_max = bi->coordOfNode[0][0];
  y_min = y_max = bi->coordOfNode[0][1];
  z_min = z_max = bi->coordOfNode[0][2];
  for (int64_t i=1; i < bi->nNode; i++) {
    double x = bi->coordOfNode[i][0];
    double y = bi->coordOfNode[i][1];
    double z = bi->coordOfNode[i][2];
    if (x > x_max) x_max = x;
    if (x < x_min) x_min = x;
    if (y > y_max) y_max = y;
    if (y < y_min) y_min = y;
    if (z > z_max) z_max = z;
    if (z < z_min) z_min = z;
  }
  double diameter = x_max-x_min;
  if (diameter < y_max-y_min) diameter = y_max-y_min;
  if (diameter < z_max-z_min) diameter = z_max-z_min;
  fprintf (stderr, "diameter: %lf\n", diameter);
  return diameter;
}


// Add a copy of object of bi to new_bi translating by (dx, dy, dz)
// Space for arrays should have been malloc-ed before calling this function
void add_object_copy (struct bem_input *bi, struct bem_input *new_bi, double dx, double dy, double dz)
{
  int64_t nodeid_offset = new_bi->nNode;
  int64_t faceid_offset = new_bi->nFace;
  //fprintf (stderr, "add_object_copy: nodeid_offset: %"PRId64", faceid_offset: %"PRId64"\n",
  //         nodeid_offset, faceid_offset);
  // Copy nodes
  new_bi->nNode += bi->nNode;
  for (int64_t i=0; i < bi->nNode; i++) {
    const int64_t i1 = i + nodeid_offset;
    new_bi->coordOfNode[i1][0] = bi->coordOfNode[i][0] + dx;
    new_bi->coordOfNode[i1][1] = bi->coordOfNode[i][1] + dy;
    new_bi->coordOfNode[i1][2] = bi->coordOfNode[i][2] + dz;
  }
  // Copy faces
  new_bi->nFace += bi->nFace;
  {
    int64_t npf = bi->nNodePerFace;
    for (int64_t i=0; i < bi->nFace; i++) {
      for (int64_t j=0; j < npf; j++)  {
	const int64_t i1 = i + faceid_offset;
	new_bi->idOfFace[i1*npf+j] = bi->idOfFace[i*npf+j] + nodeid_offset;
      }
    }
  }
  // Copy IFValue and DFValue
  {
    assert (bi->nIFValue == new_bi->nIFValue);
    int64_t nifv = bi->nIFValue;
    memcpy (new_bi->IFValue + faceid_offset * nifv, bi->IFValue,
	    sizeof(int64_t) * bi->nFace * nifv);
  }
  {
    assert (bi->nDFValue == new_bi->nDFValue);
    int64_t ndfv = bi->nDFValue;
    memcpy (new_bi->DFValue + faceid_offset * ndfv, bi->DFValue,
	    sizeof(double) * bi->nFace * ndfv);
  }
  return;
}


void place_pyramid (int layers, double distance, struct bem_input *bi, struct bem_input *new_bi)
{
  double dx, dy, dz;
  double diam = diameter (bi);
  double offset = diam * distance;
  dz = -((double)(layers)*0.5*offset);
  for (int iz=1; iz<=layers; iz++) {
    double dxy0 = - ((double)(iz-1) * 0.5 * offset);
    dy = dxy0;
    for (int iy=1; iy<=iz; iy++) {
      dx = dxy0;
      for (int ix=1; ix<=iz; ix++) {
	add_object_copy (bi, new_bi, dx, dy, dz);
	dx += offset;
      }
      dy += offset;
    }
    dz += offset;
  }
  return;
}

void place_cuboid (int nx, int ny, int nz, double distance, struct bem_input *bi, struct bem_input *new_bi)
{
  double dx, dy, dz;
  double diam = diameter (bi);
  double offset = diam * distance;
  double dx0 = - ((double)(nx)*0.5*offset);
  double dy0 = - ((double)(ny)*0.5*offset);
  double dz0 = - ((double)(nz)*0.5*offset);
  dz = dz0;
  for (int iz=0; iz<nz; iz++) {
    dy = dy0;
    for (int iy=0; iy<ny; iy++) {
      dx = dx0;
      for (int ix=0; ix<nx; ix++) {
	add_object_copy (bi, new_bi, dx, dy, dz);
	dx += offset;
      }
      dy += offset;
    }
    dz += offset;
  }
  return;
}


int main( int argc, char *argv[] )
{
  struct timeval tv1, tv2;
  char *ifile, *ofile;
  struct bem_input bi, new_bi;
  enum bi_format fmt;

  gettimeofday (&tv1, NULL);
  check_argc (2, argc, argv);
  char c = argv[1][0];
  switch (c) { 
  case 'p':
    /* Pyramid */
    {
      check_argc (5, argc, argv);
      double distance = atof(argv[2]);
      int64_t layers = atoll(argv[3]);
      ifile = argv[4];
      ofile = (argc>5)?argv[5]:make_ofile_name (ifile, argv, 1, 4);
      // Read base file
      fmt = open_and_read_bem_input (ifile, &bi, BI_AUTO);
      fprintf (stderr, "--- summary of input file -------\n");
      print_bem_input (stderr, &bi, BI_PRETTY);
      fprintf (stderr, "--- summary of input file end ---\n");
      // Calculate # copies and initialize
      int64_t ncopies = (2*layers+1)*layers*(layers+1)/6;
      initialize_bem_input (&bi, &new_bi, ncopies);
      // Generate pyramid
      place_pyramid (layers, distance, &bi, &new_bi);
    }
    break;
  case 'c':
    /* Cuboid */
    {
      check_argc (7, argc, argv);
      double distance = atof(argv[2]);
      int nx = atoi(argv[3]);
      int ny = atoi(argv[4]);
      int nz = atoi(argv[5]);
      FILE *fpin;
      ifile = argv[6];
      ofile = (argc>7)?argv[7]:make_ofile_name (ifile, argv, 1, 6);
      // Read
      fmt = open_and_read_bem_input (ifile, &bi, BI_AUTO);
      fprintf (stderr, "--- summary of input file -------\n");
      print_bem_input (stderr, &bi, BI_PRETTY);
      fprintf (stderr, "--- summary of input file end ---\n");
      // Calculate # copies and initialize
      int64_t ncopies = nx * ny * nz;
      initialize_bem_input (&bi, &new_bi, ncopies);
      // Generate cuboid
      place_cuboid (nx, ny, nz, distance, &bi, &new_bi);
    }
    break;
  default:
    fprintf (stderr, "Error: illegal command: %c\n\n", c);
    usage (argc, argv);
  }

  // When the second character of the first argument is 'b' or 't',
  // Write the output file in the Binary or Text format respectively.
  // Otherwise, written in the same format to the input file.
  if (argv[1][1] == 'b') fmt = BI_BINARY;
  if (argv[1][1] == 't') fmt = BI_TEXT;
  open_write_show (1, &new_bi, fmt, ofile);

  // Show summary
  fprintf (stderr, "--- summary of output file -------\n");
  print_bem_input (stdout, &new_bi, BI_PRETTY);
  fprintf (stderr, "--- summary of output file end ---\n");

  return 0;
}
