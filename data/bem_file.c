// The routine to load a PPOpen BEM input file
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "bem_file.h"

#define BUFSIZE 10000
#define BI_BINARY_PREAMBLE "BI_BINARY"


/* For debugging */
void print_longs (long* vals, long nval) {
  for (int i=0; i<nval; i++) {
    printf ("%ld: %ld\n", i, vals[i]); }}
void print_doubles (double* vals, long nval) {
  for (int i=0; i<nval; i++) {
    printf ("%ld: %.16g\n", i, vals[i]); }}


/* Aux functions called by print_bem_input */
// print 1 long value (val)
void print_bem_input_long (FILE* fp, long val,
			   enum bi_format fmt, const char* prt_string)
{
  switch (fmt) {
  case BI_TEXT:
    fprintf (fp, "%ld\n", val); break;
  case BI_BINARY:
    fwrite (&val, sizeof(long), 1, fp); break;
  case BI_PRETTY:
    fprintf (fp, "%s %ld\n", prt_string ,val); break;
  default:
    fprintf (stderr, "Unknown bi_format in print_bem_input_long!\n");
    return;
  }
}

// print long values (vals)
// nval is # of values. Put a newline after each nsep values in text formsts.
void print_bem_input_longs (FILE* fp, long* vals, long nval, long nsep,
			    enum bi_format fmt, const char* prt_string)
{
  switch (fmt) {
  case BI_TEXT:
    for (long i=0; i<nval; i++) {
      fprintf (fp, "%ld", vals[i]);
      fputc ((i%nsep == nsep-1 || i == nval-1)?'\n':' ', fp);
    }
    break;
  case BI_BINARY:
    fwrite (vals, sizeof(long), nval, fp);
    break;
  case BI_PRETTY:
    fprintf (fp, "%s\n", prt_string);
    print_bem_input_longs (fp, vals, nval, nsep, BI_TEXT, prt_string);
    break;
  default:
    fprintf (stderr, "Unknown bi_format in print_bem_input_longs!\n");
    return;
  }
}

// print double values (vals)
// nval is # of values. Put a newline after each nsep values in text formsts.
void print_bem_input_doubles (FILE* fp, double* vals, long nval, long nsep,
			      enum bi_format fmt, const char* prt_string)
{
  switch (fmt) {
  case BI_TEXT:
    for (long i=0; i<nval; i++) {
      fprintf (fp, "%.16g", vals[i]);
      putc ((i%nsep == nsep-1 || i == nval-1)?'\n':' ', fp);
    }
    break;
  case BI_BINARY:
    fwrite (vals, sizeof(double), nval, fp);
    break;
  case BI_PRETTY:
    fprintf (fp, "%s\n", prt_string);
    print_bem_input_doubles (fp, vals, nval, nsep, BI_TEXT, prt_string);
    break;
  default:
    fprintf (stderr, "Unknown bi_format in print_bem_input_doubles!\n");
    return;
  }
}


/* Print BEM input (pbin) to (fp) in the specified format (fmt) */
void print_bem_input (FILE* fp, struct bem_input* pbin, enum bi_format fmt)
{
  // Check whether fmt is available value.
  switch (fmt) {
  case BI_TEXT:
  case BI_BINARY:
  case BI_PRETTY:
    break;
  default:
    fprintf (stderr, "Unknown bi_format in print_bem_input!\n");
    return;
  }
  // Write preamble for BI_BINARY
  if (fmt == BI_BINARY) {
    fprintf (fp, BI_BINARY_PREAMBLE "\n");
  }
  // Print
  print_bem_input_long (fp, pbin->nNode, fmt, "Number of nodes:");
  print_bem_input_doubles (fp, (double*)pbin->coordOfNode, 3*pbin->nNode, 3,
			   fmt, "Coordinates of the nodes:");
  print_bem_input_long (fp, pbin->nFace, fmt, "Number of faces:");
  print_bem_input_long (fp, pbin->nNodePerFace, fmt, "Number of nodes for each face:");
  print_bem_input_long (fp, pbin->nIFValue, fmt,
			"Number of long values for each face:");
  print_bem_input_long (fp, pbin->nDFValue, fmt,
			"Number of double values for each face:");
  print_bem_input_longs (fp, pbin->idOfFace,
			 pbin->nNodePerFace*pbin->nFace,
			 pbin->nNodePerFace, fmt, "Node IDs forming faces:");
  print_bem_input_longs (fp, pbin->IFValue, pbin->nFace*pbin->nIFValue, 1,
			 fmt, "Long parameter values of faces:");
  print_bem_input_doubles (fp, pbin->DFValue, pbin->nFace*pbin->nDFValue, 1,
			   fmt, "Double parameter values of faces:");
}

 
/* Aux functions called by read_bem_input */
// Read 1 long value and store to (pval)
void read_bem_input_long (FILE* fp, long* pval, enum bi_format fmt)
{
  char line[BUFSIZE];
  switch (fmt) {
  case BI_TEXT: 
    fgets (line, BUFSIZE, fp);
    *pval = atol (line);
    break;
  case BI_BINARY:
    fread (pval, sizeof(long), 1, fp);
    break;
  default:
    fprintf (stderr, "Unknown bi_format in read_bem_input_long!\n");
    return;
  }
}

// Read long values and store to (vals)
// nval is # of values. Gut a newline after each nsep values in text formsts.
void read_bem_input_longs (FILE* fp, long* vals, long nval, long nsep,
				enum bi_format fmt)
{
  char line[BUFSIZE];
  switch (fmt) {
  case BI_TEXT:
    {
      long i=0;
      while (i<nval) {
	fgets (line, BUFSIZE, fp);
	char* pos = line;
	char* pos_nxt;
	for (long j = 0; i<nval && j<nsep; i++,j++) {
	  vals[i] = strtol (pos, &pos_nxt, 10);
	  pos = pos_nxt;
	}
      }
    }
    break;
  case BI_BINARY:
    fread (vals, sizeof(long), nval, fp);
    break;
  default:
    fprintf (stderr, "Unknown bi_format in read_bem_input_longs!\n");
    return;
  }
}


// Read double values and store to (vals)
// nval is # of values. Gut a newline after each nsep values in text formsts.
void read_bem_input_doubles (FILE* fp, double* vals, long nval, long nsep,
				  enum bi_format fmt)
{
  char line[BUFSIZE];
  switch (fmt) {
  case BI_TEXT:
    {
      long i=0;
      while (i<nval) {
	fgets (line, BUFSIZE, fp);
	char* pos = line;
	char* pos_nxt;
	for (long j = 0; i<nval && j<nsep; i++,j++) {
	  vals[i] = strtod (pos, &pos_nxt);
	  pos = pos_nxt;
	}
      }
    }
    break;
  case BI_BINARY:
    fread (vals, sizeof(double), nval, fp);
    break;
  default:
    fprintf (stderr, "Unknown bi_format in read_bem_input_doubles!\n");
    return;
  }
}


// Read BEM input data from fp and store it to the object pointed by pbin
// fp: file pointer of the input
// fmt: file format (BI_TEXT, BI_BINARY, BI_AUTO)
// Returns the format of input file if succeeded, returns -1 if failed.
enum bi_format read_bem_input (FILE* fp, struct bem_input* pbin, enum bi_format fmt)
{
  // Check and read preamble for BI_BINARY
  if (fmt == BI_BINARY || fmt == BI_AUTO) {
    if ( fgetc(fp) != BI_BINARY_PREAMBLE[0] ) {
      if (fmt == BI_BINARY) {
	fprintf (stderr, "This file does not look a BEM input file in the binary format in read_bem_input!\n");
	return -1;
      }
      fmt = BI_TEXT;
    } else {
      ungetc (BI_BINARY_PREAMBLE[0], fp);
      const int bufsz = sizeof(BI_BINARY_PREAMBLE "\n")+1;
      char buf[bufsz]; // BI_BINARY_PREAMBLE + \n + \0
      fgets (buf, bufsz, fp);
      if ( strncmp(buf, BI_BINARY_PREAMBLE "\n", bufsz) ) {
	fprintf (stderr, "This file does not look a BEM input file in the binary format in read_bem_input!\n");
	return -1;
      }
      if (fmt == BI_AUTO) {
	fmt = BI_BINARY;
      }
    }
  }
  // Check whether fmt is available value.
  switch (fmt) {
  case BI_TEXT:
  case BI_BINARY:
    break;
  default:
    fprintf (stderr, "Unknown bi_format in read_bem_input!\n");
    return -1;
  }
  // Get nNode
  read_bem_input_long (fp, &pbin->nNode, fmt);
  // Allocate and get coordOfNode
  pbin->coordOfNode = (double(*)[3]) malloc (pbin->nNode * 3 * sizeof(double));
  read_bem_input_doubles (fp, (double*)pbin->coordOfNode, pbin->nNode * 3, 3, fmt);
  // Get nFace
  read_bem_input_long (fp, &pbin->nFace, fmt);
  // Get nNodePerFace
  read_bem_input_long (fp, &pbin->nNodePerFace, fmt);
  // Get nIFValue and nDFValue
  read_bem_input_long (fp, &pbin->nIFValue, fmt);
  read_bem_input_long (fp, &pbin->nDFValue, fmt);
  // Allocate and get idOfFace
  pbin->idOfFace = (long*)malloc (pbin->nFace * pbin->nNodePerFace * sizeof(long));
  read_bem_input_longs (fp, pbin->idOfFace, pbin->nFace * pbin->nNodePerFace,
			     pbin->nNodePerFace, fmt);
  // Allocate and calculate coordOfFace
  pbin->coordOfFace = (double(*)[3]) malloc (pbin->nFace * 3 * sizeof(double));
  for(long i=0; i < pbin->nFace; i++){
    const long ncpf = pbin->nNodePerFace;
    double (* const con)[3] = pbin->coordOfNode;
    // Initialize coordOfFace
    pbin->coordOfFace[i][0] = pbin->coordOfFace[i][1] = pbin->coordOfFace[i][2] = 0.0;
    // Sum up all coordinates of the i-th face
    for (long j = 0; j<ncpf; j++) {
      long pid = pbin->idOfFace[ncpf*i+j];
      pbin->coordOfFace[i][0] += con[pid][0];
      pbin->coordOfFace[i][1] += con[pid][1];
      pbin->coordOfFace[i][2] += con[pid][2];
    }
    // Divide by ncpf to calcuate the average
    pbin->coordOfFace[i][0] /= (double) ncpf;
    pbin->coordOfFace[i][1] /= (double) ncpf;
    pbin->coordOfFace[i][2] /= (double) ncpf;
  }
  // Allocate and get IFValue
  pbin->IFValue = (long*) malloc (sizeof(long)* pbin->nIFValue * pbin->nFace);
  read_bem_input_longs (fp, pbin->IFValue, pbin->nIFValue * pbin->nFace, 1, fmt);
  // Allocate and get DFValue
  pbin->DFValue = (double*) malloc (sizeof(double)* pbin->nDFValue * pbin->nFace);
  read_bem_input_doubles (fp, pbin->DFValue, pbin->nDFValue * pbin->nFace, 1, fmt);
  return fmt;
}

// Load test
#ifdef TEST
int main (void)
{
  struct bem_input bem_in;
  // Read text file
  FILE* fpin = fopen ("input.txt", "r");
  read_bem_input (fpin, &bem_in, BI_AUTO);
  fclose (fpin);
  // Pretty print
  print_bem_input (stdout, &bem_in, BI_PRETTY);
  // Write as binary file
  FILE* fpout = fopen ("input.txt.bin", "w");
  print_bem_input (fpout, &bem_in, BI_BINARY);
  fclose (fpout);
  // Read the binary file
  FILE* fpin_b = fopen ("input.txt.bin", "r");
  read_bem_input (fpin_b, &bem_in, BI_AUTO);
  fclose (fpin_b);
  // Pretty print
  print_bem_input (stdout, &bem_in, BI_PRETTY);
}
#endif
