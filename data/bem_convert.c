// Convert a PPOpen BEM input file (.bin) to the text format
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "bem_file.h"


/*
Usage:
% ./bin2txt <filename>
Convert the BEM input file <filename> in the text/binary format to the binary/text format.
If <filename> is given, the result is output to the file named <filename>+".bin" or <filename>+".txt"
When <filename> is not specified, stdin and stdout are used as the input and output.
 */

double show_elapsed_time (struct timeval *tv1, struct timeval *tv2) {
  double tm;
  tm =  (tv2->tv_usec - tv1->tv_usec) / 1000000.0;   // [us]
  tm += (tv2->tv_sec - tv1->tv_sec);                 // [s]
  fprintf (stderr, "Time: %lf sec.\n", tm);
  return tm;
}

int main (int argc, char **argv)
{
  char *infile = "(standard input)";     // input file name
  char *outfile = "(standard output)";   // output file name
  FILE *fpin = stdin;     // input stream
  FILE *fpout = stdout;   // output stream
  struct timeval t_in1, t_in2, t_out1, t_out2;
  
  if (argc > 1) {
    infile = argv[1];
    if ( !(fpin = fopen (infile, "r")) ) {
      fprintf (stderr, "Input file open error!\n");
      exit (99);
    } 
  }
  
  struct bem_input bin;           // The object to which the input data are stored.
  enum bi_format fmt_i, fmt_o;    // The formats of input/output file.

  // Read
  gettimeofday (&t_in1, NULL);
  fmt_i = read_bem_input (fpin, &bin, BI_AUTO);
  gettimeofday (&t_in2, NULL);
  fclose (fpin);

  // Check the format of the input file and set the format of the output file.
  switch (fmt_i) {
  case BI_TEXT:
    fprintf (stderr, "Read the BEM input in the TEXT format from \"%s\".\n", infile);
    fmt_o = BI_BINARY;
    break;
  case BI_BINARY:
    fprintf (stderr, "Read the BEM input in the BINARY format from \"%s\".\n", infile);
    fmt_o = BI_TEXT;
    break;
  default:
    fprintf (stderr, "Read error occurred when reading %s.!\n", infile);
    exit (99);
  }
  show_elapsed_time (&t_in1, &t_in2);

  // Set output stream
  if (argc > 1) { // input file name is given
    outfile = (char*) malloc (sizeof(char)*(strlen(infile)+10));
    switch (fmt_o) {
    case BI_TEXT:
      // outfile = infile + ".txt"
      strcpy (outfile, infile); strcat (outfile, ".txt");
      break;
    case BI_BINARY:
      // outfile = infile + ".bin"
      strcpy (outfile, infile); strcat (outfile, ".bin");
      break;
    default:
      fprintf (stderr, "Unexpected error!\n");
      exit (99);
    }
    if ( !(fpout = fopen (outfile, "w")) ) {
      fprintf (stderr, "Output file open error!\n");
      exit (99);
    }
  }

  // Write
  gettimeofday (&t_out1, NULL);
  print_bem_input (fpout, &bin, fmt_o);
  gettimeofday (&t_out2, NULL);
  fclose (fpout);

  // Show result.
  switch (fmt_o) {
  case BI_TEXT:
    fprintf (stderr, "Wrote the BEM input in the TEXT format to \"%s\".\n", outfile);
    break;
  case BI_BINARY:
    fprintf (stderr, "Wrote the BEM input in the BINARY format from \"%s\".\n", outfile);
    break;
  default:
    fprintf (stderr, "Unexpected error!\n");
    exit (99);
  }
  show_elapsed_time (&t_out1, &t_out2);

  return 0;
}
