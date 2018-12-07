
// Auxillary functions
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include "bem_file.h"

double show_elapsed_time (struct timeval *tv1, struct timeval *tv2) {
  double tm;
  tm =  (tv2->tv_usec - tv1->tv_usec) / 1000000.0;   // [us]
  tm += (tv2->tv_sec - tv1->tv_sec);                 // [s]
  fprintf (stderr, "Time: %lf sec.\n", tm);
  return tm;
}


// If flag is non-zero, execute the sequence of:
// * fopen output file
// * write to the output file in the specified format (fmt)
// * Show the elapsed time and the result text 
void open_write_show (const int flag, struct bem_input *pbin,
		      const enum bi_format fmt,
		      const char* ofile)
{
  if (flag) {
    char *fmt_str;
    switch (fmt) {
    case BI_TEXT:  
      fmt_str = "TEXT";
      break;
    case BI_BINARY:
      fmt_str = "BINARY";
      break;
    case BI_VTK:
      fmt_str = "VTK";
      break;
    default:
      fprintf (stderr, "Unexpexcted fmt in open_write_show.\n");
      exit (99);
    }
    FILE *fpout;
    struct timeval tv1, tv2;
    if ( !(fpout = fopen (ofile, "w")) ) {
      fprintf (stderr, "Output file %s open error!\n", ofile);
      exit (99);
    }
    gettimeofday (&tv1, NULL);
    print_bem_input (fpout, pbin, fmt);
    gettimeofday (&tv2, NULL);
    fprintf (stderr, "Wrote the BEM input in the %s format to \"%s\".\n",
	     fmt_str, ofile);
    show_elapsed_time (&tv1, &tv2);
    fclose (fpout);
  }
}
