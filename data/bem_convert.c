// bem_convert:
// Convert a PPOpen BEM input file (.bin) to the text format
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include "bem_file.h"

#define USAGE_STRING \
"Usage: %s [-o <ofile_hd>] [-tbvs] <ifile>\n" \
"Convert the BEM input file <ifile> in the text/binary format to the binary, text, or vtk format.\n" \
"(Note: the vtk format file can be uset do visualize the data using ParaView.)\n" \
"The output file name is <ofile_hd>+\".bin\", <ofile_hd>+\".txt\", or <ofile_hd>+\".vtk\"\n" \
"If <ofile_hd> is not given, the output file name is <ifile>+\".bin\", <ifile>+\".txt\", or <ifile>+\".vtk\"\n" \
"-t: output the text file.\n" \
"-b: output the binary file.\n" \
"-v: output the vtk file.\n" \
"-s: only show of the summary of the input file.\n" \
"If none of -tbvs is specified, convert the text file to the binary file (and vice versa).\n"

// Options set by command line arguments.
struct options {
  char *ifile;       // input file name
  char *ofile_txt;   // output file name
  char *ofile_bin;   // output file name
  char *ofile_vtk;   // output file name
  int txt, bin, vtk; // Set to 1 if -t, -b, -v is specified
  int show_only;     // Set to 1 if -s is specified
};

void usage (int argc, char** argv) {
  fprintf (stderr, USAGE_STRING, argv[0]);
  return; 
}

void set_options (int argc, char** argv, struct options *popt)
{
  char* ofile_hd = NULL;
  // set default 
  popt->ofile_txt = popt->ofile_bin = popt->ofile_vtk = NULL;
  popt->txt = popt->bin = popt->vtk = 0;
  popt->show_only = 0;
  // Analyze command line arguments.
  int ch;
  while ((ch = getopt(argc, argv, "o:tbvs")) != -1) {
    switch (ch) {
    case 't':
      popt->txt = 1;
      break;
    case 'b':
      popt->bin = 1;
      break;
    case 'v':
      popt->vtk = 1;
      break;
    case 's':
      popt->show_only = 1;
      break;
    case 'o':
      {
	size_t len = strlen (optarg);
	ofile_hd = malloc (sizeof(char)*(len+1));
	strncpy (ofile_hd, optarg, len+1);
      }
      break;
    default:
      fprintf (stderr, "Illegal option: %c\n", ch);
      usage (argc, argv);
      exit (99);
      break;
    }
  }
  // Set input file name.
  if (optind >= argc) {
    fprintf (stderr, "No input file name specified.\n");
    usage (argc, argv);
    exit (99);
  }
  popt->ifile = argv[optind++];
  // Set output file names.
  if (!ofile_hd) {
    size_t len = strlen (popt->ifile);
    ofile_hd = malloc (sizeof(char)*(len+1));
    strncpy (ofile_hd, popt->ifile, len+1);
  }
  {
    size_t len = strlen (ofile_hd);
    size_t sz = len+5; // ifile + ".txt" + \n
    popt->ofile_txt = malloc (sizeof(char)*sz);
    popt->ofile_bin = malloc (sizeof(char)*sz);
    popt->ofile_vtk = malloc (sizeof(char)*sz);
    snprintf (popt->ofile_txt, sz, "%s" ".txt", ofile_hd);
    snprintf (popt->ofile_bin, sz, "%s" ".bin", ofile_hd);
    snprintf (popt->ofile_vtk, sz, "%s" ".vtk", ofile_hd);
    free (ofile_hd);
  }    
  return;
}

int main (int argc, char **argv)
{
  struct options opt;
  struct bem_input bin;           // The object to which the input data are stored.
  enum bi_format fmt_i;           // The formats of input/output file.
  
  set_options (argc, argv, &opt);

  // Read
  {
    struct timeval tv1, tv2;
    gettimeofday (&tv1, NULL);
    fmt_i = open_and_read_bem_input (opt.ifile, &bin, BI_AUTO);
    gettimeofday (&tv2, NULL);
    // Check the format of the input file and set the format of the output file.
    switch (fmt_i) {
    case BI_TEXT:
      fprintf (stderr, "Read the BEM input in the TEXT format from \"%s\".\n", opt.ifile);
      if (!(opt.txt || opt.bin || opt.vtk) && !opt.show_only ) {
	opt.bin = 1;
      }
      break;
    case BI_BINARY:
      fprintf (stderr, "Read the BEM input in the BINARY format from \"%s\".\n", opt.ifile);
      if (!(opt.txt || opt.bin || opt.vtk) && !opt.show_only ) {
	opt.txt = 1;
      }
      break;
    default:
      fprintf (stderr, "Read error occurred when reading %s.!\n", opt.ifile);
      exit (99);
    }
    show_elapsed_time (&tv1, &tv2);
  }
  
  // Write output files if needed.
  open_write_show (opt.txt, &bin, BI_TEXT,   opt.ofile_txt);
  open_write_show (opt.bin, &bin, BI_BINARY, opt.ofile_bin);
  open_write_show (opt.vtk, &bin, BI_VTK,    opt.ofile_vtk);
  // Print summary of the file to STDOUT
  print_bem_input (stdout, &bin, BI_PRETTY);

  return 0;
}
