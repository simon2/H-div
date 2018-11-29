enum bi_format {
  BI_TEXT,      // oroginal text format
  BI_BINARY,    // binary format
  BI_PRETTY,    // pretty print (only for print)
  BI_AUTO,      // automatic select (only for read)
};

struct bem_input {
  long nNode;               // # of nodes
  double (*coordOfNode)[3]; // coordOfnode[i][d] is d-th dimension of i-th coordinate
                            // where 0 < i < nNode and 0 < d < 3
  long nFace;               // # of faces
  long nNodePerFace;        // # of nodes that form each face     
  long nIFValue;            // # of integer parameter values for each face
  long nDFValue;            // # of double parameter values for each face
  long *idOfFace;           // idOfFace[i*nNodePerFace+j] is j-th node ID of i-th face
                            // where 0 < i < nFace and 0 < j < jNodePerFace
  double (*coordOfFace)[3]; // coordOfFace[i][d] is d-th dimension of i-th face's center
                            // where 0 < i < nFace and 0 < d < 3
  long *IFValue;            // IFValue[i*nFace+j] is i-th integer parameter value
                            // of j-th face where 0<i<nFace and 0<j<nIFValue
  double *DFValue;          // DFValue[i*nFace+j] is i-th double parameter value
                            // of j-th face where 0<i<nFace and 0<j<nIDValue
};

void print_bem_input (FILE* fp, struct bem_input* pbin, enum bi_format fmt);
enum bi_format read_bem_input (FILE* fp, struct bem_input* pbin, enum bi_format fmt);
