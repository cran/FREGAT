// Functions from package 'snpStats' v.1.20.0 (c) 2015 David Clayton

/* Modified for R_xlen_t 26/06/2015 */

/* 
   Read a plink .bed file as a SnpMatrix
 
*/


#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* #include "Rmissing.h" */


/* this may not be fast but seems reliable */

void skip(FILE *in, int lines, int size) {
  if (lines) {
    for (int i=0; i<lines; i++) {
      for (int j=0; j<size; j++) {
	fgetc(in);
	if (feof(in)) error("unexpected end of file");
      }
    }
  }    
}


SEXP readbed(SEXP Bed, SEXP Id, SEXP Snps, SEXP Rsel, SEXP Csel) {
  const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
  const unsigned char mask = '\x03';
  int nrow = LENGTH(Id);
  int ncol = LENGTH(Snps);
  const char *file = CHAR(STRING_ELT(Bed, 0));
  FILE *in = fopen(file, "rb");
  if (!in)
    error("Couldn't open input file: %s", file);
  unsigned char start[3];
  if (fread(start, 1, 3, in)!=3)
    error("Failed to read first 3 bytes");
  if (start[0]!='\x6C' || start[1]!='\x1B')
    error("Input file does not appear to be a .bed file (%X, %X)", 
	  start[0], start[1]);

  /* Create output object */

  SEXP Result, Dimnames, Package, Class;
  PROTECT(Result = allocMatrix(RAWSXP, nrow, ncol));
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Dimnames, 0, Id);
  SET_VECTOR_ELT(Dimnames, 1, Snps);
  setAttrib(Result, R_DimNamesSymbol, Dimnames);
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("SnpMatrix"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("FREGAT"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  SET_S4_OBJECT(Result);
  
  unsigned char *result = RAW(Result); 
  R_xlen_t ncell = (R_xlen_t)nrow*(R_xlen_t)ncol;
  memset(result, 0x00, ncell);

  /* Read in first few bytes */

  int *seek = NULL;
  int snp_major = start[2];
  int nbyte = 0;
  if (snp_major) {
    if (!isNull(Rsel))
      error("can't select by rows when .bed file is in SNP-major order");
    if (!isNull(Csel)) {
      seek = INTEGER(Csel);
      nbyte = 1 + (nrow-1)/4; 
    }
  }
  else {
    if (!isNull(Csel))
      error("can't select by columns when .bed file is in individual-major order");
    if (!isNull(Rsel)) {
      seek = INTEGER(Rsel);
      nbyte = 1 + (ncol-1)/4;
    }
  }

  if (seek) 
    skip(in, seek[0]-1, nbyte);
  int part=0, i=0, j=0;
  R_xlen_t ij = 0;
  unsigned char byte = 0x00;
  while (1) {
    if (!part) {
      byte = (unsigned char) fgetc(in);
      if (feof(in))
	error("Unexpected end of file reached");
      part = 4;
    }
    unsigned char code = byte & mask;
    byte = byte >> 2;
    part--;
    result[ij] = recode[code];
    if (snp_major) {
      ij++;
      i++;
      if (i==nrow) {
	i = part = 0;
	j++;
	if (j==ncol)
	  break;
	if (seek)
	  skip(in, seek[j]-seek[j-1]-1, nbyte);
      }
    }	
    else {
      ij += nrow;
      j++;
      if (j==ncol){
	j = part = 0;
	i++;
	if (i==nrow)
	  break;
	if (seek)
	  skip(in, seek[i]-seek[i-1]-1, nbyte);
	ij = i;
      }
    }
  }
  fclose(in);
  UNPROTECT(4);
  return Result;
}

