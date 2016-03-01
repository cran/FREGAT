// Functions from package 'snpStats' v.1.20.0 (c) 2015 David Clayton

#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <R.h>

/* One byte coding of SNP genotype allowing for uncertainty */

const unsigned char lup0[253] = 
  {1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,18, 19, 20, 21, 22, 23, 
   3, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 
   42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 
   61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
   80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 
   99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 
   114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
   129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 
   144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 
   159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 
   174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 
   189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 
   204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 
   219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 
   234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 
   249, 250, 251, 252, 253, 2};

const int lup1[253] = 
  {0, 252, 21, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 
   19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 
   39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 
   58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 
   77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
   96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 
   112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 
   127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 
   142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 
   157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 
   172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 
   187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 
   202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 
   217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 
   232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 
   247, 248, 249, 250, 251};

const double lup2[253] = 
  {0./21., 0./21., 0./21., 0./21., 0./21., 0./21., 0./21., 0./21., 0./21., 
   0./21., 0./21., 0./21., 0./21., 0./21., 0./21., 0./21., 0./21., 0./21., 
   0./21., 0./21., 0./21., 0./21., 1./21., 1./21., 1./21., 1./21., 1./21., 
   1./21., 1./21., 1./21., 1./21., 1./21., 1./21., 1./21., 1./21., 1./21., 
   1./21., 1./21., 1./21., 1./21., 1./21., 1./21., 1./21., 2./21., 2./21., 
   2./21., 2./21., 2./21., 2./21., 2./21., 2./21., 2./21., 2./21., 2./21., 
   2./21., 2./21., 2./21., 2./21., 2./21., 2./21., 2./21., 2./21., 2./21., 
   3./21., 3./21., 3./21., 3./21., 3./21., 3./21., 3./21., 3./21., 3./21., 
   3./21., 3./21., 3./21., 3./21., 3./21., 3./21., 3./21., 3./21., 3./21., 
   3./21., 4./21., 4./21., 4./21., 4./21., 4./21., 4./21., 4./21., 4./21., 
   4./21., 4./21., 4./21., 4./21., 4./21., 4./21., 4./21., 4./21., 4./21., 
   4./21., 5./21., 5./21., 5./21., 5./21., 5./21., 5./21., 5./21., 5./21., 
   5./21., 5./21., 5./21., 5./21., 5./21., 5./21., 5./21., 5./21., 5./21., 
   6./21., 6./21., 6./21., 6./21., 6./21., 6./21., 6./21., 6./21., 6./21., 
   6./21., 6./21., 6./21., 6./21., 6./21., 6./21., 6./21., 7./21., 7./21., 
   7./21., 7./21., 7./21., 7./21., 7./21., 7./21., 7./21., 7./21., 7./21., 
   7./21., 7./21., 7./21., 7./21., 8./21., 8./21., 8./21., 8./21., 8./21., 
   8./21., 8./21., 8./21., 8./21., 8./21., 8./21., 8./21., 8./21., 8./21., 
   9./21., 9./21., 9./21., 9./21., 9./21., 9./21., 9./21., 9./21., 9./21., 
   9./21., 9./21., 9./21., 9./21., 10./21., 10./21., 10./21., 10./21., 10./21.,
   10./21., 10./21., 10./21., 10./21., 10./21., 10./21., 10./21., 11./21., 
   11./21., 11./21., 11./21., 11./21., 11./21., 11./21., 11./21., 11./21., 
   11./21., 11./21., 12./21., 12./21., 12./21., 12./21., 12./21., 12./21., 
   12./21., 12./21., 12./21., 12./21., 13./21., 13./21., 13./21., 13./21., 
   13./21., 13./21., 13./21., 13./21., 13./21., 14./21., 14./21., 14./21., 
   14./21., 14./21., 14./21., 14./21., 14./21., 15./21., 15./21., 15./21., 
   15./21., 15./21., 15./21., 15./21., 16./21., 16./21., 16./21., 16./21., 
   16./21., 16./21., 17./21., 17./21., 17./21., 17./21., 17./21., 18./21., 
   18./21., 18./21., 18./21., 19./21., 19./21., 19./21., 20./21., 20./21., 
   21./21};

const double lup3[253] = 
  {0./21., 1./21., 2./21., 3./21., 4./21., 5./21., 6./21., 7./21., 8./21., 
   9./21., 10./21., 11./21., 12./21., 13./21., 14./21., 15./21., 16./21., 
   17./21., 18./21., 19./21., 20./21., 21./21., 0./21., 1./21., 2./21., 
   3./21., 4./21., 5./21., 6./21., 7./21., 8./21., 9./21., 10./21., 11./21., 
   12./21., 13./21., 14./21., 15./21., 16./21., 17./21., 18./21., 19./21., 
   20./21., 0./21., 1./21., 2./21., 3./21., 4./21., 5./21., 6./21., 7./21., 
   8./21., 9./21., 10./21., 11./21., 12./21., 13./21., 14./21., 15./21., 
   16./21., 17./21., 18./21., 19./21., 0./21., 1./21., 2./21., 3./21., 4./21., 
   5./21., 6./21., 7./21., 8./21., 9./21., 10./21., 11./21., 12./21., 13./21., 
   14./21., 15./21., 16./21., 17./21., 18./21., 0./21., 1./21., 2./21., 3./21., 
   4./21., 5./21., 6./21., 7./21., 8./21., 9./21., 10./21., 11./21., 12./21., 
   13./21., 14./21., 15./21., 16./21., 17./21., 0./21., 1./21., 2./21., 3./21., 
   4./21., 5./21., 6./21., 7./21., 8./21., 9./21., 10./21., 11./21., 12./21., 
   13./21., 14./21., 15./21., 16./21., 0./21., 1./21., 2./21., 3./21., 4./21., 
   5./21., 6./21., 7./21., 8./21., 9./21., 10./21., 11./21., 12./21., 13./21., 
   14./21., 15./21., 0./21., 1./21., 2./21., 3./21., 4./21., 5./21., 6./21., 
   7./21., 8./21., 9./21., 10./21., 11./21., 12./21., 13./21., 14./21., 0./21.,
   1./21., 2./21., 3./21., 4./21., 5./21., 6./21., 7./21., 8./21., 9./21., 
   10./21., 11./21., 12./21., 13./21., 0./21., 1./21., 2./21., 3./21., 4./21., 
   5./21., 6./21., 7./21., 8./21., 9./21., 10./21., 11./21., 12./21., 0./21., 
   1./21., 2./21., 3./21., 4./21., 5./21., 6./21., 7./21., 8./21., 9./21., 
   10./21., 11./21., 0./21., 1./21., 2./21., 3./21., 4./21., 5./21., 6./21., 
   7./21., 8./21., 9./21., 10./21., 0./21., 1./21., 2./21., 3./21., 4./21., 
   5./21., 6./21., 7./21., 8./21., 9./21., 0./21., 1./21., 2./21., 3./21., 
   4./21., 5./21., 6./21., 7./21., 8./21., 0./21., 1./21., 2./21., 3./21., 
   4./21., 5./21., 6./21., 7./21., 0./21., 1./21., 2./21., 3./21., 4./21., 
   5./21., 6./21., 0./21., 1./21., 2./21., 3./21., 4./21., 5./21., 0./21., 
   1./21., 2./21., 3./21., 4./21., 0./21., 1./21., 2./21., 3./21., 0./21., 
   1./21., 2./21., 0./21., 1./21., 0./21.};

/* Mapping of posterior probabilities into one byte code */

unsigned char post2g(double pAB, double pBB) {
  double pAA = 1.0 - pAB - pBB;
  /* Round to integers in (0,21) */
  pAA *= 21.;
  pAB *= 21.;
  pBB *= 21.;
  int iAA = floor(pAA+0.499999);
  int iAB = floor(pAB+0.499999);
  int iBB = floor(pBB+0.499999);
  int tot = iAA+iAB+iBB;
  if (tot!=21) { /* Then simple rounding doesn't do it */
    /* Remainders */
    pAA -= iAA*21.;
    pAB -= iAB*21.;
    pBB -= iBB*21.;
    if (tot<21) {
      if (pAB>pBB && pAB>pAA)
	iAB++;
      else if (pBB>pAB && pBB>pAA)
	iBB++;
      else
	iAA++;
    }
    else {
      if (pAB<pBB && pAB<pAA)
	iAB--;
      else if (pBB<pAB && pBB<pAA)
	iBB--;
      else
	iAA--;
    }
  }
  if ((iAA + iAB + iBB)!=21)
    error("Bug -- illegal sum");
  /* Simple sequential coding 0:252 */
  int first = 253 - ((22-iAB)*(23-iAB))/2 + iBB;
  /* Recode to 1:253, so that vertices are at 1, 2, 3 */
  return lup0[first]; 
}

/* Mapping of one byte code into posterior uncertainty */

int g2post(unsigned char code, double *pAA, double *pAB, double *pBB) {
  if (!code || code>253)
    return(0);
  int first = lup1[code-1];
  *pAB = lup2[first];
  *pBB = lup3[first];
  *pAA = 1.0-(*pAB)-(*pBB);
  return(1);
}

/* Mapping of posterior mean genotype into one byte code */

unsigned char mean2g(double mean, int maxE) {
  if (mean<0.0 || mean>2.0)
    return((unsigned char) 0);
  if (mean==0.0 || mean==2.0)
    return((unsigned char) (1+mean));
  if (maxE) { /* Maximum entropy assignment */
    /* Probabilities are theta, theta*phi, and theta*phi^2
       Phi is the positive root of quadratic */
    double m1 = mean - 1.0;
    double phi = (sqrt(4.0 - 3.0*m1*m1) + m1)/(2*(1 - m1));
    double phi2 = phi*phi;
    double theta = 1.0/(1.0 + phi + phi2);
    return(post2g(theta*phi, theta*phi2));
  }
  else { /* Least uncertain */
    if (mean>=1.0)
      return(post2g(2.0-mean, mean-1.0));
    else 
      return(post2g(mean, 0.0));
  }
}

/* One byte code into posterior mean */

double g2mean(unsigned char code){
  if (!code || code>253)
    return(-1.0);
  if (code<4)
    return( (double) (code-1));
  int first = lup1[code-1];
  return(lup2[first]+2*lup3[first]);
}

/* One byte code into additive and dominance dummy variables */

int g2ad(unsigned char code, double *xadd, double *xdom) {
  if (!code || code>253) 
    return(1); /* Fault indicator */
  if (code<4) { /* Known genotype */
    *xadd = (double) (code-1);
    *xdom = (double) (code==3);
  }
  else { /* Uncertain genotype */
    int first = lup1[code-1];
    *xadd = lup2[first]+2*lup3[first];
    *xdom = lup3[first];
  }
  return(0);
}

