/* History: Apr 20 2021 Initial coding
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define MYDEBUG 0
#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}
#define TWOPI 6.2831853071795862
#define LOG2 0.69314718055994529
#define LOG0ARG 1.0e-300
#define LOG0 -690.77552789821368
#define MISSING 9
#define NUMERICZERO 1.0e-300
#define LARGENUMERIC 1.0e300
#define LARGEINT 999999
#define LARGENEGINT -999999
#define MYMIN(a, b) (((a) < (b)) ? (a) : (b))
#define MYMAX(a, b) (((a) > (b)) ? (a) : (b))

#define IARG_DEBUG 0
#define IARG_data_nr 1
#define IARG_data_nc 2
#define IARG_freqVecLen 3
#define IARG_maxiter 4
#define IARG_print 5
#define IARG_kLen 6

#define DARG_r_0 0
#define DARG_p_00 1
#define DARG_p_11 2
#define DARG_pi_0 3
#define DARG_pi_1 4
#define DARG_epsilon 5
#define DARG_machine_eps 6

struct mystr {
  int DEBUG;
  int *data;            /* matrix passed in by rows */
  int data_nr;          /* # of rows */
  int data_nc;          /* # of cols */
  int **data_rowPtr;    /* Pointers to the rows of data */
  int *freqVec;         /* frequency vector passed in */
  int freqVecLen;       /* length of freqVec */
  int freqVecSum;       /* sum(freqVec) */
  double r_0;
  double p_00;
  double p_11;
  double pi_0;
  double pi_1;
  double epsilon;       /* Stopping tolerance for EM alg */
  double machine_eps;   /* Machine epsilon for determining if numbers are equal */
  int maxiter;          /* Max iterations for EM alg */
  int print;
  int converged;
  
  int *k;               /* vector of states */
  int kLen;             /* length of vector k */
  double *init_prob;    /* initial probabilities, set to pointer passed in */
  double **trans;       /* transition matrix */
  double **forwardMat;
  double **backwardMat;
  int **classMat;       /* classification as a lookup table */  

  
  char *data_rowMaxEqual0;  /* 0-1 vector if the max value of a row is 0 */
  char *data_rowMinEqual1;  /* 0-1 vector if the min value of a row is 1 */

  /* scratch objects */
  double **tmpMat;
  double *tmpVec;
  double *tmpVec2;
  double **trans_tmp;   /* temp matrix for computing trans */

  /* For local decode call */
  int *localDecodeVec;  /* Points to input (return) object (ret_vec) */
  int random_draw;      /* 0 or 1 */

};
typedef struct mystr MYSTR;


/* Used for debugging */
/*
static void writeVec(fid, vec, len, ncol)
FILE *fid;
double *vec;
int len, ncol;
{
  int i, m=0;

  for (i=0; i<len-1; i++) {
    fprintf(fid, "%20.11g\t", vec[i]);
    m++;
    if (ncol && (m == ncol)) {
      m = 0;
      fprintf(fid, "\n");
    }
  }
  fprintf(fid, "%20.11g\n", vec[len-1]);

} 
*/

void print_dVec(vec, n, n2, name)
double *vec;
int n, n2;
char name[10];
{
  int i, m=0;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %g ", vec[i]);
    m++;
    if (n2 && (m == n2)) {
      Rprintf("\n");
      m = 0;
    }
  }
  Rprintf("\n \n");
}

void print_iVec(vec, n, n2, name)
int *vec;
int n, n2;
char name[10];
{
  int i, m=0;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %d ", vec[i]);
    m++;
    if (n2 && (m == n2)) {
      Rprintf("\n");
      m = 0;
    }
  }
  Rprintf("\n \n");
}

void print_dMat(mat, nr, nc, name)
double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) Rprintf(" %g ", mat[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}
void print_sumAbsVec(vec, n) 
double *vec;
int n;
{
  int i;
  double sum=0.0;
  for (i=0; i<n; i++) sum += fabs(vec[i]);
  Rprintf("absSum=%g %40.8f\n", sum, sum);

}
void print_sumAbsMat(mat, nr, nc) 
double **mat;
int nr, nc;
{
  int i, j;
  double sum=0.0;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) sum += fabs(mat[i][j]);
  }
  Rprintf("absSum=%g %40.8f\n", sum, sum);

}


static void setVecToVal(vec, n, val)
double *vec, val;
int n;
{
  int i;
  for (i=0; i<n; i++) vec[i] = val;

}
static void setMatToVal(mat, nr, nc, val)
double **mat, val;
int nr, nc;
{
  int i, j;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) mat[i][j] = val;
  }
}

static char * cVec_alloc(n, initFlag, initVal)
int n, initFlag;
char initVal;
{
  int i;
  char *ret, *p;

  ret = (char *) malloc(n*sizeof(char));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} 

/* Function to allocate memory for a double vector */
static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */

/* Function to allocate memory for an int vector */
static int * iVec_alloc(n, initFlag, initVal)
int n, initFlag, initVal;
{
  int i, *ret, *p;

  ret = (int *) malloc(n*sizeof(int));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: iVec_alloc */

/*
static double ** ptrdVec_alloc(n, initFlag)
int n, initFlag;
{
  int i;
  double **ret, **p;

  ret = (double **) malloc(n*sizeof(double *));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = NULL;
  }

  return(ret);

} 
*/

static int ** ptriVec_alloc(n, initFlag)
int n, initFlag;
{
  int i, **ret, **p;

  ret = (int **) malloc(n*sizeof(int *));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = NULL;
  }

  return(ret);

} 

/* Function to allocate a double matrix */
static double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: dMat_alloc */

/* Function to allocate a double matrix */
static int ** iMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag, initVal;
{
  int **mat, **ptr, i;

  mat = (int **) malloc(nrow*sizeof(int *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = iVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: iMat_alloc */


/* Function to allocate a double array (3-dim) */
/*
static double *** dArray_alloc(n1, n2, n3, initFlag, initVal)
int n1, n2, n3, initFlag;
double initVal;
{
  double ***mat;
  int i, j;

  mat = (double ***) malloc(n1*sizeof(double **));
  CHECK_MEM(mat);
  for (i=0; i<n1; i++) {
    mat[i] = (double **) malloc(n2*sizeof(double *));
    CHECK_MEM(mat[i]);
    for (j=0; j<n2; j++) mat[i][j] = dVec_alloc(n3, initFlag, initVal);
  }
     
  return(mat);

}  */

/* Function to free a matrix */
static void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

/* Function to free an aray (3d) */
/*
static void array_free(x, n1, n2)
void ***x;
int n1, n2;
{
  int i, j;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) free(x[i][j]);
    free(x[i]);
  }
  free(x);

}  */


/* Function to fill in a matrix from a vector (by column) */
static void fillMat(vec, nr, nc, out)
double *vec, **out;
int nr, nc;
{
  int i, j, col=0, ii;

  ii = 0;
  for (j=0; j<nc; j++) {
    for (i=0; i<nr; i++) {
      out[i][col] = vec[ii];
      ii++;
    }
    col++;
  }

} 

/*
static double maxAbsDiff(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double val, maxval=-1.0, *p1, *p2;

  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) {
    val = fabs(*p1 - *p2);
    if (val > maxval) maxval = val;
  }

  return(maxval);

} 
*/

static void copy_dVec(from, into, n)
double *from, *into;
int n;
{
  int i;
  double *p1, *p2;

  for (i=0, p1=from, p2=into; i<n; i++, p1++, p2++) *p2 = *p1;

} /* END: copy_dVec */

static double sumvec(vec, n)
double *vec;
int n;
{
  int i;
  double ret=0.0, *p;

  for (i=0, p=vec; i<n; i++, p++) ret += *p;

  return(ret);

}

static int sumivec(vec, n)
int *vec, n;
{
  int i, ret=0, *p;

  for (i=0, p=vec; i<n; i++, p++) ret += *p;

  return(ret);

}

static double dotprod(vec1, vec2, n)
double *vec1, *vec2;
int n;
{
  int i;
  double ret=0.0, *p1, *p2;

  for (i=0, p1=vec1, p2=vec2; i<n; i++, p1++, p2++) ret += *p1 * *p2;

  return(ret);

}

static void vecprod(vec1, vec2, n, ret)
double *vec1, *vec2, *ret;
int n;
{
  int i;
  double *p1, *p2, *p3;

  for (i=0, p1=vec1, p2=vec2, p3=ret; i<n; i++, p1++, p2++, p3++) *p3 = *p1 * *p2;

}

static void vecprod_di(vec1, vec2, n, ret)
double *vec1, *ret;
int *vec2, n;
{
  int i, *p2;
  double *p1, *p3;

  for (i=0, p1=vec1, p2=vec2, p3=ret; i<n; i++, p1++, p2++, p3++) *p3 = *p1 * *p2;

}

static void putColInVec(mat, col, n, ret)
double **mat, *ret;
int col, n;
{
  int i;

  for (i=0; i<n; i++) ret[i] = mat[i][col];

}

/*
static void logvec(vec, n, ret)
double *vec, *ret;
int n;
{
  int i;
  double *p1, *p2;

  for (i=0, p1=vec, p2=ret; i<n; i++, p1++, p2++) *p2 = log(*p1);

}
static void reverseCumSum(vec, n, ret)
double *vec, *ret; 
int n;
{
  int i;
  double sum;

  if (ret == vec) error("INTERNAL CODING ERROR in reverseCumSum");
  sum    = sumvec(vec, n);

  ret[0] = sum;
  for (i=1; i<n; i++) ret[i] = ret[i-1] - vec[i-1]; 

}
*/

static int maxiVecMiss(x, n)
int *x, n;
{
  int i, ret=LARGENEGINT, xi;

  for (i=0; i<n; i++) {
    xi = x[i];
    if ((xi != MISSING) && (xi > ret)) ret = xi;  
  }

  return(ret);

}

static int miniVecMiss(x, n)
int *x, n;
{
  int i, ret=LARGEINT, xi;

  for (i=0; i<n; i++) {
    xi = x[i];
    if ((xi != MISSING) && (xi < ret)) ret = xi;  
  }

  return(ret);

}

static void divideVecBy(x, n, value, ret)
double *x, value, *ret;
int n;
{
  int i;
  double div, *p, *pret;

  div = 1.0/value;
  for (i=0, p=x, pret=ret; i<n; i++, p++, pret++) *pret = *p * div;

}

static void putMatinRowVec(mat, nr, nc, ret)
double **mat, *ret;
int nr, nc;
{
  int i, j, k=0;

  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) ret[k++] = mat[i][j];
  }

}

/**************************************************************************
***************************************************************************
***************************************************************************
***************************************************************************
***************************************************************************/

static int areEqual(x, y, eps)
double x, y, eps;
{

  if (fabs(x - y) <= MYMIN(fabs(x), fabs(y))*eps) {
    return(1);
  } else {
    return(0);
  }

}

static void forwardIter(mystr, time, individual, ret)
MYSTR *mystr;
int time, individual;
double **ret;
{
  if (mystr->DEBUG) Rprintf("Begin: forwardIter\n");

  int i, j, kLen=mystr->kLen, *k=mystr->k, *dat, dati, **classMat=mystr->classMat;
  double **trans=mystr->trans, *vec, *vec0, *init_prob=mystr->init_prob;
  double *tmpVec=mystr->tmpVec;

  dat  = mystr->data_rowPtr[individual];
  dati = dat[0];
  vec  = ret[0];

  if (dati != MISSING) {
    for (i=0; i<kLen; i++) vec[i] = init_prob[i]*classMat[dati][k[i]];
  } else {
    copy_dVec(init_prob, vec, kLen);
  }

  if (time > 0) {
    for (i=1; i<=time; i++) {
      dati = dat[i];
      vec  = ret[i];
      vec0 = ret[i - 1];
      if (dati != MISSING) {
        for (j=0; j<kLen; j++) {
          putColInVec(trans, j, kLen, tmpVec);
          vec[j] = dotprod(vec0, tmpVec, kLen)*classMat[dati][k[j]];
        }
      } else {
        for (j=0; j<kLen; j++) {
          putColInVec(trans, j, kLen, tmpVec);
          vec[j] = dotprod(vec0, tmpVec, kLen);
        }
      }
    }
  }
  if (mystr->DEBUG) Rprintf("End: forwardIter\n");

}

static void backwardIter(mystr, time, individual, ret)
MYSTR *mystr;
int time, individual;
double **ret;
{
  if (mystr->DEBUG) Rprintf("Begin: backwardIter\n");

  int i, j, kLen=mystr->kLen, length=mystr->data_nc, ip1;
  int *dat, dati, **classMat=mystr->classMat;
  double **trans=mystr->trans,  *vec, *vec0;
  double eps=mystr->machine_eps, *tmpVec=mystr->tmpVec;
 
  dat = mystr->data_rowPtr[individual];
  vec = ret[length-1];
  setVecToVal(vec, kLen, 1.0);

  if (time < length - 1) {
    for (i=length-2; i>=time; i--) {
      ip1  = i + 1;
      dati = dat[ip1];
      vec  = ret[i];
      vec0 = ret[ip1];
     
      if (dati != MISSING) {
        vecprod_di(vec0, classMat[dati], kLen, tmpVec);
        for (j=0; j<kLen; j++) vec[j] = dotprod(tmpVec, trans[j], kLen);
      } else {
        if (areEqual(sumvec(vec0, kLen), (double) kLen, eps)) {
          setVecToVal(vec, kLen, 1.0);
        } else {
          for (j=0; j<kLen; j++) vec[j] = dotprod(vec0, trans[j], kLen);   
        }
      }
    }
  }

  if (mystr->DEBUG) Rprintf("Begin: backwardIter\n");
}

static double calcIndLikelihood(mystr, individual)
MYSTR *mystr;
int individual;
{
  if (mystr->DEBUG) Rprintf("Begin: calcIndLikelihood\n");

  double fsum, **tmpMat=mystr->tmpMat, ret, cons;
  double pi_0=mystr->pi_0, pi_1=mystr->pi_1;
  int kLen=mystr->kLen, nc=mystr->data_nc;
  char *maxEq0=mystr->data_rowMaxEqual0, *minEq1=mystr->data_rowMinEqual1;

  forwardIter(mystr, nc-1, individual, tmpMat);
  fsum = sumvec(tmpMat[nc-1], kLen);
  cons = 1.0 - pi_0 - pi_1;

  if (maxEq0[individual]) {
    ret = pi_0 + cons*fsum;
  } else if (minEq1[individual]) {
    ret = pi_1 + cons*fsum;
  } else {
    ret = cons*fsum;
  }

  if (mystr->DEBUG) Rprintf("End: calcIndLikelihood\n");
  return(ret);

}

static double calcLikelihoodPattern(mystr)
MYSTR *mystr;
{
  if (mystr->DEBUG) Rprintf("Begin: calcLikelihoodPattern\n");

  double ret=0.0, tmp;
  int i, *freqVec=mystr->freqVec;

  for (i=0; i<mystr->data_nr; i++) {
    tmp  = calcIndLikelihood(mystr, i);
    ret += log(tmp)*freqVec[i];
  }

  if (mystr->DEBUG) Rprintf("End: calcLikelihoodPattern\n");
  return(ret);

}

/* Updates init_prob in mystr */
static void calcInitialPattern(mystr)
MYSTR *mystr;
{
  double *init_prob=mystr->init_prob, cons, like;
  int i, j, kLen=mystr->kLen, *freqVec=mystr->freqVec;
  double *tmpVec=mystr->tmpVec, tmp, fvi, *tmpVec2=mystr->tmpVec2;
  double **forwardMat=mystr->forwardMat, **backwardMat=mystr->backwardMat;

  cons = 1.0 - mystr->pi_0 - mystr->pi_1;
  setVecToVal(tmpVec2, kLen, 0.0);

  for (i=0; i<mystr->data_nr; i++) {
    like  = calcIndLikelihood(mystr, i);  /* calls forwardIter using tmpMat */
    fvi   = freqVec[i];

    forwardIter(mystr, 0, i, forwardMat);
    backwardIter(mystr, 0, i, backwardMat);   

    tmp = cons/like;
    vecprod(forwardMat[0], backwardMat[0], kLen, tmpVec); /* forward[[1]] * backward[[1]] in R code */
    for (j=0; j<kLen; j++) tmpVec2[j] += tmpVec[j]*tmp*fvi;
  }
  tmp = sumvec(tmpVec2, kLen);
  for (j=0; j<kLen; j++) init_prob[j] = tmpVec2[j]/tmp;

}

static void normalizeRows(mat, nr, nc, ret)
double **mat, **ret;
int nr, nc;
{
  int i;
  double sum;

  for (i=0; i<nr; i++) {
    sum = sumvec(mat[i], nc);
    if (fabs(sum) > NUMERICZERO) divideVecBy(mat[i], nc, sum, ret[i]); 
  }

}

/* Updates trans in mystr */
static void calcTransitionPattern(mystr)
MYSTR *mystr;
{
  double cons, **trans=mystr->trans, fvi, num, denom, **trans_tmp=mystr->trans_tmp, frac;
  double **forwardMat=mystr->forwardMat, **backwardMat=mystr->backwardMat;
  int ind, time, nr=mystr->data_nr, nc=mystr->data_nc, kLen=mystr->kLen, **classMat=mystr->classMat;
  int init_state, new_state, tm1, cls, *freqVec=mystr->freqVec, **rowPtr=mystr->data_rowPtr, *datarow, dati;

  cons = 1.0 - mystr->pi_0 - mystr->pi_1;
  setMatToVal(trans_tmp, kLen, kLen, 0.0);

  for (ind=0; ind<nr; ind++) {
    denom = calcIndLikelihood(mystr, ind);  /* updates tmpMat */
    if (denom < NUMERICZERO) continue;

    datarow = rowPtr[ind];
    fvi     = freqVec[ind];
    for (time=1; time<nc; time++) {
      dati = datarow[time];
      tm1  = time - 1;
      forwardIter(mystr, tm1, ind, forwardMat);
      backwardIter(mystr, time, ind, backwardMat);
      for (init_state=0; init_state<kLen; init_state++) {
        for (new_state=0; new_state<kLen; new_state++) {
          cls = classMat[dati][new_state]; /* do not subtract 1 */
          if (cls) {
            num  = forwardMat[tm1][init_state]*backwardMat[time][new_state]*trans[init_state][new_state];
            num *= cons;
            frac = num/denom;
            trans_tmp[init_state][new_state] += frac*fvi; 
          }
        }
      }
    }  
  }   

  normalizeRows(trans_tmp, kLen, kLen, trans);

}

static void calcStayer(mystr, pi_0, pi_1)
MYSTR *mystr;
double *pi_0, *pi_1;
{
  double pi0, pi1, sum0, sum1, like;
  int i, fvi, *freqVec=mystr->freqVec, freqVecSum=mystr->freqVecSum;
  char *rowMaxEq0=mystr->data_rowMaxEqual0, *rowMinEq1=mystr->data_rowMinEqual1, maxEq0, minEq1;

  pi0         = *pi_0;
  pi1         = *pi_1;
  sum0        = 0.0;
  sum1        = 0.0;
  mystr->pi_0 = pi0;
  mystr->pi_1 = pi1;

  for (i=0; i<mystr->data_nr; i++) {
    maxEq0  = rowMaxEq0[i];
    minEq1  = rowMinEq1[i];
    if (maxEq0 || minEq1) {
      fvi     = freqVec[i];
      like    = calcIndLikelihood(mystr, i);
      if (maxEq0) {
        sum0 += (pi0/like)*fvi;
      } else if (minEq1) {
        sum1 += (pi1/like)*fvi;
      } 
    } 

  }

  *pi_0 = sum0/((double) freqVecSum);
  *pi_1 = sum1/((double) freqVecSum);

}

static void getRowPtrs(rowMatAsVec, nr, nc, ret)
int *rowMatAsVec, **ret;
int nr, nc;
{
  int i;
  
  for (i=0; i<nr; i++) ret[i] = &rowMatAsVec[i*nc];

} 

static void setMinMaxVecs(mystr)
MYSTR *mystr;
{
  int i, nr=mystr->data_nr, *datarow, **rowPtr=mystr->data_rowPtr;
  int tmp1, tmp2, nc=mystr->data_nc;
  char *v0=mystr->data_rowMaxEqual0, *v1=mystr->data_rowMinEqual1;

  for (i=0; i<nr; i++) {
    datarow = rowPtr[i];
    tmp1    = maxiVecMiss(datarow, nc);
    if (tmp1 == 0) {
      v0[i] = 1;
    } else {
      v0[i] = 0;
    }
    tmp2    = miniVecMiss(datarow, nc);

    if (tmp2 == 1) {
      v1[i] = 1;
    } else {
      v1[i] = 0;
    }
  }

}

static void setClassMat(mystr)
MYSTR *mystr;
{
  int i, j, kLen=mystr->kLen, **classMat=mystr->classMat;
  
  /* classMat was initialized to -9999 */
  for (i=0; i<kLen; i++) {
    for (j=0; j<kLen; j++) {
      if (i == j) {
        classMat[i][j] = 1;
      } else {
        classMat[i][j] = 0;
      }
    }
  }
  i = MISSING;
  for (j=0; j<kLen; j++) classMat[i][j] = 1;

}


static void mystr_init(mystr, iargs, dargs, data, freqVec, k)
MYSTR *mystr;
int *iargs, *k, *data, *freqVec;
double *dargs;
{
  int debug, klen, nr, len;
  
  debug               = iargs[IARG_DEBUG];
  if (debug) Rprintf("Begin: mystr_init\n");

  if (iargs[IARG_kLen] >= MISSING) error("ERROR: Reset MISSING macro");

  mystr->converged    = 0;
  mystr->DEBUG        = debug;
  mystr->data_nr      = iargs[IARG_data_nr];
  mystr->data_nc      = iargs[IARG_data_nc];
  mystr->freqVecLen   = iargs[IARG_freqVecLen];
  mystr->maxiter      = iargs[IARG_maxiter];
  mystr->print        = iargs[IARG_print];
  mystr->kLen         = iargs[IARG_kLen];

  mystr->r_0          = dargs[DARG_r_0];
  mystr->p_00         = dargs[DARG_p_00];
  mystr->p_11         = dargs[DARG_p_11];
  mystr->pi_0         = dargs[DARG_pi_0];
  mystr->pi_1         = dargs[DARG_pi_1];
  mystr->epsilon      = dargs[DARG_epsilon];
  mystr->machine_eps  = dargs[DARG_machine_eps];

  mystr->data         = data;
  mystr->freqVec      = freqVec;
  mystr->freqVecSum   = sumivec(freqVec, mystr->data_nr);
  mystr->k            = k;

  /* Get the row pointers to data */
  mystr->data_rowPtr = ptriVec_alloc(mystr->data_nr, 1);
  getRowPtrs(data, mystr->data_nr, mystr->data_nc, mystr->data_rowPtr);

  klen               = mystr->kLen;
  mystr->trans       = dMat_alloc(klen, klen, 0, 0.0);
  mystr->trans_tmp   = dMat_alloc(klen, klen, 0, 0.0);
  mystr->forwardMat  = dMat_alloc(mystr->data_nc, klen, 0, 0.0);
  mystr->backwardMat = dMat_alloc(mystr->data_nc, klen, 0, 0.0);

  nr = mystr->data_nr;
  mystr->data_rowMaxEqual0 = cVec_alloc(nr, 0, 0);
  mystr->data_rowMinEqual1 = cVec_alloc(nr, 0, 0);
  setMinMaxVecs(mystr);

  mystr->tmpVec  = dVec_alloc(klen, 0, 0);
  mystr->tmpVec2 = dVec_alloc(klen, 0, 0);
  len            = MYMAX(mystr->data_nc, klen);
  mystr->tmpMat  = dMat_alloc(len, klen, 0, 0);

  mystr->classMat  = iMat_alloc(MISSING+1, klen, 1, -9999);
  setClassMat(mystr);
  
  if (debug) Rprintf("End: mystr_init\n");

}

static void mystr_free(mystr)
MYSTR *mystr;
{
  int len;

  if (mystr->DEBUG) Rprintf("Begin: mystr_free\n");

  if (mystr->data_rowPtr) free(mystr->data_rowPtr);
  if (mystr->trans) matrix_free((void **) mystr->trans, mystr->kLen);
  if (mystr->trans_tmp) matrix_free((void **) mystr->trans_tmp, mystr->kLen);
  if (mystr->forwardMat) matrix_free((void **) mystr->forwardMat, mystr->data_nc);
  if (mystr->backwardMat) matrix_free((void **) mystr->backwardMat, mystr->data_nc);
  if (mystr->data_rowMaxEqual0) free(mystr->data_rowMaxEqual0);
  if (mystr->data_rowMinEqual1) free(mystr->data_rowMinEqual1);
  if (mystr->tmpVec) free(mystr->tmpVec);
  if (mystr->tmpVec2) free(mystr->tmpVec2);
  len = MYMAX(mystr->data_nc, mystr->kLen);
  if (mystr->tmpMat) matrix_free((void **) mystr->tmpMat, len);
  len = MISSING + 1;
  if (mystr->classMat) matrix_free((void **) mystr->classMat, len);

  if (mystr->DEBUG) Rprintf("End: mystr_free\n");

} 
 
static int checkStop(ll0, ll1, iter, mystr)
double ll0, ll1;
int iter;
MYSTR *mystr;
{
  if (mystr->DEBUG) Rprintf("Begin: checkStop\n");

  double v1;
  int ret=0;

  v1  = fabs(ll0 - ll1);
  if (v1 < mystr->epsilon) ret=1;
  
  if (mystr->print > 1) Rprintf("Iter=%d diff=%g\n", iter, v1);
  
  if (mystr->DEBUG) Rprintf("End: checkStop\n");

  return(ret);

} /* END: checkStop */

static void EM_alg(mystr)
MYSTR *mystr;
{
  double loglike0, loglike1, *init_prob=mystr->init_prob;
  double **trans=mystr->trans;
  int iter=1, conv=0;

  if (mystr->DEBUG) Rprintf("Begin: EM_alg\n");
 
  init_prob[0] = mystr->r_0;
  init_prob[1] = 1.0 - init_prob[0];
  trans[0][0]  = mystr->p_00;
  trans[1][1]  = mystr->p_11;
  trans[0][1]  = 1.0 - trans[0][0];
  trans[1][0]  = 1.0 - trans[1][1];

  loglike0 = calcLikelihoodPattern(mystr);
  calcInitialPattern(mystr);
  calcTransitionPattern(mystr);
  calcStayer(mystr, &(mystr->pi_0), &(mystr->pi_1));
  loglike1 = calcLikelihoodPattern(mystr);
  if (checkStop(loglike0, loglike1, iter, mystr)) return;

  for (iter=2; iter<=mystr->maxiter; iter++) {
    loglike0 = loglike1;
    calcInitialPattern(mystr);
    calcTransitionPattern(mystr);
    calcStayer(mystr, &(mystr->pi_0), &(mystr->pi_1));
    loglike1 = calcLikelihoodPattern(mystr);
    if (checkStop(loglike0, loglike1, iter, mystr)) {
      conv = 1;
      break;
    }
  }
  mystr->converged = conv;
  if (mystr->print) {
    if (conv) {
      Rprintf("EM algorithm converged in %d iterations\n", iter);
    } else {
      Rprintf("EM algorithm did not converge after %d iterations\n", mystr->maxiter);
    }
  }

  if (mystr->DEBUG) Rprintf("End: EM_alg\n");

} /* END: EM_alg */


void MY_C_EM(iargs, dargs, data, freqVec, kVec, ret_conv, ret_init, ret_trans, ret_pi0, ret_pi1)
int *iargs, *kVec, *ret_conv, *data, *freqVec;
double *dargs, *ret_init, *ret_trans, *ret_pi0, *ret_pi1;
{
  MYSTR mystr;

  mystr_init(&mystr, iargs, dargs, data, freqVec, kVec);
  mystr.init_prob = ret_init;
  
  EM_alg(&mystr);

  /* return values */
  *ret_conv = mystr.converged;
  *ret_pi0  = mystr.pi_0;
  *ret_pi1  = mystr.pi_1;
  putMatinRowVec(mystr.trans, mystr.kLen, mystr.kLen, ret_trans);

  mystr_free(&mystr);

  return;

}  /* END: C_EM_ */

static int * localDecodeInd(mystr, ind, retptr)
MYSTR *mystr;
int ind, *retptr;
{
  int nc=mystr->data_nc, random_draw=mystr->random_draw, i, kLen=mystr->kLen, ival;
  double **forwardMat=mystr->forwardMat, **backwardMat=mystr->backwardMat, tmp, *fvec, *bvec, dprod;
  double pi_0_current, pi_1_current, pi_0=mystr->pi_0, pi_1=mystr->pi_1, denom, zero_prob, one_prob;
  char *rowMaxEq0=mystr->data_rowMaxEqual0, *rowMinEq1=mystr->data_rowMinEqual1;
  
  forwardIter(mystr, nc-1, ind, forwardMat);
  backwardIter(mystr, 0, ind, backwardMat);

  if (rowMinEq1[ind]) {
    pi_0_current = 0.0;
    pi_1_current = pi_1;
  } else if (rowMaxEq0[ind]) {
    pi_0_current = pi_0;
    pi_1_current = 0.0;
  } else {
    pi_0_current = 0.0;
    pi_1_current = 0.0;
  }

  for (i=0; i<nc; i++) {
    tmp      = 1.0 - pi_0_current - pi_1_current;
    fvec     = forwardMat[i];
    bvec     = backwardMat[i];
    dprod    = dotprod(fvec, bvec, kLen);
    denom    = pi_0_current + pi_1_current + tmp*dprod;

    one_prob = (tmp*(fvec[1]*bvec[1]) + pi_1_current)/denom;
    if (random_draw) {
      ival = rbinom(1, one_prob);
    } else {
      zero_prob = (tmp*(fvec[0] *bvec[0]) + pi_0_current)/denom;
      if (zero_prob > one_prob) {
        ival = 0;
      } else {
        ival = 1;
      }
    }
    *retptr = ival;
    retptr++;
  }

  return(retptr);

} /* END: localDecodeInd */

static void localDecode(mystr)
MYSTR *mystr;
{
  int i, *retptr=mystr->localDecodeVec;

  for (i=0; i<mystr->data_nr; i++) {
    retptr = localDecodeInd(mystr, i, retptr);
  }

} /* END: localDecode */



void MY_C_LocalDecode(iargs, dargs, data, freqVec, kVec, random_draw, init_prob, transVec, ret_vec)
int *iargs, *kVec, *data, *freqVec, *ret_vec, *random_draw;
double *dargs, *init_prob, *transVec;
{
  MYSTR mystr;

  /* For random number generation */
  GetRNGstate();

  mystr_init(&mystr, iargs, dargs, data, freqVec, kVec);
  mystr.localDecodeVec = ret_vec;
  mystr.random_draw    = *random_draw;
  mystr.init_prob      = init_prob;
  fillMat(transVec, mystr.kLen, mystr.kLen, mystr.trans);

  localDecode(&mystr);

  /* return values, ret_vec set above */

  mystr_free(&mystr);
  PutRNGstate();

  return;

}  /* END: MY_C_LocalDecode */



