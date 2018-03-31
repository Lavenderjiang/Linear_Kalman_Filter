/*
 *
 * Adapted from Simon D. Levy's TinyEKF
 * https://github.com/simondlevy/TinyEKF
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Cholesky-decomposition matrix-inversion code, adapated from
   http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */


static int choldc1(double * a, double * p, int n) {
    int i,j,k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i*n+j];
            for (k = i - 1; k >= 0; k--) {
                sum -= a[i*n+k] * a[j*n+k];
            }
            if (i == j) {
                if (sum <= 0) {
                    return 1; 
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j*n+i] = sum / p[i];
            }
        }
    }

    return 0; 
}


static int choldcsl(double * A, double * a, double * p, int n) 
{
    int i,j,k; double sum;
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            a[i*n+j] = A[i*n+j];
    if (choldc1(a, p, n)) return 1;
    for (i = 0; i < n; i++) {
        a[i*n+i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j*n+k] * a[k*n+i];
            }
            a[j*n+i] = sum / p[j];
        }
    }

    return 0; 
}


//matrix inverse using cholesky decomposition, a is the output
static int cholsl(double * A, double * a, double * p, int n) 
{
    int i,j,k;
    if (choldcsl(A,a,p,n)) return 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i*n+j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i*n+i] *= a[i*n+i];
        for (k = i + 1; k < n; k++) {
            a[i*n+i] += a[k*n+i] * a[k*n+i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i*n+j] += a[k*n+i] * a[k*n+j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i*n+j] = a[j*n+i];
        }
    }

    return 0; 
}


/*
// From Rosetta's code
// https://rosettacode.org/wiki/Cholesky_decomposition
static void cholesky_inv(double * A, double * Ainv, int n) {

    for (int i = 0; i < n; i++)
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += Ainv[i * n + k] * Ainv[j * n + k];
            Ainv[i * n + j] = (i == j) ?
                           sqrt(A[i * n + i] - s) :
                           (1.0 / Ainv[j * n + j] * (A[i * n + j] - s));
        }
 
    return L;
} */

static void zeros(double * a, int m, int n)
{
    int j;
    for (j=0; j<m*n; ++j)
        a[j] = 0;
}

#ifdef DEBUG
static void dump(double * a, int m, int n, const char * fmt)
{
    int i,j;

    char f[100];
    sprintf(f, "%s ", fmt);
    for(i=0; i<m; ++i) {
        for(j=0; j<n; ++j)
            printf(f, a[i*n+j]);
        printf("\n");
    }
}
#endif

/* C <- A * B */
static void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols)
{
    int i, j,l;

    for(i=0; i<arows; ++i)
        for(j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for(l=0; l<acols; ++l)
                c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
        }
}

static void mulvec(double * a, double * x, double * y, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i) {
        y[i] = 0;
        for(j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

static void transpose(double * a, double * at, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

/* A <- A + B */
static void accum(double * a, double * b, int m, int n)
{        
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] += b[i*n+j];
}

/* C <- A + B */
static void add(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] + b[j];
}


/* C <- A - B */
static void sub(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] - b[j];
}

static void negate(double * a, int m, int n)
{        
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] = -a[i*n+j];
}

static void mat_addeye(double * a, int n)
{
    int i;
    for (i=0; i<n; ++i)
        a[i*n+i] += 1;
}

/* TinyEKF code ------------------------------------------------------------------- */

#include "cFuncs.h"

typedef struct {

    double * x; // state vector, size Nsta * 1
    double * F; // state transition matrix, size Nsta * Nsta
    double * Ft; // transpose of F, size Nsta * Nsta
    double * w; // noise vector, size Nsta * 1
    double * P; // prediction covariance matrix, size Nsta * Nsta
    double * Q; // covariance matrix associated with w, size Nsta * Nsta

    double * z; // observation vector, size Mobs * 1
    double * v; // obervation noise vector, size Mobs * 1
    double * H; // observation transformation matrix, size Mobs * Nsta
    double * Ht; // transpose of H, size Nsta * Mobs
    double * R; // observation covarianc ematrix, size Mobs * Mobs
    double * K; // Kalman gain, size Nsta * Mobs

    /* temporary storage */
    double * tmp0; // For Eq.1, size Nsta * 1
    double * tmp1; // For storing P*Ft, size Nsta * Nsta
    double * tmp2; // For storing Eq.2, size Nsta * Nsta
    double * tmp3; // For storing P*Ht, size Nsta * Mobs
    double * tmp4; // For storing H*P*Ht, size Mobs * Mobs
    double * tmp5; // For storing inverse of H*P*Ht, size Mobs * Mobs
    double * tmp6; // For storing H*x, size Mobs * 1
    double * tmp7; // For storing H*P, size Mobs * Nsta
    double * tmp8; // For matrix inverse

} kalman_t;


//original code took void* v as an argument
//use m,n to initialize the data sturcture
static void unpack(void* v, kalman_t* ftr, int n, int m)
{
    /* skip over two ints (n,m) in data structure */
    char* cptr = (char*)v;
    cptr += 2*sizeof(int);
    double * dptr = (double *)cptr;

    //allocate memory for matrices and vectors
    ftr->x = dptr;
    dptr += n;
    ftr->F = dptr;
    dptr += n*n;
    ftr->Ft = dptr;
    dptr += n*n;
    ftr->w = dptr;
    dptr += n;
    ftr->P = dptr;
    dptr += n*n;
    ftr->Q = dptr;
    dptr += n*n;

    ftr->z = dptr;
    dptr += m;
    ftr->v = dptr;
    dptr += m;
    ftr->H = dptr;
    dptr += m*n;
    ftr->Ht = dptr;
    dptr += n*m;
    ftr->R = dptr;
    dptr += m*m;
    ftr->K = dptr;
    dptr += n*m;

    ftr->tmp0 = dptr;
    dptr += n;
    ftr->tmp1 = dptr;
    dptr += n*n;
    ftr->tmp2 = dptr;
    dptr += n*n;
    ftr->tmp3 = dptr;
    dptr += n*m;
    ftr->tmp4 = dptr;
    dptr += m*m;
    ftr->tmp5 = dptr;
    dptr += m*m;
    ftr->tmp6 = dptr;
    dptr += m;
    ftr->tmp7 = dptr;
    dptr += m*n;
    ftr->tmp8; //for matrix inverse
}

void kalman_init(void* v, int n, int m)
{
  /* retrieve n, m and set them in incoming data structure */
  int * ptr = (int *)v;
  *ptr = n;
  ptr++;
  *ptr = m;

  /* unpack rest of incoming structure for initlization */
  kalman_t ftr; //filter
  unpack(v, &ftr, n, m); 

  // zero-out matrices
  zeros(ftr.P,n,n);
  zeros(ftr.K,n,m);
  zeros(ftr.Q,n,n); 
}

int kalman_step(void* v, double* z)
{

  /* unpack incoming structure */

  int * ptr = (int *)v;
  int n = *ptr;
  ptr++;
  int m = *ptr;

  kalman_t ftr;
  unpack(v, &ftr, n, m);

  //Eq.1
  mulmat(ftr.F, ftr.x, ftr.tmp0, n, n, 1);
  add(ftr.tmp0, ftr.w, ftr.x, n);

  //Eq.2
  transpose(ftr.F, ftr.Ft, n, n);
  mulmat(ftr.P, ftr.Ft, ftr.tmp1, n, n, n);
  mulmat(ftr.F, ftr.tmp1, ftr.tmp2, n, n, n);
  add(ftr.tmp2, ftr.Q, ftr.P, n*n);

  //Eq.3 Observation vector z is updated by the argument
  /*
  mulmat(ftr.H, ftr.x, ftr.tmp6, m, n, 1);
  add(ftr.tmp6, ftr.v, ftr.z);*/
  
  //Eq.4
  transpose(ftr.H, ftr.Ht, m, n);
  mulmat(ftr.P, ftr.Ht, ftr.tmp3, n, n, m);
  mulmat(ftr.H, ftr.tmp3, ftr.tmp4, m, n, m);
  accum(ftr.tmp4, ftr.R, m, m);
  if (cholsl(ftr.tmp4, ftr.tmp5, ftr.tmp8, m) == 1) return 1;
  mulmat(ftr.Ht, ftr.tmp5, ftr.tmp3, n, m, m);
  mulmat(ftr.P,ftr.tmp3,ftr.K, n, n, m); //Kalman gain

  //Eq.5
  mulmat(ftr.H, ftr.x, ftr.tmp6, m, n, 1);
  negate(ftr.tmp6, m, 1);
  add(ftr.tmp6, z, z, m);
  mulmat(ftr.K, ftr.tmp6, ftr.tmp0, n, m, 1);
  accum(ftr.x, ftr.tmp0, n, 1);

  //Eq.6
  mulmat(ftr.H, ftr.P, ftr.tmp7, m, n, n);
  mulmat(ftr.K, ftr.tmp7, ftr.tmp1, n, m, n);
  negate(ftr.tmp1, n, n);
  accum(ftr.P, ftr.tmp1, n, n);
  return 0;

}

