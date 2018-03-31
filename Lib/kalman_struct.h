/* 
 * data structure for kalman filter
 * include this file after using #define for Nsta and Mobs
 */

typedef struct {
    int n;
    int m;

    double x[Nsta]; // state vector, size Nsta * 1
    double F[Nsta][Nsta]; // state transition matrix, size Nsta * Nsta
    double Ft[Nsta][Nsta]; // transpose of F, size Nsta * Nsta
    double w[Nsta]; // noise vector, size Nsta * 1
    double P[Nsta][Nsta]; // prediction covariance matrix, size Nsta * Nsta
    double Q[Nsta][Nsta]; // covariance matrix associated with w, size Nsta * Nsta

    double z[Mobs]; // observation vector, size Mobs * 1
    double v[Mobs]; // obervation noise vector, size Mobs * 1
    double H[Mobs][Nsta]; // observation transformation matrix, size Mobs * Nsta
    double Ht[Nsta][Mobs]; // transpose of H, size Nsta * Mobs
    double R[Mobs][Mobs]; // observation covarianc ematrix, size Mobs * Mobs
    double K[Nsta][Mobs]; // Kalman gain, size Nsta * Mobs

    /* temporary storage */
    double tmp0[Nsta]; // For Eq.1, size Nsta * 1
    double tmp1[Nsta][Nsta]; // For storing P*Ft, size Nsta * Nsta
    double tmp2[Nsta][Nsta]; // For storing Eq.2, size Nsta * Nsta
    double tmp3[Nsta][Mobs]; // For storing P*Ht, size Nsta * Mobs
    double tmp4[Mobs][Mobs]; // For storing H*P*Ht, size Mobs * Mobs
    double tmp5[Mobs][Mobs]; // For storing inverse of H*P*Ht, size Mobs * Mobs
    double tmp6[Mobs]; // For storing H*x, size Mobs * 1
    double tmp7[Mobs][Nsta]; // For storing H*P, size Mobs * Nsta

} kalman_t;



