/*
 * Adapted from Simon D. Levy's TinyEKF
 */

#include <stdio.h>
#include <stdlib.h>
#include "kalman_struct.h" //include kalman struct datatype

// Support both Arduino and command-line versions
#ifndef MAIN
extern "C" {
#endif

    void kalman_init(void *, int, int);
    int kalman_step(void *, double *);

#ifndef MAIN
}
#endif

/**
 * A header-only class for the Extended Kalman Filter.  Your implementing class should #define the constant N and 
 * and then #include <TinyEKF.h>  You will also need to implement a model() method for your application.
 */
class Kalman {

    private:

        kalman_t ftr; // struct that encapsulates all the necessary information
                         // hidden from the user for safety purpose

    protected:

        /**
          * The current state.
          */
        double * x;

        /**
         * Initializes a TinyEKF object.
         */
        Kalman() { 
            kalman_init(&this->ftr, Nsta, Mobs); 
            this->x = this->ftr.x; 
        }

        /**
         * Deallocates memory for a TinyEKF object.
         */
        ~Kalman() { }

        /**
         * Implement this function for your linear Kalman model.
         * @param F is the state transition matrix
         * @param H is the observation transformation matrix
         * @param w is the noise vector associated with state transition
         * @param v is the noise vector associated with observation
         */
        virtual void model(double F[Nsta][Nsta], double H[Mobs][Nsta], double w[Nsta], double v[Mobs]) = 0;

        /**
         * Sets the specified value of the prediction error covariance. <i>P<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setP(int i, int j, double value) 
        { 
            this->ftr.P[i][j] = value; 
        }

        /**
         * Sets the specified value of the process noise covariance. <i>Q<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setQ(int i, int j, double value) 
        { 
            this->ftr.Q[i][j] = value; 
        }

        /**
         * Sets the specified value of the observation noise covariance. <i>R<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setR(int i, int j, double value) 
        { 
            this->ftr.R[i][j] = value; 
        }

        /*
        //default zeros
        void setB(int i, int j, double value) 
        { 
            this->ftr.R[i][j] = value; 
        }
        */

    public:

        /**
         * Returns the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @return state value at index
         */
        double getX(int i) 
        { 
            return this->ftr.x[i]; 
        }

        /**
         * Sets the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @param value value to set
         */
        void setX(int i, double value) 
        { 
            this->ftr.x[i] = value; 
        }

        /**
          Performs one step of the prediction and update.
         * @param z observation vector, length <i>m</i>
         * @return true on success, false on failure caused by non-positive-definite matrix.
         */
        bool step(double * z) 
        { 
            this->model(this->ftr.F, this->ftr.H, this->ftr.w, this->ftr.v); 

            return kalman_step(&this->ftr, z) ? false : true; // return 1 indicates
                                                           // failed inversion
        }
};