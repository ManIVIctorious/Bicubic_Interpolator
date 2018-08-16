
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

int BicubicInterpolation(double* *v, int* nq, double h, int n_interpoints){

    int i, j, k, l, m, n;
    int n_points_new;
    double sum;
    double base_displacement;
    double * invM  = NULL;
    double * v_new = NULL;
    double * f     = NULL;
    double * x     = NULL;

// number of resulting data points and their base displacement
    n_points_new = ((nq[0] - 1) * (n_interpoints + 1) + 1) * ((nq[1] - 1) * (n_interpoints + 1) + 1);
    base_displacement = 1.0 / (n_interpoints + 1);

// allocate memory for inverse of M
    invM = calloc(16*16, sizeof(double));

// Fill matrix invM
//{{{
    invM[  0] =  1.0;
    invM[ 20] =  1.0;
    invM[ 32] = -3.0;
    invM[ 33] =  3.0;
    invM[ 36] = -2.0;
    invM[ 37] = -1.0;
    invM[ 48] =  2.0;
    invM[ 49] = -2.0;
    invM[ 52] =  1.0;
    invM[ 53] =  1.0;
    invM[ 72] =  1.0;
    invM[ 92] =  1.0;
    invM[104] = -3.0;
    invM[105] =  3.0;
    invM[108] = -2.0;
    invM[109] = -1.0;
    invM[120] =  2.0;
    invM[121] = -2.0;
    invM[124] =  1.0;
    invM[125] =  1.0;
    invM[128] = -3.0;
    invM[130] =  3.0;
    invM[136] = -2.0;
    invM[138] = -1.0;
    invM[148] = -3.0;
    invM[150] =  3.0;
    invM[150] =  3.0;
    invM[156] = -2.0;
    invM[158] = -1.0;
    invM[160] =  9.0;
    invM[161] = -9.0;
    invM[162] = -9.0;
    invM[163] =  9.0;
    invM[164] =  6.0;
    invM[165] =  3.0;
    invM[166] = -6.0;
    invM[167] = -3.0;
    invM[168] =  6.0;
    invM[169] = -6.0;
    invM[170] =  3.0;
    invM[171] = -3.0;
    invM[172] =  4.0;
    invM[173] =  2.0;
    invM[174] =  2.0;
    invM[175] =  1.0;
    invM[176] = -6.0;
    invM[177] =  6.0;
    invM[178] =  6.0;
    invM[179] = -6.0;
    invM[180] = -3.0;
    invM[181] = -3.0;
    invM[182] =  3.0;
    invM[183] =  3.0;
    invM[184] = -4.0;
    invM[185] =  4.0;
    invM[186] = -2.0;
    invM[187] =  2.0;
    invM[188] = -2.0;
    invM[189] = -2.0;
    invM[190] = -1.0;
    invM[191] = -1.0;
    invM[192] =  2.0;
    invM[194] = -2.0;
    invM[200] =  1.0;
    invM[202] =  1.0;
    invM[212] =  2.0;
    invM[214] = -2.0;
    invM[220] =  1.0;
    invM[222] =  1.0;
    invM[222] =  1.0;
    invM[224] = -6.0;
    invM[225] =  6.0;
    invM[226] =  6.0;
    invM[227] = -6.0;
    invM[228] = -4.0;
    invM[229] = -2.0;
    invM[230] =  4.0;
    invM[231] =  2.0;
    invM[232] = -3.0;
    invM[233] =  3.0;
    invM[234] = -3.0;
    invM[235] =  3.0;
    invM[236] = -2.0;
    invM[237] = -1.0;
    invM[238] = -2.0;
    invM[239] = -1.0;
    invM[240] =  4.0;
    invM[241] = -4.0;
    invM[242] = -4.0;
    invM[243] =  4.0;
    invM[244] =  2.0;
    invM[245] =  2.0;
    invM[246] = -2.0;
    invM[247] = -2.0;
    invM[248] =  2.0;
    invM[249] = -2.0;
    invM[250] =  2.0;
    invM[251] = -2.0;
    invM[252] =  1.0;
    invM[253] =  1.0;
    invM[254] =  1.0;
    invM[255] =  1.0;
//}}}


// allocate memory for new value array and
//  the f and x arrays required for the calculation of the interpolation functions
    v_new = malloc(n_points_new * sizeof(double));  // set to (*v), old (*v) is freed
    f = malloc(16 * sizeof(double));                // freed
    x = malloc(16 * sizeof(double));                // freed


// calculate first derivatives
    double *fx  = calloc(nq[0]*nq[1], sizeof(double));  // freed
    double *fy  = calloc(nq[0]*nq[1], sizeof(double));  // freed
    double *fxy = calloc(nq[0]*nq[1], sizeof(double));  // freed


// calculate first x derivatives and store them in fx
    for(i = 0; i < nq[0]; ++i){
        for(j = 1; j < (nq[1]-1); ++j){
            fx[i*nq[1] + j] = ( (*v)[i*nq[1] + (j+1)] - (*v)[i*nq[1] + (j-1)] ) / (2*h);
        }
    }

// calculate first y derivatives and store them in fy
    for(i = 1; i < (nq[0]-1); ++i){
        for(j = 0; j < nq[1]; ++j){
            fy[i*nq[1] + j] = ((*v)[(i+1)*nq[1] + j] - (*v)[(i-1)*nq[1] + j]) / (2*h);
        }
    }

// calculate cross xy derivatives and store them in fxy
    for(i = 1; i < (nq[0]-1); ++i){
        for(j = 1; j < (nq[1]-1); ++j){
            fxy[i*nq[1] + j] =  ( (*v)[(i-1)*nq[1] + (j-1)] - (*v)[(i-1)*nq[1] + (j+1)] - (*v)[(i+1)*nq[1] + (j-1)] + (*v)[(i+1)*nq[1] + (j+1)] ) / (4*h*h);
        }
    }


    for(i = 0; i < (nq[0]-1); ++i){
        for(j = 0; j < (nq[1]-1); ++j){

        // Matrix multiplication:   invM * f = x  (m×n * n×1 = m×1)
            for(m = 0; m < 16; ++m){
                x[m] = 0.0;

                x[m] += invM[m*16 +  0] * (*v)[  i  *nq[0] +   j  ];
                x[m] += invM[m*16 +  1] * (*v)[  i  *nq[0] + (j+1)];
                x[m] += invM[m*16 +  2] * (*v)[(i+1)*nq[0] +   j  ];
                x[m] += invM[m*16 +  3] * (*v)[(i+1)*nq[0] + (j+1)];

                x[m] += invM[m*16 +  4] *   fx[  i  *nq[0] +   j  ];
                x[m] += invM[m*16 +  5] *   fx[  i  *nq[0] + (j+1)];
                x[m] += invM[m*16 +  6] *   fx[(i+1)*nq[0] +   j  ];
                x[m] += invM[m*16 +  7] *   fx[(i+1)*nq[0] + (j+1)];

                x[m] += invM[m*16 +  8] *   fy[  i  *nq[0] +   j  ];
                x[m] += invM[m*16 +  9] *   fy[  i  *nq[0] + (j+1)];
                x[m] += invM[m*16 + 10] *   fy[(i+1)*nq[0] +   j  ];
                x[m] += invM[m*16 + 11] *   fy[(i+1)*nq[0] + (j+1)];

                x[m] += invM[m*16 + 12] *  fxy[  i  *nq[0] +   j  ];
                x[m] += invM[m*16 + 13] *  fxy[  i  *nq[0] + (j+1)];
                x[m] += invM[m*16 + 14] *  fxy[(i+1)*nq[0] +   j  ];
                x[m] += invM[m*16 + 15] *  fxy[(i+1)*nq[0] + (j+1)];
            }

        // calculate bicubic polynomial and store it to its appropriate position in v_new
        //  v_old[i*nq[1]+j] = v_new[i*(n_interpoints+1) * ((nq[1]-1)*(n_interpoints+1)+1) + j*(n_interpoints+1)]
            for(m = 0; m < (n_interpoints+2); ++m){
                for(n = 0; n < (n_interpoints+2); ++n){

                    for(k = 0, sum = 0.0; k < 4; ++k){
                        for(l = 0; l < 4; ++l){
                            sum += x[k*4 + l] * pow(base_displacement*n, l) * pow(base_displacement*m, k);
                        }
                    }

                    v_new[i*(n_interpoints+1) * ((nq[1]-1)*(n_interpoints+1)+1) + j*(n_interpoints+1) + m*((nq[1]-1)*(n_interpoints+1)+1) + n] = sum;
                }
            }
        }
    }

// set all points from old v to their appropriate position in v_new
    for(i = 0; i < nq[0]; ++i){
        for(j = 0; j < nq[1]; ++j){
            v_new[ i*(n_interpoints+1) * ((nq[1]-1)*(n_interpoints+1)+1) + j*(n_interpoints+1) ] = (*v)[i*nq[1] + j];
        }
    }

// free old v and point it to new one
    free((*v)); (*v) = v_new;
    free(fx);    fx  = NULL;
    free(fy);    fy  = NULL;
    free(fxy);   fxy = NULL;

    free(f);      f  = NULL;
    free(x);      x  = NULL;

    return n_points_new;
}

int InputFunction(char* inputfile, double** *q, int* nq, double* *v, double** *mu, int dimension, int dipole_flag);
int BicubicInterpolation(double* *v, int* nq, double h, int n_interpoints);

int main(int argc, char **argv){
    
    if(argc < 1){
        fprintf(stderr,
            "\n (-) Error: Please specify input file"
            "\n     Aborting...\n\n"
        );
        exit(1);
    }

    char * inputfile = argv[1];
    int n_inter = 1;

    int i, j;
    int    * nq = calloc(2, sizeof(int));
    double *  v = malloc(sizeof(double));
    double ** q = malloc(2 * sizeof(double*));
              q[0] = malloc(sizeof(double));
              q[1] = malloc(sizeof(double));


// input file
    InputFunction(inputfile, &q, nq, &v, NULL, 2, 0);
    double dq = q[1][1] - q[1][0];


// Print old values
//  to standard error
    fprintf(stderr, "\nOld Values:\n");
    for(i = 0; i < nq[0]*nq[1]; ++i){
        if(i % nq[1] == 0) fprintf(stderr, "\n");
        fprintf(stderr, "\t%2d  |  % lf\t% lf\t% lf\n", i, q[0][i], q[1][i], v[i]);
    }

//  to standard out    
    fprintf(stdout, "\nOld Values:\n");
    for(i = 0; i < nq[0]; ++i){
        for(j = 0; j < nq[1]; ++j){
            fprintf(stdout, "\t% 7.3lf", v[i*nq[1] + j]);
        }
        fprintf(stdout, "\n");
    }


// interpolation routine
    BicubicInterpolation(&v, nq, dq, n_inter);


// Print new values
//  to standard error
    fprintf(stderr, "\nNew Values:\n");
    for(i = 0; i < ((nq[0] - 1) * (n_inter + 1) + 1); ++i){
        for(j = 0; j < ((nq[1] - 1) * (n_inter + 1) + 1); ++j){
            fprintf(stderr, "\t%2d", i*((nq[1] - 1) * (n_inter + 1) + 1) + j);
            fprintf(stderr, "  |  ");
            fprintf(stderr, "\t% 7.3lf", dq*i/(n_inter+1));
            fprintf(stderr, "\t% 7.3lf", dq*j/(n_inter+1));
            fprintf(stderr, "  |  ");
            fprintf(stderr, "\t% 7.3lf", v[i*((nq[1] - 1) * (n_inter + 1) + 1) + j]);
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
    }

//  to standard out
    fprintf(stdout, "\nNew Values:\n");
    for(i = 0; i < ((nq[0] - 1) * (n_inter + 1) + 1); ++i){
        for(j = 0; j < ((nq[1] - 1) * (n_inter + 1) + 1); ++j){
            fprintf(stdout, "\t% 7.3lf", v[i*((nq[1] - 1) * (n_inter + 1) + 1) + j]);
        }
        fprintf(stdout, "\n");
    }

    return 0;
}
