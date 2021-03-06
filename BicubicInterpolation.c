
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int BicubicInterpolation(double* *v, int* nq, double h, int n_interpoints);

int BicubicInterpolation(double* *v, int* nq, double h, int n_interpoints){

    int i, j, k, l, m, n;
    int n_points_new;
    double sum;
    double base_displacement;
    double * invM  = NULL;  // freed
    double *fx     = NULL;  // freed
    double *fy     = NULL;  // freed
    double *fxy    = NULL;  // freed
    double * x     = NULL;  // freed
    double * v_new = NULL;  // set to (*v), old (*v) is freed

// number of resulting data points and their base displacement
    n_points_new = ((nq[0] - 1) * (n_interpoints + 1) + 1) * ((nq[1] - 1) * (n_interpoints + 1) + 1);
    base_displacement = 1.0 / (n_interpoints + 1);

// allocate memory for inverse of M, the new value array and the x array
//  being the result of the equation system invM * [f(1:4); fx(1:4); fy(1:4); fxy(1:4)]
    invM = calloc(16*16, sizeof(double));
    v_new = malloc(n_points_new * sizeof(double));
    x = malloc(16 * sizeof(double));

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

// allocate memory for the x-derivative (fx), y-derivative (fy) and cross-derivative (fxy)
//  arrays, each containing the respecting derivative for every original point
    fx  = calloc(nq[0]*nq[1], sizeof(double));
    fy  = calloc(nq[0]*nq[1], sizeof(double));
    fxy = calloc(nq[0]*nq[1], sizeof(double));


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
            fxy[i*nq[1] + j] =  (
                  (*v)[(i-1)*nq[1] + (j-1)] - (*v)[(i-1)*nq[1] + (j+1)]
                - (*v)[(i+1)*nq[1] + (j-1)] + (*v)[(i+1)*nq[1] + (j+1)]
            ) / (4*h*h);
        }
    }

// above derivatives only fill the grid body, leaving either the top and bottom,
//  and/or right and left edges to be zero. By using forward and backward differences
//  these points can be approximated as well:

// calculate forward and backward x-derivatives for left and right edges, respectively
    for(i = 0; i < nq[0]; ++i){
        fx[i*nq[1]             ] = ( (*v)[i*nq[1] + 1        ] - (*v)[i*nq[1]            ] ) / h;
        fx[i*nq[1] + (nq[1]-1) ] = ( (*v)[i*nq[1] + (nq[1]-1)] - (*v)[i*nq[1] + (nq[1]-2)] ) / h;
    }


// calculate forward and backward y-derivatives for top and bottom edges, respectively
    for(j = 0; j < nq[1]; ++j){
        fy[                  j ] = ( (*v)[nq[1]           + j] - (*v)[                  j] ) / h;
        fy[(nq[0]-1)*nq[1] + j ] = ( (*v)[(nq[0]-1)*nq[1] + j] - (*v)[(nq[0]-2)*nq[1] + j] ) / h;
    }


// calculate forward and backward cross xy-derivatives
//  for left and right edges:   ([-1;0;1]/(2*dy)) * ([-1,1]/dx) = [1,-1; 0,0; -1,1] / (2*dx*dy)
    for(i = 1; i < (nq[0]-1); ++i){
        fxy[i*nq[1]            ] = (
                                     (*v)[(i-1)*nq[1]            ] - (*v)[(i-1)*nq[1] +     1    ]
                                   - (*v)[(i+1)*nq[1]            ] + (*v)[(i+1)*nq[1] +     1    ]
                                   ) / (2*h*h);

        fxy[i*nq[1] + (nq[1]-1)] = (
                                     (*v)[(i-1)*nq[1] + (nq[1]-2)] - (*v)[(i-1)*nq[1] + (nq[1]-1)]
                                   - (*v)[(i+1)*nq[1] + (nq[1]-2)] + (*v)[(i+1)*nq[1] + (nq[1]-1)]
                                   ) / (2*h*h);
    }

//  for top and bottom edges:   ([-1;1]/dy) * ([-1,0,1]/(2*dx)) = [1,0,-1; -1,0,1] / (2*dx*dy)
    for(j = 1; j < (nq[1]-1); ++j){
        fxy[                  j] = (
                                     (*v)[                  (j-1)] - (*v)[                  (j+1)]
                                   - (*v)[nq[1]           + (j-1)] + (*v)[nq[1]           + (j+1)]
                                   ) / (2*h*h);

        fxy[(nq[0]-1)*nq[1] + j] = (
                                     (*v)[(nq[0]-2)*nq[1] + (j-1)] - (*v)[(nq[0]-2)*nq[1] + (j+1)]
                                   - (*v)[(nq[0]-1)*nq[1] + (j-1)] + (*v)[(nq[0]-1)*nq[1] + (j+1)]
                                   ) / (2*h*h);
    }

// calculate cross xy-derivatives for the corner points
    fxy[    0          ] = (
                        (*v)[0]                 - (*v)[1]
                      - (*v)[nq[1]]             + (*v)[nq[1] + 1]
                    ) / (h*h);

    fxy[(nq[1]-1)      ] = (
                        (*v)[(nq[1]-2)]         - (*v)[(nq[1]-1)]
                      - (*v)[nq[1] + (nq[1]-2)] + (*v)[nq[1] + (nq[1]-1)]
                    ) / (h*h);

    fxy[(nq[0]-1)*nq[1]] = (
                        (*v)[(nq[0]-2)*nq[1]            ] - (*v)[(nq[0]-2)*nq[1] +     1    ]
                      - (*v)[(nq[0]-1)*nq[1]            ] + (*v)[(nq[0]-1)*nq[1] +     1    ]
                    ) / (h*h);

    fxy[nq[0]*nq[1]-1  ] = (
                        (*v)[(nq[0]-2)*nq[1] + (nq[1]-2)] - (*v)[(nq[0]-2)*nq[1] + (nq[1]-1)]
                      - (*v)[(nq[0]-1)*nq[1] + (nq[1]-2)] + (*v)[(nq[0]-1)*nq[1] + (nq[1]-1)]
                    ) / (h*h);



// start interpolation procedure
    for(i = 0; i < (nq[0]-1); ++i){
        for(j = 0; j < (nq[1]-1); ++j){

        // Matrix multiplication:   invM * f = x  (m??n * n??1 = m??1)
            for(m = 0; m < 16; ++m){
                x[m] = 0.0;

                x[m] += invM[m*16 +  0] * (*v)[  i  *nq[1] +   j  ];
                x[m] += invM[m*16 +  1] * (*v)[  i  *nq[1] + (j+1)];
                x[m] += invM[m*16 +  2] * (*v)[(i+1)*nq[1] +   j  ];
                x[m] += invM[m*16 +  3] * (*v)[(i+1)*nq[1] + (j+1)];

                x[m] += invM[m*16 +  4] *   fx[  i  *nq[1] +   j  ];
                x[m] += invM[m*16 +  5] *   fx[  i  *nq[1] + (j+1)];
                x[m] += invM[m*16 +  6] *   fx[(i+1)*nq[1] +   j  ];
                x[m] += invM[m*16 +  7] *   fx[(i+1)*nq[1] + (j+1)];

                x[m] += invM[m*16 +  8] *   fy[  i  *nq[1] +   j  ];
                x[m] += invM[m*16 +  9] *   fy[  i  *nq[1] + (j+1)];
                x[m] += invM[m*16 + 10] *   fy[(i+1)*nq[1] +   j  ];
                x[m] += invM[m*16 + 11] *   fy[(i+1)*nq[1] + (j+1)];

                x[m] += invM[m*16 + 12] *  fxy[  i  *nq[1] +   j  ];
                x[m] += invM[m*16 + 13] *  fxy[  i  *nq[1] + (j+1)];
                x[m] += invM[m*16 + 14] *  fxy[(i+1)*nq[1] +   j  ];
                x[m] += invM[m*16 + 15] *  fxy[(i+1)*nq[1] + (j+1)];
            }

        // calculate bicubic polynomial and store it to its appropriate position in v_new
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
//  v_old[i*nq[1]+j] = v_new[i*(n_interpoints+1) * ((nq[1]-1)*(n_interpoints+1)+1) + j*(n_interpoints+1)]
    for(i = 0; i < nq[0]; ++i){
        for(j = 0; j < nq[1]; ++j){
            v_new[ i*(n_interpoints+1) * ((nq[1]-1)*(n_interpoints+1)+1) + j*(n_interpoints+1) ] = (*v)[i*nq[1] + j];
        }
    }

// free all arrays and v_old, then point v to v_new
    free((*v)); (*v) = v_new;
    free(fx);    fx  = NULL;
    free(fy);    fy  = NULL;
    free(fxy);   fxy = NULL;
    free(x);      x  = NULL;
    free(invM); invM = NULL;

    return n_points_new;
}
