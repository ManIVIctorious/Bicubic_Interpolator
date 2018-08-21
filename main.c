
#include <stdio.h>
#include <stdlib.h>

int InputFunction(char* inputfile, double** *q, int* nq, double* *v, double** *mu, int dimension, int dipole_flag);
int BicubicInterpolation(double* *v, int* nq, double h, int n_interpoints);

int main(int argc, char **argv){

    if(argc != 2 && argc != 3){
        fprintf(stderr,
            "\nSynopsis:"
            "\n\t%s <input-file> [number-of-interpolation-points]"
            "\n\n", argv[0]
        );
        exit(0);
    }

    int    n_inter   = 1;
    char * inputfile = argv[1];

    if(argc == 3){
        n_inter = atoi(argv[2]);
    }


    int       i, j, nq[2];
    double    dq;
    double *  v = NULL;
    double ** q = NULL;


// file input
    v = malloc(    sizeof(double));
    q = malloc(2 * sizeof(double*));
    if(v == NULL || q == NULL){
        fprintf(stderr,
            "\n (-) Error in initial memory allocation of *v and/or *q"
            "\n     Aborting...\n\n"
        );
        exit(1);
    }
    q[0] = NULL;
    q[1] = NULL;

    InputFunction(inputfile, &q, nq, &v, NULL, 2, 0);
    dq = q[1][1] - q[1][0];


// run interpolation function
    BicubicInterpolation(&v, nq, dq, n_inter);


// print new values
    printf("N");
    printf("\t%d", ((nq[0] - 1) * (n_inter + 1) + 1));
    printf("\t%d", ((nq[1] - 1) * (n_inter + 1) + 1));
    printf("\n");

    for(i = 0; i < ((nq[0] - 1) * (n_inter + 1) + 1); ++i){
        for(j = 0; j < ((nq[1] - 1) * (n_inter + 1) + 1); ++j){

            printf("\t% 16.12lf", q[0][0] + dq/(n_inter+1)*i);
            printf("\t% 16.12lf", q[1][0] + dq/(n_inter+1)*j);
            printf("\t% 16.12lf", v[i*((nq[1] - 1) * (n_inter + 1) + 1) + j]);
            printf("\n");
        }
        printf("\n");
    }


// free memory
    free(v); v = NULL;
    free(q[0]); q[0] = NULL;
    free(q[1]); q[1] = NULL;
    free(q);    q    = NULL;

    return 0;
}
