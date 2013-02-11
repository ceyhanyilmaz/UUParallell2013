#include <stdio.h>

typedef struct {
    int nrows, ncols;
    float *data;
} matrix_t;

float constC(int i, int j, matrix_t A, matrix_t B) {
    int k;
    float res = 0.0;
    for (k = 0; k < A.nrows; k++) {
        res += A.data[i + (k*A.nrows)]*B.data[j*B.ncols + k];
    }
    return res;
}

int main() {
    matrix_t A;
    A.nrows = 5;
    A.ncols = 5;
    A.data = (float *)malloc(A.ncols*A.nrows*sizeof(float));
    int k;
    
    for (k = 0; k < A.ncols*A.nrows; k++) {
        A.data[k] = k % 9;
    }   

    matrix_t B;
    B.nrows = 5;
    B.ncols = 5;
    B.data = (float *)malloc(B.ncols*B.nrows*sizeof(float));

    for (k = 0; k < B.ncols*B.nrows; k++) {
        B.data[k] = k % 7;
    }   

    matrix_t C;
    C.nrows = 5;
    C.ncols = 5;
    C.data = (float *)malloc(C.ncols*C.nrows*sizeof(float));

    int i, j;
    if (A.nrows != B.ncols) {
        return -1;
    }

    for (i = 0; i < A.nrows; i++) {
        for(j = 0; j < B.ncols; j++) {
            C.data[i + j*C.ncols] = constC(i,j,A,B);
            printf("C_(%d,%d) = %f\n", i,j,C.data[i + j*C.ncols]);
        }
    }
    
    return 0;
}