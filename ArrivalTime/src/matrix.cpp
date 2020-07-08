#include <stdio.h>
#include "matrix.h"

void cov(float matrix[][CHS], int idx0, int idx1, float cov[CHS][CHS], float coff) {

    int i, j, k;
    float sum;

    /*
    for(i = 0; i < CHS; i++) {
        for(j = 0; j < CHS; j++) {
            cov[i][j] = 0.0;
        }
    }*/

    for(i = 0; i < CHS; i++) {
        for(j = 0; j < CHS; j++) {
            sum = 0.0;
            for(k = idx0; k < idx1; k++) {
                sum += matrix[k][i]*matrix[k][j];
            }
            cov[i][j] = sum/coff;
        }
    }

}
