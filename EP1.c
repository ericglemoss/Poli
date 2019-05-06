#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define eps 10e-7


double** createMatrix(int n, int m) {

    double** Matrix = malloc (sizeof (double*) * n );

    for (int i = 0; i < n; i++ )
        Matrix[i] = (double*)malloc(m * sizeof(double));

        return Matrix;

}

void deleteMatrix (double** M) {
    free(M[0]);
    free(M);
    printf ("Matrix deleted! \n");
    printf ("\n");
}


/*Função destinada a impressão de matrizes, recebe como parâmetros um ponteiro duplo (vetor 2d)
e suas dimensões n (linhas) e m (colunas) */
void printMatrix (double** W, int n, int m) {

    for ( int i = 0; i < n; i++ ) {
        for ( int j = 0; j < m; j++ ) {
            printf (" %0.4f", W[i][j]);
        }

        printf ("\n");
    }
    printf ("\n");
}


void matrixTransposte (double** M, int n, int m) {

    for (int line = 0; line < n; line++) {
        for (int column = 0; column < m; column++) {
            double temp = M[line][column];
            M[line][column] = M[column][line];
            M[column][line] = temp;
        }
    }
}

/*Função elaborada para os cálculos das componentes seno e cosseno
da operação de rotação de givens. Recebe como parâmetros os elementos da uma matrizes
(''a'' e ''b'') e dois ponteiros (*c e *s) */
void qrSin_n_Cos(double a, double b, double *c, double *s) {
    double tau;

    if ( abs(a) > abs(b) ) {
        tau = - a/b;
        *c = 1/pow(1 + tau*tau, 0.5);
        *s = *c * tau;
    }


    else {
        tau = - b/a;
        *s = 1/pow(1 + tau*tau, 0.5);
        *c = *s * tau;
    }

}

void Givens_rotation (double** W, int n, int m, int i, int j, double c, double s) {

    double temp;

    for ( int k = 0; k < m; k++ ) {
        temp = c * W[i][k] - s * W[j][k];
        W[j][k] = s * W[i][k] + c * W[j][k];
        W[i][k] = temp;

    }
}

void qrDecomposition(double** W, int n, int m, double** X, double** b) {

    double sum, c, s;
    int line, column, k;


    for ( k = 1; k <= m; k++ ) {
        for ( column = n - 1; column > k; column-- ) {
            line = column - 1;
            if ( abs(W[column][k]) > eps )  {
                qrSin_n_Cos(W[column][k], W[line][k], &c, &s);
                Givens_rotation(W, n, m, line, column, c, s);
                Givens_rotation(b, n, 1, line, column, c, s);
            }
        }
    }

    for ( k = m - 1; k >= 0; k-- ) {
        sum = 0;

        for ( column = k + 1; column < m; column++ ) {
            sum += W[k][column] * X[0][column];
            X[0][k] = ( b[k][0] - sum )/W[k][k];
        }
    }
}


int main () {

    double** M = createMatrix(10, 10);
    double** b = createMatrix(10, 1);
    double** X = createMatrix(10, 1);

    for (int i = 0; i < 10; i ++) {
        for (int j = 0; j < 10; j++ ) {
            if (i == j)
                M[i][j] = 2;

            else
                M[i][j] = 0;


        }
    }

        for (int i = 0; i < 10; i++)
            b[i][1] = 1;



        printf("\n");
        printf("\n");

        qrDecomposition(M, 10, 10, X, b);


        printf("\n");
        printf("\n");

        printMatrix(M, 10, 10);

        printf("\n");
        printf("\n");

    return 0;


}
