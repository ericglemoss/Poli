#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define eps 10e-7

typedef struct Matrix { //Estrutura de dados criada para simular uma matriz
    int n, m;//Linhas e colunas da matriz, respectivamente
    double** v; //Ponteiro duplo que será utilizado para alocar dinamicamente uma matriz
} Matrix;

Matrix* createMatrix (int n, int m){ //Aloca dinamicamente uma matriz (ou seja, cria uma matriz)
    Matrix* x;

    x = malloc( ( sizeof( Matrix ) ) );
    x->v = malloc( sizeof (double *) * n);
    x->v[0] = calloc( sizeof(double), m*n);

    for (int i = 0; i < n; i++) {
        x->v[i] = x->v[0] + m * i;
    }

    x->n = n;
    x->m = m;

    return x;

}

void deleteMatrix (Matrix* M) { //Desaloca (deleta) uma matriz criada
    free(M->v[0]);
    free(M->v);
    free(M);
}

void printMatrix (Matrix* M) { //Imprime uma matriz
    for (int i = 0; i < M->n; i++) {
        for (int j = 0; j < M->m; j++){
            printf (" %8.3f", M->v[i][j]);
        }
        printf ("\n");
    }
    printf ("\n");
}

void matrixTransposte (Matrix* M) { //Calcula a tranposta de uma matriz
    for (int i = 0; i < M->n; i++){
        for (int j = 0; j < M->m; j++){
            double temp = M->v[i][j];
            M->v[i][j] = M->v[j][i];
            M->v[j][i] = temp;
        }
    }

}



void qrRotation ( double** v, int n, int m, int i, int j, int k ) { //Aplicação do método rotação de givens
    double tau, c, s, aux;

    if ( abs( v[i][k] ) > abs( v[j][k] ) ){
        tau = - ( v[j][k]/v[i][k] );
        c = 1/( pow(1 + tau*tau, 0.5) );
        s = tau*c;

    }

    else {
        tau = - ( v[i][k] / v[j][k] );
        s = 1/( pow(1 + tau*tau, 0.5) );
        c = tau*s;

    }

    for ( int counter = k; counter <= m - 1; k++ ){
        aux = c * v[i][k] - s * v[j][k];
        v[j][k] = s * v[i][k] + c * v[j][k];
        v[i][k] = aux;
    }

}


Matrix* solveSystem (Matrix* R, Matrix* b) { //Resolve um sistema ->exclusivamente<- triangular superior
    Matrix * X = createMatrix(R->n, 1);
    double sum;

    if (R->v[R->n - 1][R->m - 1]) != eps
        X[R->n - 1][1] = b[R->n - 1][1]/R->v[R->n - 1][R->m - 1];

    else
        return NULL; //Sistema indeterminado


    for ( int i = R->m - 1; i <= 0; i-- ) {
        sum = 0;
        for ( int j = R->m - 1; j <  )
    }




        }
    }


}

void Rot_givens ( double** W, int n, int m ) {




}
