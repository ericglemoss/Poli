#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define eps 10e-5

void NormalizeColumns (double** W, int n, int p) {

  double sum, S[p];

    for ( int j = 0 ; j < p; j++ ) {
    sum = 0;
      for ( int i = 0; i < n; i++ ) {
        sum += pow(W[i][j], 2);
      }
      if(sum > 0){
        for ( int k = 0; k < n; k++ ) {
            W[k][j] /= sqrt(sum);
        }
      }
  }

}

void DefineMax (double** H, int p, int m) {
  for(int i = 0; i < p; i++){
    for (int j = 0; j < m; j++){
      if (H[i][j] < 0)
        H[i][j] = 0;
    }
  }
}

void RandomInit (double ** W, int lines, int columns) {
  for ( int i = 0; i < lines; i++ ) {
      for ( int j = 0; j < columns; j++ ) {
          W[i][j] = (double)rand();
      }
  }

}

double** createMatrix(int n, int m) {

    double** Matrix = malloc (sizeof (double*) * n );

    for (int i = 0; i < n; i++ )
        Matrix[i] = (double*)malloc(m * sizeof(double));
    for (int i = 0; i < n; i++){
      for (int j = 0; j < m; j++){
        Matrix[i][j] = 0;
      }
    }
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

                    printf (" %8.4f", W[i][j]);
        }

        printf ("\n");
    }
    printf ("\n");
}

double** matrixTransposte (double** M, int n, int m) {

    double** T = createMatrix(m, n);

    for (int line = 0; line < n; line++) {
        for (int column = 0; column < m; column++) {
            T[column][line] = M[line][column];
        }
    }

    return T;
}

/*Função elaborada para os cálculos das componentes seno e cosseno
da operação de rotação de givens. Recebe como parâmetros os elementos da uma matrizes
(''a'' e ''b'') e dois ponteiros (*c e *s) */
void qrSin_n_Cos(double Wik, double Wjk, double *c, double *s) {
    double tau;

    if ( fabs(Wik) > fabs(Wjk) ) {
        tau = - Wjk/Wik;
        *c = 1/pow(1 + tau*tau, 0.5);
        *s = *c * tau;
    }


    else {
        tau = - Wik/Wjk;
        *s = 1/pow(1 + tau*tau, 0.5);
        *c = *s * tau;
    }
}

void Givens_rotation ( double** W, int n, int m, int i, int j, double c, double s) {

    double temp;

    for ( int k = 0; k < m; k++ ) {
        temp = c * W[i][k] - s * W[j][k];
        W[j][k] = s * W[i][k] + c * W[j][k];
        W[i][k] = temp;
    }
}

void qrDecomposition(double** W, int n, int m, double** X, double** b) {

    double sum, c, s, tau;
    int i, j, k;

    for ( k = 0; k < m ; k++ ) {
        for ( j = n-1; j > k ; j-- ) {
            i = j - 1;
            if ( fabs(W[j][k]) > eps ) {
                qrSin_n_Cos(W[i][k], W[j][k], &c, &s);
                Givens_rotation(W, n, m, i, j, c, s);
                Givens_rotation(b, n, 1, i, j, c, s);
            }
        }
    }

    for ( k = m - 1; k >= 0; k-- ) {
      sum = 0;
      for ( j = k; j < m; j++ )
        sum += W[k][j] * X[j][0];
      X[k][0] = (b[k][0] - sum )/W[k][k];
    }
}

double** solveMultipleSystems (double** W, int n, int p, int m, double** H, double** A) {
    double sum, c, s;

    for ( int k = 0; k < p; k++ ) {
        for ( int j = n - 1; j > k; j-- ) {
            int i = j - 1;
            if (W[j][k] != 0 ) {
              qrSin_n_Cos(W[i][k], W[j][k], &c, &s);
              Givens_rotation(W, n, p, i, j, c, s);
              Givens_rotation(A, n, m, i, j, c, s);
            }
        }
    }
    for ( int k = p - 1; k >= 0; k--  ) {
    for ( int j = 0; j < m; j++  ) {
      sum = 0;
      for (int i = k ; i < p ; i++){
        sum += W[k][i] * H[i][j];
      }
      H[k][j] = (A[k][j] - sum)/W[k][k];
    }
  }
}

void CreatecopiaA(double** B, double** A, int n,int m){
    for (int i = 0; i < n; i++){
      for (int j = 0; j < m; j++)
        B[i][j] = A[i][j];
    }
}
double** MatrixMultiplication (double** W, double** H, int n, int m, int p) {
  double** A = createMatrix(m, n);

     for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < m; j++ )
            for ( int k = 0; k < p; k++ )
                A[i][j] += W[i][k] * H[k][j];

  return A;
}

void MMQ_Alternados(double** W, double** A, double** H, int n, int p, int m){
    RandomInit(W, n, p);
    int itmax = 0;
    double sum;
    double* S;
    double Er = 1, E_atual = 1, E_ant = 0;
    double** B   = createMatrix(n, m);
    double** A_aprox = createMatrix(n,m);
    double** A_t = createMatrix(m, n);
    double** H_t = createMatrix(m, p);
    double** W_t = createMatrix(p, n);
    double soma_erro = 0;

    CreatecopiaA(B, A, n, m);
    while(itmax < 100 & Er > 0.000001){
      E_atual = 0;
      H = createMatrix(p, m);
      W_t = createMatrix(p, n);
      CreatecopiaA(A, B, n, m);
      itmax += 1;
      NormalizeColumns(W, n, p);
      solveMultipleSystems(W, n, p, m, H, A);
      DefineMax(H, p, m);
      H_t = matrixTransposte(H, p, m);
      A_t = matrixTransposte(B, n, m);
      solveMultipleSystems(H_t, m, p, n, W_t, A_t);
      W = matrixTransposte(W_t, p, n);
      DefineMax(W, n, p);
      A_aprox = MatrixMultiplication(W, H, n, m, p);
      printMatrix(A_aprox, n, m);
      for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++ )
            E_atual += pow(B[i][j] - A_aprox[i][j], 2 );
      }
      Er = fabs(E_atual-E_ant)/E_atual;
      E_ant = E_atual;
    }
    deleteMatrix(A_t);
    deleteMatrix(H_t);
    deleteMatrix(W_t);
    deleteMatrix(B);
    deleteMatrix(A_aprox);
}

void normalizarRGB(double **A, int n, int m){
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++)
        A[i][j] = A[i][j]/255;
    }
}

void Ler_arquivo(double **W, int n, int m, char *arquivo){
  FILE *fpointer = fopen(arquivo, "r");
  char p, q;

  for(int k = 0; k < n; k++){
      for(int j = 0; p != '\r' && p != '\n' && p != EOF; j++ ){
        fscanf(fpointer,"%d%c", &q, &p);
        if(j < m)
          W[k][j] = (double) q/255.0;
      }
      p = '\0';
  }
  fclose(fpointer);
}

int main () {

        int ndig_treino = 100;
        int p = 5;
        double** H = createMatrix(p, ndig_treino);
        double** A0 = createMatrix(784, ndig_treino);
        double** W0 = createMatrix(784, p);
        Ler_arquivo(A0, 784, ndig_treino, "train_dig0");
        MMQ_Alternados(W0, A0, H, 784, p, ndig_treino);
        deleteMatrix(A0);
        deleteMatrix(H);
  return 0;
}
