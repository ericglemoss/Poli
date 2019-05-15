#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define eps 10e-5

void deleteMatrix (double** M, double n) {
    int i;

    for(i = 0; i < n; i++)
      free(M[i]);
    free(M);

    printf ("Matrix deleted! \n");
    printf ("\n");
}

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
  srand(12345);
  for ( int i = 0; i < lines; i++ ) {
      for ( int j = 0; j < columns; j++ ) {
          W[i][j] = (double)rand();
      }
  }
}

double** CreateMatrix(int n, int m) {

    double** Matrix = malloc (sizeof (double*) * n );

    for (int i = 0; i < n; i++ ){
        Matrix[i] = (double*)malloc(m * sizeof(double));
    }
 return Matrix;
}

/*Função destinada a impressão de matrizes, recebe como parâmetros um ponteiro duplo (vetor 2d)
  e suas dimensões n (linhas) e m (colunas) */
void printMatrix (double** W, int n, int m) {

    for ( int i = 0; i < n; i++ ) {
      for ( int j = 0; j < m; j++ ) {
          printf (" %8.3f", W[i][j]);
      }
      printf ("\n");
    }
    printf ("\n");
}

double** matrixTransposte (double** M, int n, int m) {

    double** T = CreateMatrix(m, n);

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

    #pragma omp parallel for
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

void solveMultipleSystems (double** W, int n, int p, int m, double** H, double** A) {
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

double** zeros(int n, int m){
  double** W = CreateMatrix(n, m);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++) {
      W[i][j] = 0.000;
    }

  }
  return W;
}

double** MatrixMultiplication (double** W, double** H, int n, int m, int p) {
  double** A = CreateMatrix(n, m);
  double soma;

     for ( int i = 0; i < n; i++ ){
        for ( int j = 0; j < m; j++ ){
            soma = 0;
            for ( int k = 0; k < p; k++ ){
                soma += W[i][k]*H[k][j];
                A[i][j] = soma;
            }
        }
     }
  return A;
}

double** MMQ_Alternados(double** A, int n, int p, int m){
    double** W = CreateMatrix(n, p);
    double **H = CreateMatrix(p, m);
    RandomInit(W, n, p);
    int itmax = 0;
    double sum;
    double* S;
    double Er = 1, E_atual = 1, E_ant = 0;
    double** B   = CreateMatrix(n, m);
    double** A_aprox = CreateMatrix(n, m);
    double** A_t = CreateMatrix(m, n);
    double** H_t = CreateMatrix(m, p);
    double** W_t = CreateMatrix(p, n);
    double soma_erro = 0;

    CreatecopiaA(B, A, n, m);
    while(itmax < 100 & Er > 0.00001){
      E_atual = 0;
      H = zeros(p, m);
      W_t = zeros(p, n);
      CreatecopiaA(A, B, n, m);
      itmax += 1;
      printf("%d\n", itmax);
      NormalizeColumns(W, n, p);
      solveMultipleSystems(W, n, p, m, H, A);
      DefineMax(H, p, m);
      H_t = matrixTransposte(H, p, m);
      A_t = matrixTransposte(B, n, m);
      solveMultipleSystems(H_t, m, p, n, W_t, A_t);
      W = matrixTransposte(W_t, p, n);
      DefineMax(W, n, p);
      if(itmax == 99){
        printMatrix(W, n,p);
        printMatrix(H, p, m);
      }
      A_aprox = MatrixMultiplication(W, H, n, m, p);
      printf("oi nunao\n");
      for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++ )
            E_atual += pow(B[i][j] - A_aprox[i][j], 2 );
      }
      Er = fabs(E_atual-E_ant)/E_atual;
      E_ant = E_atual;
      printf("%.5f\n", Er);
    }
    deleteMatrix(A_t, m);
    deleteMatrix(H_t, m);
    deleteMatrix(W_t, p);
    deleteMatrix(B, n);
    deleteMatrix(A_aprox, n);
    deleteMatrix(H, p);
    return W;
}

double** Ler_arquivo(int n, int m, char *arquivo){
  FILE *fpointer = fopen(arquivo, "r");
  double** W = CreateMatrix(n, m);
  char d;
  int q;
  for(int k = 0; k < n; k++){
      for(int j = 0; d != '\r' && d != '\n' && d != EOF; j++ ){
        fscanf(fpointer,"%d%c", &q, &d);
        if(j < m)
          W[k][j] = (double)q/255.0;
      }
      d = '\0';
  }
  fclose(fpointer);
  return W;
}

int main () {

        int ndig_treino = 4000;
        int n_test = 10000;
        int p = 15;
        int i;
        double** A0 = CreateMatrix(784, ndig_treino);
        double** W0 = CreateMatrix(784, p);
      //  double
        //FILE*fpointer = fopen("test_index.txt", "r");
      //  int *gabarito = (int*) malloc(n_test * sizeof(int));
      //  int *resultados = (int*) malloc(n_test * sizeof(int));
      //  double *erros = (double*) malloc(n_test * sizeof(double));
        A0 = Ler_arquivo(784, ndig_treino, "train_dig0.txt");
        W0 = MMQ_Alternados(A0, 784, p, ndig_treino);
        deleteMatrix(A0, 784);

        double** A1 = CreateMatrix(784, ndig_treino);
        double** W1 = CreateMatrix(784, p);
        A1 = Ler_arquivo(784, ndig_treino, "train_dig1.txt");
        W1 = MMQ_Alternados(A1, 784, p, ndig_treino);
        deleteMatrix(A1, 784);

        double** A2 = CreateMatrix(784, ndig_treino);
        double** W2 = CreateMatrix(784, p);
        A2 = Ler_arquivo(784, ndig_treino, "train_dig2.txt");
        W2 = MMQ_Alternados(A2, 784, p, ndig_treino);
        deleteMatrix(A2, 784);


        double** A3 = CreateMatrix(784, ndig_treino);
        double** W3 = CreateMatrix(784, p);
        A3 = Ler_arquivo(784, ndig_treino, "train_dig3.txt");
        W3 = MMQ_Alternados(A3, 784, p, ndig_treino);
        deleteMatrix(A3, 784);


        double** A4 = CreateMatrix(784, ndig_treino);
        double** W4 = CreateMatrix(784, p);
        A4 = Ler_arquivo(784, ndig_treino, "train_dig4.txt");
        W4 = MMQ_Alternados(A4, 784, p, ndig_treino);
        deleteMatrix(A4, 784);


        double** A5 = CreateMatrix(784, ndig_treino);
        double** W5 = CreateMatrix(784, p);
        A5 = Ler_arquivo(784, ndig_treino, "train_dig5.txt");
        W5 = MMQ_Alternados(A5, 784, p, ndig_treino);
        deleteMatrix(A5, 784);


        double** A6 = CreateMatrix(784, ndig_treino);
        double** W6 = CreateMatrix(784, p);
        A6 = Ler_arquivo(784, ndig_treino, "train_dig6.txt");
        W6 = MMQ_Alternados(A6, 784, p, ndig_treino);
        deleteMatrix(A6, 784);


        double** A7 = CreateMatrix(784, ndig_treino);
        double** W7 = CreateMatrix(784, p);
        A7 = Ler_arquivo(784, ndig_treino, "train_dig7.txt");
        W7 = MMQ_Alternados(A7, 784, p, ndig_treino);
        deleteMatrix(A7, 784);


        double** A8 = CreateMatrix(784, ndig_treino);
        double** W8 = CreateMatrix(784, p);
        A8 = Ler_arquivo(784, ndig_treino, "train_dig8.txt");
        W8 = MMQ_Alternados(A8, 784, p, ndig_treino);
        deleteMatrix(A8, 784);


        double** A9 = CreateMatrix(784, ndig_treino);
        double** W9 = CreateMatrix(784, p);
        A9 = Ler_arquivo(784, ndig_treino, "train_dig9.txt");
        W9 = MMQ_Alternados(A9, 784, p, ndig_treino);
        deleteMatrix(A9, 784);



        printf("fim\n");
  return 0;
}
