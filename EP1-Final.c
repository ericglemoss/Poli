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

  double sum;

    for ( int j = 0 ; j < p; j++ ) {
    sum = 0;
      for ( int i = 0; i < n; i++ ) {
        sum += pow(W[i][j], 2);
      }
      if(sum > eps){
        for ( int k = 0; k < n; k++ ) {
            W[k][j] /= sqrt(sum);
        }
      }
  }
}

void DefineMax (double** H, int p, int m) {
  for(int i = 0; i < p; i++){
    for (int j = 0; j < m; j++){
      if (H[i][j] < eps)
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

    double sum, c, s;
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
      W[i][j] = 0;
    }
  }
  return W;
}

double** MatrixMultiplication (double** W, double** H, int n, int m, int p) {
  double** A = CreateMatrix(n, m);
  A = zeros(n, m);
  int i, j, k;

     for (i = 0; i < n; i++ ){
        for (j = 0; j < m; j++ ){
            for (k = 0; k < p; k++ ){
            A[i][j] += W[i][k]*H[k][j];
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
    double Er = 1, E_atual, E_ant = 0;
    double** B   = CreateMatrix(n, m);
    double** A_aprox = CreateMatrix(n, m);
    double** A_t = CreateMatrix(m, n);
    double** H_t = CreateMatrix(m, p);
    double** W_t = CreateMatrix(p, n);

    CreatecopiaA(B, A, n, m);
    while((itmax < 100) & (Er > 0.0001)){
      E_atual = 0;
      H = zeros(p, m);
      W_t = zeros(p, n);
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
      for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++ )
            E_atual += pow(B[i][j] - A_aprox[i][j], 2);
      }
      Er = fabs(E_atual-E_ant)/E_atual;
      E_ant = E_atual;
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
  W = zeros(n, m);
  char d = 'a';
  int q = 0;
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

double** Treinar_matriz(int n, int m, int p, char *arquivo){
    double** A = CreateMatrix(n, m);
    double** W = CreateMatrix(n, p);

    A = Ler_arquivo(n, m, arquivo);
    W = MMQ_Alternados(A, n, p, m);
    deleteMatrix(A, n);

    return W;
}

int main () {

        int ndig_treino = 100;
        int n_test = 10000;
        int p = 5;
        int i;
        double** W[10];
        for(i = 0; i < 10; i++)
          W[i] = CreateMatrix(ndig_treino, p);

        double** A = CreateMatrix(784, n_test);
        double** B = CreateMatrix(784, n_test);

        W[0] = Treinar_matriz(784, ndig_treino, p, "train_dig0.txt");
        W[1] = Treinar_matriz(784, ndig_treino, p, "train_dig1.txt");
        W[2] = Treinar_matriz(784, ndig_treino, p, "train_dig2.txt");
        W[3] = Treinar_matriz(784, ndig_treino, p, "train_dig3.txt");
        W[4] = Treinar_matriz(784, ndig_treino, p, "train_dig4.txt");
        W[5] = Treinar_matriz(784, ndig_treino, p, "train_dig5.txt");
        W[6] = Treinar_matriz(784, ndig_treino, p, "train_dig6.txt");
        W[7] = Treinar_matriz(784, ndig_treino, p, "train_dig7.txt");
        W[8] = Treinar_matriz(784, ndig_treino, p, "train_dig8.txt");
        W[9] = Treinar_matriz(784, ndig_treino, p, "train_dig9.txt");

        int *gabarito = (int*) malloc(n_test * sizeof(int));
        int *resultados = (int*) malloc(n_test * sizeof(int));
        double *erros = (double*) malloc(n_test * sizeof(double));
        FILE*fpointer = fopen("test_index.txt", "r");
        for(i = 0; i < n_test; i++){
          erros[i] = INFINITY;
          fscanf(fpointer, "%d", &gabarito[i]);
        }

        fclose(fpointer);

        double** H[10];
        for(i = 0; i < 10; i++)
          H[i] = CreateMatrix(p, n_test);

        int k = 0, j = 0;
        A = Ler_arquivo(784, n_test, "test_images.txt");
        CreatecopiaA(B, A, 784, n_test);
        double **WxH = CreateMatrix(784, n_test);
        WxH = zeros(784, n_test);
        double Erro = 0, Erro2 = 0, somatoria = 0;

        for (i = 0; i < 10; i++){
          CreatecopiaA(A, B, 784, n_test);
          solveMultipleSystems(W[i], 784, p, n_test, H[i], A);

          WxH = MatrixMultiplication(W[i], H[i], 784, n_test, p);

          for (k = 0; k < n_test; k++){
            somatoria = 0;
            for (j = 0; j < 784; j++){
              Erro2 = A[j][k] - WxH[j][k];
              somatoria += (Erro2*Erro2);
            }
            Erro = sqrt(somatoria);
            if (erros[k] > Erro){
              erros[k] = Erro;
              resultados[k] = i;
            }
          }
        }
        int acertos = 0;
        int digito_correto[10];
        int digito_total[10];

        for(i = 0; i < 10; i++){
          digito_correto[i] = 0;
          digito_total[i] = 0;
        }

        for(i = 0; i < n_test; i++){
          digito_total[gabarito[i]]++;
          if(gabarito[i] == resultados[i]){
            acertos++;
            digito_correto[gabarito[i]]++;
          }
        }

        for (i = 0; i < 10; i++)
          printf("Digito %d: %6d/%5d corretos - %lf%%\n", i, digito_correto[i], digito_total[i], (double)digito_correto[i]/digito_total[i] * 100);

        printf("%6d/%5d corretos - %lf%%\n\n", acertos, n_test, (double) acertos/n_test * 100);
        for(i = 0; i < 10; i++){
          deleteMatrix(W[i], 784);
          deleteMatrix(H[i], p);
        }
        deleteMatrix(A, 784);
        deleteMatrix(B, 784);
        deleteMatrix(WxH, 784);


        free(gabarito);
        free(resultados);
        free(erros);

        printf("fim\n");
  return 0;

}