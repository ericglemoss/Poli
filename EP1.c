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

      for ( int k = 0; k < n; k++ ) {
          W[k][j] /= sqrt(sum);
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

                      printf (" %8.3f", W[i][j]);
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

      int i, j, k;
      double sum, c, s;

      for ( k = 0; k < p; k++ ) {
          for ( j = n - 1; j > k; j-- ) {
              i = j - 1;
              if ( fabs(W[j][k]) > eps ) {
                qrSin_n_Cos(W[i][k], W[j][k], &c, &s);
                Givens_rotation(W, n, p, i, j, c, s);
                Givens_rotation(A, n, m, i, j, c, s);
              }
          }
      }
      for ( j = 0; j < m; j++  ) {
      for ( k = p - 1; k >= 0; k-- ) {
        sum = 0;
        for (i = k ; i < p ; i++){
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

  void MMQ_Alternados(double** W, double** A, double** H, int n, int p, int m){
      RandomInit(W, n, p);
    //  printMatrix(W, 3, 2);
      int itmax = 0;
      double sum;
      double* S;
      double Er = 1, E_atual = 0, E_ant = 0;
      double** B   = createMatrix(n, m);
      double** A_t = createMatrix(m, n);
      double** H_t = createMatrix(m, p);
      double** W_t = createMatrix(p, n);
      double soma_erro = 0;

      CreatecopiaA(B, A, n, m);
      while(itmax < 20 & Er > eps){
        CreatecopiaA(A, B, n, m);
        itmax += 1;
        printf("%d\n", itmax );
        printf("W\n");
        printMatrix(W,n,p);
        NormalizeColumns(W, n, p);
        printf("W\n");
        printMatrix(W,n,p);
      //  printMatrix(W, 3, 2);
        solveMultipleSystems(W, n, p, m, H, A);
        printf("H\n");
        printMatrix(H,p,m);
        DefineMax(H, p, m);
        printf("A\n");
        printMatrix(A,n,m);
        for (int i = 0; i < m; i++){
          for (int k = p; k < n; k++ )
            soma_erro += pow(A[k][i], 2);
        }
        E_atual = sqrt(soma_erro);
        Er = (fabs(E_atual-E_ant)/E_atual);
        printf("% .2f\n", Er);
        H_t = matrixTransposte(H, p, m);
        A_t = matrixTransposte(B, n, m);
        printMatrix(H_t, m, p);
        W_t = matrixTransposte(W, n, p);
        printMatrix(W_t, p, n);
        solveMultipleSystems(H_t, m, p, n, W_t, A_t);
        printMatrix(H_t, m, p);
        printMatrix(W_t, p, n);
        W = matrixTransposte(W_t, p, n);
        DefineMax(W, n, p);
        E_ant = E_atual;
        }

  }

  int main () {

          double** A = createMatrix(3, 3);
          double** W = createMatrix(3, 2);
          double** H = createMatrix(2, 3);

          A[0][0] = 0.3;
          A[1][0] = 0.5;
          A[2][0] = 0.4;
          A[0][1] = 0.6;
          A[1][1] = 0;
          A[2][1] = 0.8;
          A[0][2] = 0;
          A[1][2] = 1;
          A[2][2] = 0;

          printf("\n");

        //  printMatrix(A, 3, 3);

          MMQ_Alternados(W, A, H, 3, 2, 3);

          printf("\n");
  
          printMatrix(W, 3, 2);
          printMatrix(H, 2, 3);

          printf("\n");

          deleteMatrix(A);


          printf("\n");
          printf("\n");

    return 0;
  }
