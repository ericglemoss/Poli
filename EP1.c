#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPS 10e-5

void CopyMatrix(double** B, double** A, int n, int m) {
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < m; j++ )
            B[i][j] = A[i][j];
}

void NormalizeVector (double** W, int n, int p) {

  double sum;

  for ( int j = 0 ; j < p; j++ ) {
    sum = 0;
    for ( int i = 0; i < n; i++ )
      sum += pow(W[i][j], 2);

    for ( int k = 0; k < n; k++ )
        W[k][j] /= sqrt(sum);

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

double** CreateMatrix(int n, int m) {

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

double** MatrixMultiplication (double** W, double** H, int n, int m, int p) {
    double** Result = CreateMatrix(n, m);

    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < m; j++ )
            for ( int k = 0; k < p; k++ )
                Result[i][j] += W[i][k] * H[k][j];

    return Result;
}

double** CalculateError (double** A_original, double** A_calculada, int n, int m) {
    double** Error = CreateMatrix(n, m);
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < m; j++ )
            Error[i][j] =   ( ( A_original[i][j] - A_calculada[i][j] )/A_original[i][j] ) * 100;

    return Error;
}

void DeleteMatrix (double** M) {
    free(M[0]);
    free(M);
    printf ("Matrix deleted! \n");
    printf ("\n");
}

void PrintMatrix (double** W, int n, int m) {

    for ( int i = 0; i < n; i++ ) {
        for ( int j = 0; j < m; j++ ) {

                    printf (" %8.5f", W[i][j]);
        }

        printf ("\n");
    }
    printf ("\n");
}

double** MatrixTransposte (double** M, int n, int m) {

    double** T = CreateMatrix(m, n);

    for (int line = 0; line < n; line++) {
        for (int column = 0; column < m; column++) {
            T[column][line] = M[line][column];
        }
    }

    return T;
}

void QRSin_N_Cos(double Wik, double Wjk, double *c, double *s) {
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

void GivensRotation ( double** W, int n, int m, int i, int j, double c, double s) {

    double temp;

    for ( int k = 0; k < m; k++ ) {
        temp = c * W[i][k] - s * W[j][k];
        W[j][k] = s * W[i][k] + c * W[j][k];
        W[i][k] = temp;
    }
}

void QRDecomposition(double** W, int n, int m, double** X, double** b) {

    double sum, c, s, tau;
    int i, j, k;

    for ( k = 0; k < m ; k++ ) {
        for ( j = n-1; j > k ; j-- ) {
            i = j - 1;
            if ( fabs(W[j][k]) > EPS ) {
                QRSin_N_Cos(W[i][k], W[j][k], &c, &s);
                GivensRotation(W, n, m, i, j, c, s);
                GivensRotation(b, n, 1, i, j, c, s);
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

double** SolveMultipleSystems (double** W, int n, int p, int m, double** H, double** A) {

    int i, j, k;
    double sum, c, s;

    for ( k = 0; k < p; k++ ) {
        for ( j = n - 1; j > k; j-- ) {
            i = j - 1;
            if ( fabs(W[j][k]) > EPS ) {
              QRSin_N_Cos(W[i][k], W[j][k], &c, &s);
              GivensRotation(W, n, p, i, j, c, s);
              GivensRotation(A, n, m, i, j, c, s);
            }
        }
    }

    for ( k = p - 1; k >= 0; k-- ) {
      for ( j = 0; j < m; j++ ) {
        sum = 0;
        for (i = k ; i < p; i++)
          sum += W[k][i] * H[i][j];
        H[k][j] = (A[k][j] - sum)/W[k][k];
      }
    }
}

void AlternatingLeastSquares(double** W, double** A, double** H, int n, int p, int m){
    RandomInit(W, n, p);
  //  PrintMatrix(W, 3, 2);
    int itmax = 0;
    double sum;
    double* S;
    double Er = 1, E_atual = 0, E_ant = 0;
    double** B;
    double** A_t = CreateMatrix(n, m);
    double** H_t = CreateMatrix(m, p);
    double** W_t = CreateMatrix(p, n);
    double soma_erro = 0;

    //CopyMatrix(B, A, m, n);

    PrintMatrix(B, m, n);
    while(itmax < 2 & Er > EPS){
    PrintMatrix(B, m, n);

      itmax += 1;
      printf("%d\n", itmax );
      NormalizeVector(W, n, p);
      SolveMultipleSystems(W, n, p, m, H, A);
      DefineMax(H, p, m);
      PrintMatrix(H, p, m);

      for (int i = 0; i < m; i++)
        for (int k = p; k < n; k++ )
          soma_erro += pow(A[k][i], 2);


      E_atual = sqrt(soma_erro);
      Er = (fabs(E_atual-E_ant)/E_atual);

      H_t = MatrixTransposte(H, p, m);
      A_t = MatrixTransposte(B, n, m);
     // PrintMatrix(H, p, m);
      //PrintMatrix(H_t, m, p);
      W_t = MatrixTransposte(W, n, p);
      //PrintMatrix(W_t, p, n);
      SolveMultipleSystems(H_t, m, p, n, W_t, A_t);
      //PrintMatrix(W_t, p, n);
      W = MatrixTransposte(W_t, p, n);
      DefineMax(W, n, p);
      //PrintMatrix(W, n, p);
      E_ant = E_atual;
      printf(" %g \n", Er );
      }

}


int main () {

        double** A = CreateMatrix(20, 3);
        double** W = CreateMatrix(20, 17);
        double** H = CreateMatrix(17, 3);


        for ( int i = 0; i < 20; i++ )
            for ( int j = 0; j < 17; j++ ) {
                if ( abs(i - j) <= 4 )
                    W[i][j] = (double) 1/(i + j + 1);

                else
                    W[i][j] = 0;
            }


        for ( int i = 0; i < 20; i++ ) {
            A[i][0] = 1;
            A[i][1] = i;
            A[i][2] = 2*i - 1;
        }

        PrintMatrix(A, 20, 3);

        printf("\n");

        PrintMatrix(W, 20, 17);

        SolveMultipleSystems(W, 20, 20, 17, H, A);

        PrintMatrix(H, 17, 3);

        printf("\n");

        DeleteMatrix(A);
        DeleteMatrix(W);
        DeleteMatrix(H);

        printf("\n");
        printf("\n");

  return 0;
}
