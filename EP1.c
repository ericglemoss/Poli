#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 10e-6

float DefineMax (float hij) {
  if (hij > 0)
    return hij;
  return 0;
}

void RandomInit (float ** W, int lines, int columns) {
  for ( int i = 0; i < lines; i++ ) {
      for ( int j = 0; j < columns; j++ ) {
          W[i][j] = rand();
      }
  }

}

float** CreateMatrix (int n, int m) {
  float ** Matrix = malloc(sizeof(float * ) * n);
  for (int i = 0; i < n; i++)
    Matrix[i] = (float * ) malloc(m * sizeof(float));

    for (int i = 0; i < n; i ++) {
        for (int j = 0; j < m; j++) {
            Matrix[i][j] = 0.0;
        }
    }
  return Matrix;
}

void DeleteMatrix(float ** M) {
  free(M[0]);
  free(M);
  printf("Matrix deleted! \n");
  printf("\n");
}

void PrintMatrix(float ** W, int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      printf(" %1.0f", W[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

float** MatrixTransposte(float ** M, int n, int m) {

  float** T = CreateMatrix(m, n);

  for (int line = 0; line < n; line++) {
    for (int column = 0; column < m; column++) {
      T[line][column] = M[column][line];
    }
  }
  return T;
}

void QRSinNCos(float Wjk, float Wik, float * c, float * s) {
  float tau;

  if (fabs(Wik) > fabs(Wjk)) {
    tau = -Wjk / Wik;
    * c = 1 / sqrt(1 + tau * tau);
    * s = * c * tau;
  } else {
    tau = -Wik / Wjk;
    * s = 1 / sqrt(1 + tau * tau);
    * c = * s * tau;
  }

}

void GivensRotation(float ** W, int n, int m, int i, int j, float c, float s) {
  float temp;

  for (int k = 0; k <= m - 1; k++) {
    temp = c * W[i][k] - s * W[j][k];
    W[j][k] = s * W[i][k] + c * W[j][k];
    W[i][k] = temp;
  }

}

void QRDecomposition(float** W, int n, int m, float** X, float** b) {

  float sum, c, s, tau;
  int i, j, k;

  for (k = 0; k < m ; k++) {
    for (j = n - 1; j > k; j--) {
      i = j - 1;
      if (fabs(W[j][k]) > EPS) {
            QRSinNCos(W[j][k], W[i][k], & c, & s);
            GivensRotation(W, n, m, i, j, c, s);
            GivensRotation(b, n, 1, i, j, c, s);
        }
      }
    }


  for (k = n - 1; k >= 0; k--) {
    sum = 0;
    for (j = k; j < m; j++) {
      sum += W[k][j] * X[j][0];
      X[k][0] = (double)(b[k][0] - sum) / W[k][k];
    }
  }

}

float** SolveMultipleSystems(float ** W, float ** A, int n, int p, int m) {

  int i, j, k;
  float sum, c, s;
  float **H = CreateMatrix(p, m);

  for (k = 0; k < p; k++) {
    for (j = n - 1; j < k; j--) {
      i = j - 1;
      if (fabs(W[j][k]) > EPS) {
        QRSinNCos(W[j][k], W[i][k], & c, & s);
        GivensRotation(W, n, p, i, j, c, s);
        GivensRotation(A, n, m, i, j, c, s);
      }

    }
  }

  for (k = p - 1; k >= 0; k--) {
    sum = 0;
    for (j = 0; j < m; j++) {
      sum += W[k][i] * H[i][j];
      H[k][j] = (double)(A[k][j] - sum) / W[k][k];
    }
  }

  return H;

}

/*void AlternantingLeastSquares(float ** A, int n, int m, int p) {

  float ** W = CreateMatrix(n, p);
  RandomInit(W);

  float error_norm;
  int iterations;
  float ** A_copy = A;

  for ()

} */

int main() {

    float** M = CreateMatrix(14, 14);
    float** A = CreateMatrix(14, 3);
    float** H = CreateMatrix(14, 3);

    for (int i = 0; i < 14; i ++) {
        for (int j = 0; j < 14; j++ ) {
            if (abs(i-j) == 1)
                M[i][j] = 1;

            if (i == j)
                M[i][j] = 2;

            if ( abs(i - j) > 1)
              M[i][j] = 0;
        }
    }

    for (int i = 0; i < 14; i++){
        A[i][0] = 1;
        A[i][1] = i;
        A[i][2] = 2*i + 1;
    }


        printf("\n");
        printf("\n");

        PrintMatrix(M, 14, 14);
        PrintMatrix(A, 14, 3);


        //SolveMultipleSystems(M, 14, 14, 3, H, A);

        PrintMatrix(M, 14, 14);


        printf("\n");

        PrintMatrix(A, 14, 3);
        printf("\n");

        PrintMatrix(H, 14, 3);

        printf("\n");

        DeleteMatrix(M);


        printf("\n");
        printf("\n");

    return 0;

  }
