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
      printf(" %8.5f", W[i][j]);
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

void QRDecomposition(float ** W, int n, int m, float ** X, float ** b) {

  float sum, c, s, tau;
  int i, j, k;

  for (k = 0; k <= m - 1; k++) {
    for (j = n - 1; j > k; j--) {
      i = j - 1;
      if (fabs(W[j][k]) > eps) {
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
  float ** H = CreateMatrix(p, m);

  for (k = 0; k < p; k++) {
    for (j = n - 1; j < k; j--) {
      i = j - 1;
      if (fabs(W[j][k]) > EPS) {
        QRSinNCos(W[j][k], W[i][k], & c, & s);
        Givens_rotation(W, n, p, i, j, c, s);
        Givens_rotation(A, n, m, i, j, c, s);
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

void AlternantingLeastSquares(float ** A, int n, int m, int p) {

  float ** W = CreateMatrix(n, p);
  Random_Init(W);

  float error_norm;
  int iterations;
  float ** A_copy = A;

  for ()

}

int main() {

  float ** M = CreateMatrix(64, 64);
  float ** b = CreateMatrix(64, 1);
  float ** X = CreateMatrix(64, 1);

  DeleteMatrix(b);
  DeleteMatrix(X);
  DeleteMatrix(M);

  return 0;

}
