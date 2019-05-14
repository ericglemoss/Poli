#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double** criar(int n, int m){
  double** matriz = malloc (sizeof (double*) * n );

  for (int i = 0; i < n; i++ )
      matriz[i] = (double*)malloc(m * sizeof(double));

      return matriz;

}

void Print_matriz(double ** e, int n, int m){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      printf("%0.4f ", e[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

double** matrixTransposte (double** M, int n, int m) {

    double** T = criar(m, n);

    for (int line = 0; line < n; line++) {
        for (int column = 0; column < m; column++) {
            T[line][column] = M[column][line];
        }
    }

    return T;
}

void Givens_c_s (double x, double y, double *c, double *s){
  double tau ;

  if (abs(x) > abs(y)){
    tau = -y/x;
    *c = 1/(sqrt(1+tau*tau));
    *s = tau**c;
  }
  else{
    tau = -x/y;
    *s = 1/(sqrt(1+tau*tau));
    *c = tau**s;
  }
}

void Rot_givens (double**w, int n, int m, int i, int j, double c, double s){
  double aux;
  int r;

  for(r = 0; r < m; r++){
    aux = c*w[i][r] - s*w[j][r];
    w[j][r] = s*w[i][r] + c*w[j][r];
    w[i][r] = aux;
  }
}

void QR_fatoracao(double** w, int n, int m, double** x, double** b){
    double somatoria, c, s;
    int i, j, k;

    for(int k = 0; k < m; k++){
      for(int j = n -1; j > k; j--){
        i = j-1;
        if (w[j][k] != 0){
          Givens_c_s(w[j][k], w[i][k], &c, &s);
          Rot_givens(w, n, m, i, j, c, s);
          Rot_givens(b, n, 1, i, j, c, s);
        }
      }
    }

    for( k = m - 1; k >= 0; k--){
      somatoria = 0;

      for( j = k+1; k < m; k++){
        somatoria = somatoria + w[k][j]*x[0][j];
      }
      x[0][k] = (b[k][0] - somatoria)/w[k][k];
  }
}

int main(){
  double** M = criar(10, 10);
    double** b = criar(10, 1);
    double** X = criar(10, 1);
    double** Q = criar(10, 10);
    double *c, *s;

    for (int i = 0; i < 10; i ++) {
        for (int j = 0; j < 10; j++ ) {
            if (i == j)
                M[i][j] = 1;

            else
                M[i][j] = 1;


        }
    }

        for (int i = 0; i < 10; i++)
            b[i][1] = 0;



        printf("\n");
        printf("\n");

        Givens_c_s(10, 2, c, s);

        printf("%d\n", *c);
        printf("%d\n", *s);

        Print_matriz(M, 10, 10);

        QR_fatoracao(M, 10, 10, X, b );

        Print_matriz(X, 10, 1);


        printf("\n");
        printf("\n");

        Print_matriz(M, 10, 10);

        Q = matrixTransposte(M, 10, 10);

        Print_matriz(M, 10, 10);

        printf("\n");

        Print_matriz(Q, 10, 10);

        printf("\n");
        printf("\n");

    return 0;
  }
