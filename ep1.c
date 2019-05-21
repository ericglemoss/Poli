                                            /*PRIMEIRO EXERCÍCIO PROGRAMA DE CÁLCULO NUMÉRICO*/
                            /*----------------------------------------------------------------------------
                              |  NOME: Eric Gonçalves Lemos                                               |
                              |  NUSP: 9358916                                                            |
                              |  Curso: Engenharia Elétrica                                               |
                              |  Linguagem utilizada: C                                                   |
                              ----------------------------------------------------------------------------*/

/*As bibliotecas utilizadas na elaboração do código são naturais ao GCC.
As bibliotecas <stdio.h> e <stdlib.h> permitem a utilizações de comandos simples, tais como o de input e output "printf" e "scanf".
A biblioteca <math.h> foi utilizada por ter definida em seu cabeçalho funções matemáticas úteis, tais como "abs" (módulo de inteiro),
"fabs" (módulo de um número real) e "pow" (cálculo de uma potência). A biblioteca <time.h> foi utilizada para poder ter acesso a função time(), que retorna
um horário do computador, para a geração de uma seed utilizada em uma função que gera números aleatórios. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//A constante "eps" faz referência ao erro de 10e-5 definido na segunda tarefa.
#define eps 10e-5

/*--------------------------------------------------- FUNÇÕES AUXILIARES PARA OPERAÇÕES COM MATRIZES ----------------------------------------*/

/* Função utilizada para alocar dinamicamente uma matriz (ponteiro duplo) e inicializá-la com zeros. Esse método recebe como parâmetros o número de
 linhas (inteiro "n") e de colunas (inteiro "m") desejados. Nesse processo, utiliza-se a função de alocação dinâmica de memoria "malloc". */
double** CreateMatrix(int n, int m) {
    double** Matrix = malloc (sizeof (double*) * n );

    for (int i = 0; i < n; i++) {
        Matrix[i] = (double*)malloc(m * sizeof(double));
    }

    for (int k = 0; k < n; k++) {
        for (int j = 0; j < m; j++) {
            Matrix[k][j] = 0.0;
        }
    }

 return Matrix;

}


/* Função utilizada para desalocar uma matriz criada. Esse método tem como parâmetros a matriz (um ponteiro duplo) e seu número de linhas (inteiro "n").
Para que haja o desalocamento, utilizou-se a função "free". O processo consiste no desalocamento de cada linha e, por fim, do ponteiro.*/
void DeleteMatrix(double** M, double n) {

    int i;

    for(i = 0; i < n; i++) {
      free(M[i]);
    }
    free(M);
}



/* Função utilizada para a impressão de matrizes (ponteiros duplos). Recebe como parâmetros uma matriz (ponteiro duplo "W"), seu número de linhas
(inteiro "n") e seu número de colunas (inteiro "m"). Seu funcionamento é simples: a matriz é percorrida e cada elemento é impresso com uma precisão de
8 dígitos significativos, sendo 3 após a vírgula. */
void PrintMatrix (double** W, int n, int m) {

    for ( int i = 0; i < n; i++ ) {
      for ( int j = 0; j < m; j++ ) {
          printf (" %8.5f", W[i][j]);
      }
      printf ("\n");
    }
    printf ("\n");
}



/* Função utilizada para determinar a transposta de uma matriz. Recebe como parâmetros a matriz que terá sua transposta calculada (ponteiro duplo "M"),
seu número de linhas (inteiro "n") e número de colunas (inteiro "m") e matriz recipiente da transposta ("K"). O método atua da seguinte forma: no escopo da função, é criada uma nova matriz de
dimensões "m x n" (matriz transposta) que será retornada. Em conseguinte, a matriz "M" é percorrida e e seus elementos "Mji" são passados ao novo ponteiro
duplo criado. Por fim, retorna-se a transposta. */
 void MatrixTransposte (double** M, double** N , int n, int m) {

    for (int line = 0; line < n; line++) {
        for (int column = 0; column < m; column++) {
            N[column][line] = M[line][column];
        }
    }

}



/* Função destinada a criar uma matriz de elementos nulos. Recebe como parâmetros: linhas e colunas da matriz de interesse (inteiros "n" e "m", respec.).
Seu funcionamento está baseado na utilização de um método previamente criado ("CreateMatrix") e na inicialização em zero dos elementos. */
double** MakeElementsZero(int n, int m){
  double** W = CreateMatrix(n, m);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++) {
      W[i][j] = 0;
    }
  }
  return W;
}



/* Função destinada a retornar o produto de duas matrizes (W*H), assumindo que as dimensões matemáticas exigidas sejam válidas. Recebe como parâmetros: duas
matrizes (ponteiros duplos "W" e "H") e suas dimensões (inteiros "n", "m", "p"). A matriz W possui dimensão nxp, enquanto a matriz H dimensão pxm.
O funcionamento do método é simples: respeitando-se as dimensões fornecidas, as matrizes são percorridas e seus elementos multiplicados, somados e
guardados em sua posição correspondente na matriz resultante. */
void MatrixMultiplication (double** W, double** H, double** WxH, int n, int m, int p) {
  int i, j, k;
  for (i = 0; i < n; i++ ){
    for (j = 0; j < m; j++ ){
        for (k = 0; k < p; k++ ){
            WxH[i][j] += W[i][k]*H[k][j];
        }
     }
    }

}



/* Função que toma uma matriz de referência e a copia em outra. Recebe como parâmetros: a matriz de referência (ponteiro duplo "A") e suas dimensões
(inteiros "n" e "m") e o recipiente para o qual será copiada (matriz "B"). O funcionamento desse método tem como base percorrer a matriz A e copiar
seus elementos para a matriz B. */
void CreateCopy(double** B, double** A, int n, int m){
    for (int i = 0; i < n; i++){
      for (int j = 0; j < m; j++)
        B[i][j] = A[i][j];
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------*/


/*------------------------------------------ FUNÇÕES DESTINADAS À PRIMEIRA TAREFA DO EXERCÍCIO PROGRAMA---------------------------------------------*/


/* Função utilizada para o cálculo do seno ("s") e cosseno ("c") que serão utilizados na aplicação da Rotação de Givens. Recebe como parâmetros
dois elementos da matriz onde a rotação está sendo aplicada ("Wik" e "Wjk") e dois ponteiros nos quais serão armazenados os valores do seno ("*s")
e do cosseno ("*c"). O funcionamento deste método segue aquele determinado pelas diretrizes do exercício programa. A nomeclatura dos parâmetros
busca seguir de forma fiel o pdf da disciplina, a fim de facilitar a leitura do código. */
void CalculateSinNCos(double Wik, double Wjk, double *c, double *s) {
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




/* Função de aplicação da Rotação de Givens. Recebe como parâmetros: a matriz onde o procedimento matemático será aplicado (ponteiro duplo "W"), seu número de
linhas (inteiro "n") e colunas (inteiro "m"), a coluna e linha onde ocorrerá as multiplicações (inteiros "i" e "j", respec.) e os valores dos cossenos
e senos (reais "c" e "s", respectivamente). O funcionamento dessa função segue àquele definido pelas diretrizes do exercício. */
void GivensRotation ( double** W, int n, int m, int i, int j, double c, double s) {

    double temp;

    for ( int k = 0; k < m; k++ ) {
        temp = c * W[i][k] - s * W[j][k];
        W[j][k] = s * W[i][k] + c * W[j][k];
        W[i][k] = temp;
    }
}



/* Função de decomposição QR, resolve um sistema Wx = b (onde "b" é um vetor coluna). Recebe como parâmetros: a matriz que será decomposta (ponteiro duplo "W"), seu número de linhas e colunas (inteiros "n" e "m", respec.),
uma matriz onde serão armazenadas as respostas do sistema que será resolvido e um vetor de resultados (ponteiro duplo "b"). O funcionamento deste método baseia-se na
aplicação sucessiva de Rotações de Givens na matriz de coeficientes (W) e de independente (b). Esse procedimento resulta em uma matriz triangular superior em W, que
então é utilizada para o cálculo do vetor resposta ("X"). */
void QRDecomposition(double** W, int n, int m, double** X, double** b) {

    double sum, c, s;
    int i, j, k;

    for ( k = 0; k < m ; k++ ) {
        for ( j = n-1; j > k ; j-- ) {
            i = j - 1;
            if ( fabs(W[j][k]) > eps ) {
                CalculateSinNCos(W[i][k], W[j][k], &c, &s);
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




/* Função para a resolução de múltiplos sistemas, também baseada na aplicação sucessiva da Rotação de Givens. Ou seja, resolve um problema do tipo
WH = A. Como não é possível obter uma solução exata, esse método tem como objetivo minimizar o erro ||A-WH||. Recebe como parâmetros: uma matriz
de coeficientes (ponteiro duplo "W") e seu número de linhas e colunas (inteiros "n" e "p", respec.); um vetor solução (ponteiro duplo "H") e suas
dimensões (inteiros "p" e "m", respec.) e uma matriz independente (ponteiro duplo "A") de dimensões "n" por "m". O funcionamento dessa função tem
como base a aplicação sucessiva da Rotação de Givens nas matrizes W e A. Esse processo transforma a matriz W em triangular superior, tornando o problema
em outro de fácil resolução. */
void SolveMultipleSystems (double** W, int n, int p, int m, double** H, double** A) {
    double sum, c, s;

    for ( int k = 0; k < p; k++ ) {
        for ( int j = n - 1; j > k; j-- ) {
            int i = j - 1;
            if (W[j][k] != 0 ) {
              CalculateSinNCos(W[i][k], W[j][k], &c, &s);
              GivensRotation(W, n, p, i, j, c, s);
              GivensRotation(A, n, m, i, j, c, s);
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


/* Função destinada a resumir a primeira tarefa, mostrando a aplicação dos algoritmos desenvolvidos nos itens do pdf do exercício programa. Esse método
funciona exclusivamente na chamada de funções anteriormente definidas e explicadas. */
void FirstTask() {

    printf("                                        ______________________________________                                       \n");
    printf("--------------------------------------- |PRIMEIRA TAREFA DO EXERCICIO PROGRAMA| ----------------------------------- \n \n \n ");
    printf("Os resultados da primeira tarefa serao demonstrados abaixo. Por questao de espaco, serao demonstrados apenas os vetores solucao. \n \n ");



    double** W = CreateMatrix(64, 64); //WX = b; Matrizes do item A
    double** b = CreateMatrix(64, 1);
    double** X = CreateMatrix(64, 1);

    for (int i = 0; i < 64; i++) {
        for (int j = 0; j < 64; j++) {
            if (i == j)
                W[i][j] = 2;

            if (abs(i - j) == 1)
                W[i][j] = 1;

            if (abs(i - j) > 1)
                W[i][j] = 0;
        }
        b[i][0] = 1;
    }


    QRDecomposition(W, 64, 64, X, b);

    printf("A matriz W do item A apos sucessivas aplicacoes da Rotacao de Givens sera: \n \n");
    PrintMatrix(W, 64, 64);

    printf("A matriz b do item A apos sucessivas aplicacoes da Rotacao de Givens sera: \n \n" );
    PrintMatrix(b, 64, 1);

    printf ("O vetor solucao do item A sera: \n \n");
    PrintMatrix(X, 64, 1);

    DeleteMatrix(W, 64);
    DeleteMatrix(b, 64);
    DeleteMatrix(X, 64);

    printf("------------------------------------------------------------------------ \n \n");


    double** K = CreateMatrix(20, 17); //KY = c; Matrizes do item B
    double** Y = CreateMatrix(17, 1);
    double** c = CreateMatrix(20, 1);

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 17; j++) {
            if (abs(i - j) <= 4)
                K[i][j] = (double)1/(i + j + 1);

            if (abs(i - j) > 4)
                K[i][j] = 0;
        }
        c[i][0] = 1;
    }


    QRDecomposition(K, 20, 17, Y, c);

    printf("A matriz W do item B apos sucessivas aplicacoes da Rotacao de Givens sera: \n \n");
    PrintMatrix(K, 20, 17);

    printf("A matriz b do item B apos sucessivas aplicacoes da Rotacao de Givens sera: \n \n");
    PrintMatrix(c, 20, 1);

    printf("O vetor solucao do item B sera: \n \n");
    PrintMatrix(Y, 17, 1);

    DeleteMatrix(K, 20);
    DeleteMatrix(Y, 17);
    DeleteMatrix(c, 20);

    printf("------------------------------------------------------------------------ \n \n");



    double** A = CreateMatrix(64, 3); //ZH = A; Matrizes do item C
    double** H = CreateMatrix(64, 3);
    double** Z = CreateMatrix(64, 64);

    for (int i = 0; i < 64; i++) {
        for (int j = 0; j < 64; j++) {
            if (i == j)
                Z[i][j] = 2;

            if (abs(i - j) == 1)
                Z[i][j] = 1;

            if (abs(i - j) > 1)
                Z[i][j] = 0;

        }
        A[i][0] = 1;
        A[i][1] = i;
        A[i][2] = 2*i - 1;
    }


    SolveMultipleSystems(Z, 64, 64, 3, H, A);

    printf("A matriz W do item C apos a utilizacao da funcao de resolucao de sistemas multiplos sera: \n \n");
    PrintMatrix(Z, 64, 64);

    printf("A matriz b do item C apos a utilizacao da funcao de resolucao de sistemas multiplos sera: \n \n");
    PrintMatrix(A, 64, 3);

    printf("O vetor solucao do item C sera: \n \n");
    PrintMatrix(H, 64, 3);

    DeleteMatrix(A, 64);
    DeleteMatrix(H, 64);
    DeleteMatrix(Z, 64);

    printf("------------------------------------------------------------------------ \n \n");


    double** B = CreateMatrix(20, 3); //Nv = B; Matrizes item D;
    double** N = CreateMatrix(20, 17);
    double** v = CreateMatrix(17, 3);

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 17; j++) {
            if (abs(i - j) <= 4)
                N[i][j] = (double)1/(i + j + 1);

            else
                N[i][j] = 0;

        }
        B[i][0] = 1;
        B[i][1] = i;
        B[i][2] = 2*i - 1;
    }


    SolveMultipleSystems(N, 20, 17, 3, v, B);

    printf("A matriz W do item D apos a utilizacao da funcao de resolucao de sistemas multiplos sera: \n \n");
    PrintMatrix(N, 20, 17);

    printf("A matriz b do item D apos a utilizacao da funcao de resolucao de sistemas multiplos sera: \n \n");
    PrintMatrix(v, 17, 3);

    printf("O vetor solucao do item D sera: \n");
    PrintMatrix(v, 17, 3);

    DeleteMatrix(B, 20);
    DeleteMatrix(N, 20);
    DeleteMatrix(v, 17);


}


/*-----------------------------------------------------------------------------------------------------------------------------------------------------*/


/*------------------------------------------------- FUNÇÕES DESTINADAS À SEGUNDA TAREFA DO EXERCÍCIO PROGRAMA------------------------------------------*/



/* Função utilizada para a normalização de colunas de uma matriz. Esse método tem como parâmetros uma matriz (ponteiro duplo "W"), seu número de linhas (inteiro "n")
e de colunas (inteiro "p"). Nesse processo, dada uma coluna, soma-se os quadrados de seus elementos e, por fim, os divide individualmente pela raiz da soma.
Foram tomados os devidos cuidados para casos onde a soma das colunas é nula. */
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


/* Função utilizada  para definir o max{hij, 0}. Ela tem como objetivo eliminar termos negativos de uma matriz. Esse método recebe como
parâmetros uma matriz (ponteiro duplo "H"), seu número de linhas (inteiro "p") e seu número de colunas (inteiro "m"). Nesse exercício programa,
serão considerados zero todos os números muito pequenos (no contexto, menores que eps). */
void DefineMax (double** H, int p, int m) {
  for(int i = 0; i < p; i++){
    for (int j = 0; j < m; j++){
      if (H[i][j] < eps)
        H[i][j] = 0;
    }
  }
}


/* Função utilizada para iniciar aleatoriamente uma matriz, ou seja, cada elemento da mesma é um número real gerado aleatoriamente. Nesse processo,
foi utilizada a função "srand()" que, dada uma seed, gera um número. Recebe como parâmetros uma matriz (ponteiro duplo "W"), seu número de
linhas (inteiro "lines") e de colunas (inteiro "columns"). */
void RandomInit (double** W, int lines, int columns) {
  srand(time(NULL));
  for ( int i = 0; i < lines; i++ ) {
      for ( int j = 0; j < columns; j++ ) {
          W[i][j] = (double)rand();
      }
  }
}



/* Função destinada à aplicação da fatoração por matrizes não negativas por meio do MMQ Alternado. Ou seja, dada uma matriz A, a fatoramos em um produto de outra duas, W e H.
O procedimento sequencial aqui instruído segue o definido pelos professores da disciplina no pdf do exercício programa. Recebe como parâmetros: uma matriz
(ponteiro duplo "A"), suas dimensões (inteiros "n" e "m") e um inteiro "p", que será uma das dimensões da matriz fatorada. */
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

    CreateCopy(B, A, n, m);
    while((itmax < 100) & (Er > 0.0001)){
      E_atual = 0;
      H = MakeElementsZero(p, m);
      W_t = MakeElementsZero(p, n);
      CreateCopy(A, B, n, m);
      itmax += 1;
      NormalizeColumns(W, n, p);
      SolveMultipleSystems(W, n, p, m, H, A);
      DefineMax(H, p, m);
      MatrixTransposte(H, H_t p, m);
      MatrixTransposte(B, A_t, n, m);
      SolveMultipleSystems(H_t, m, p, n, W_t, A_t);
      MatrixTransposte(W_t, W, p, n);
      DefineMax(W, n, p);
     MatrixMultiplication(W, H, A_aprox, n, m, p);
      for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++ )
            E_atual += pow(B[i][j] - A_aprox[i][j], 2);
      }
      Er = fabs(E_atual-E_ant)/E_atual;
      E_ant = E_atual;
    }
    DeleteMatrix(A_t, m);
    DeleteMatrix(H_t, m);
    DeleteMatrix(W_t, p);
    DeleteMatrix(B, n);
    DeleteMatrix(A_aprox, n);
    DeleteMatrix(H, p);
    return W;
}

void SecondTask() {

    printf("                                        ______________________________________                                       \n");
    printf("--------------------------------------- |SEGUNDA TAREFA DO EXERCICIO PROGRAMA| ----------------------------------- \n \n \n ");
    printf("Os resultados da segunda tarefa serao demonstrados abaixo.  \n \n ");

    double** A = CreateMatrix(3, 3);
    double** W = CreateMatrix(3, 2);

    A[0][0] = 0.3; A[0][1] = 0.6; A[0][2] = 0.0;
    A[1][0] = 0.5; A[1][1] = 0.0; A[1][2] = 1.0;
    A[2][0] = 0.4; A[2][1] = 0.8; A[2][2] = 0.0;

    printf("A matriz A inserida no processo de fatoracao por MMQ alternados foi: \n \n");
    PrintMatrix(A, 3, 3);


    printf("As matriz H e W resultantes do processo de fatoracao serao: \n \n");
    W = MMQ_Alternados(A, 3, 2, 3);
    PrintMatrix(W, 3, 2);

}



/*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*------------------------------ FUNÇÕES DESTINADAS À TAREFA PRINCIPAL DO EXERCÍCIO PROGRAMA E DO APRENDIZADO DA MÁQUINA--------------------------------------------------------*/



/* Função destinada a ler um arquivo, pegar seu conteúdo e colocá-lo em uma matriz. Esse método foi criado para facilitar a captura de dados para treinar a máquina
dos arquivos .txt fornecidos pelos professores da disciplina. Recebe como parâmetros: as dimensões da matriz a ser criada (inteiros "n" e "m" que são as linhas e colunas, respec.) e
um ponteiro para o arquivo que será lido. */
double** Ler_arquivo(int n, int m, char *arquivo){
  FILE *fpointer = fopen(arquivo, "r");
  double** W = CreateMatrix(n, m);
  W = MakeElementsZero(n, m);
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



/* Função destinada ao treino de um dígito, ou seja, em seu escopo há o recebimento de um arquivo de dados nos quais serão aplicadas as metodologias definidas no pdf do exercício
programa para o treino da máquina. Recebe como parâmetros: as dimensões das matrizes envolvidas no processo (inteiros "n", "m" e "p") e um ponteiro para o arquivo .txt utilizado.
O funcionamento desse método baseia-se solenemente na utilização de outros criados. */
double** Treinar_matriz(int n, int m, int p, char *arquivo){
    double** A = CreateMatrix(n, m);
    double** W = CreateMatrix(n, p);

    A = Ler_arquivo(n, m, arquivo);
    W = MMQ_Alternados(A, n, p, m);
    DeleteMatrix(A, n);

    return W;
}

int main (){

         int ndig_treino = 4000;
        int n_test = 10000;
        int p = 15;
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
        CreateCopy(B, A, 784, n_test);
        double **WxH = CreateMatrix(784, n_test);
        WxH = MakeElementsZero(784, n_test);
        double Erro = 0, Erro2 = 0, somatoria = 0;

        for (i = 0; i < 10; i++){
          CreateCopy(A, B, 784, n_test);
          SolveMultipleSystems(W[i], 784, p, n_test, H[i], A);

          MatrixMultiplication(W[i], H[i], WxH, 784, n_test, p);

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
          DeleteMatrix(W[i], 784);
          DeleteMatrix(H[i], p);
        }
        DeleteMatrix(A, 784);
        DeleteMatrix(B, 784);
        DeleteMatrix(WxH, 784);


        free(gabarito);
        free(resultados);
        free(erros);

  return 0;

}
