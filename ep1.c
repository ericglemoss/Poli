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
o tempo decorrido desde 01/01/1970. Essa última função foi utilizada para a geração de seeds para determinadas funções do programa. */
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


/* Função utilizada para desalocar uma matriz criada. Esse método tem como parâmetros a matriz (um ponteiro duplo)
e seu número de linhas (inteiro "n"). Para que haja o desalocamento, utilizou-se a função "free". O processo consiste
 no desalocamento de cada linha e, por fim, do ponteiro.*/
void DeleteMatrix(double** M, double n) {

    int i;

    for (i = 0; i < n; i++) {
      free(M[i]);
    }

    free(M);
}



/* Função utilizada para a impressão de matrizes (ponteiros duplos). Recebe como parâmetros uma matriz (ponteiro duplo "W"),
seu número de linhas (inteiro "n") e seu número de colunas (inteiro "m"). Seu funcionamento é simples: a matriz é percorrida
e cada elemento é impresso com uma precisão de 9 dígitos significativos, sendo 5 após a vírgula. */
void PrintMatrix(double** W, int n, int m) {

    for (int i = 0; i < n; i++ ) {
      for (int j = 0; j < m; j++ ) {
          printf (" %9.5f", W[i][j]);
      }
      printf ("\n");
    }
    printf ("\n");
}



/* Função utilizada para determinar a transposta de uma matriz. Recebe como parâmetros a matriz que terá sua transposta calculada
(ponteiro duplo "M"),seu número de linhas (inteiro "n") e número de colunas (inteiro "m") e a matriz recipiente da transposta ("N").
O método atua da seguinte forma: no escopo da função percorre-se a matriz M e fornece cada posição simétrica a um determinado
elemento ao recipiente da transposta. Ou seja, transformamos as linhas da matriz M em colunas na matriz N. */
 void MatrixTransposte(double** M, double** N , int n, int m) {

    for (int line = 0; line < n; line++) {
        for (int column = 0; column < m; column++) {
            N[column][line] = M[line][column];
        }
    }

}



/* Função destinada a criar uma matriz de elementos nulos. Recebe como parâmetros: uma matriz e suas linhas e colunas
(inteiros "n" e "m", respec.). No escopo da função, o ponteiro duplo é percorrido e seus elementos são sobreescritos com zero. */
void MakeElementsZero(double** W, int n, int m){
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++) {
      W[i][j] = 0.0;
    }
  }
}



/* Função destinada a retornar o produto de duas matrizes (W*H), assumindo que as dimensões matemáticas exigidas sejam válidas.
Recebe como parâmetros: duas matrizes (ponteiros duplos "W" e "H") e um recipiente ("WxH") e suas dimensões (inteiros "n", "m", "p").
A matriz W possui dimensão nxp, enquanto a matriz H dimensão pxm. O funcionamento do método é simples: respeitando-se as dimensões
fornecidas, as matrizes são percorridas e seus elementos multiplicados, somados eguardados em sua posição correspondente na matriz resultante. */
void MatrixMultiplication(double** W, double** H, double** WxH, int n, int m, int p) {

  int i, j, k;
  double sum;

  for (i = 0; i < n; i++){
      for (j = 0; j < m; j++){
        sum = 0;
        for (k = 0; k < p; k++){
            sum += W[i][k] * H[k][j];
            WxH[i][j] = sum;
        }
     }
  }
}



/* Função que toma uma matriz de referência e a copia em outra. Recebe como parâmetros: a matriz de referência (ponteiro duplo "A")
e suas dimensões (inteiros "n" e "m") e o recipiente para o qual será copiada (matriz "B"). O funcionamento desse método tem
como base percorrer a matriz A e copiar seus elementos para a matriz B. */
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

    if (fabs(Wik) > fabs(Wjk)) {
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




/* Função de aplicação da Rotação de Givens. Recebe como parâmetros: a matriz onde o procedimento matemático
 será aplicado (ponteiro duplo "W"),seu número de linhas (inteiro "n") e colunas (inteiro "m"),
 a coluna e linha onde ocorrerá as multiplicações (inteiros "i" e "j", respec.) e os valores dos cossenos
  e senos (reais "c" e "s", respectivamente).  O funcionamento dessa função segue àquele definido pelas diretrizes do exercício. */
void GivensRotation(double** W, int n, int m, int i, int j, double c, double s) {

    double temp;

    for (int k = 0; k < m; k++) {
        temp = c * W[i][k] - s * W[j][k];
        W[j][k] = s * W[i][k] + c * W[j][k];
        W[i][k] = temp;
    }
}



/* Função de decomposição QR, resolve um sistema Wx = b (onde "b" é um vetor coluna). Recebe como parâmetros: a matriz que
será decomposta (ponteiro duplo "W"), seu número de linhas e colunas (inteiros "n" e "m", respec.),
uma matriz onde serão armazenadas as respostas do sistema que será resolvido e um vetor de resultados (ponteiro duplo "b").
O funcionamento deste método baseia-se na aplicação sucessiva de Rotações de Givens na matriz de coeficientes (W) e de independente (b).
Esse procedimento resulta em uma matriz triangular superior em W, que então é utilizada para o cálculo do vetor resposta ("X"). */
void QRDecomposition(double** W, int n, int m, double** X, double** b) {

    double sum, c, s;
    int i, j, k;

    for (k = 0; k < m ; k++) {
        for (j = n - 1; j > k ; j--) {
            i = j - 1;
            if (fabs(W[j][k]) > eps) {
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
void SolveMultipleSystems(double** W, int n, int p, int m, double** H, double** A) {
    double sum, c, s;

    for (int k = 0; k < p; k++) {
        for (int j = n - 1; j > k; j--) {
            int i = j - 1;
            if (W[j][k] != 0) {
              CalculateSinNCos(W[i][k], W[j][k], &c, &s);
              GivensRotation(W, n, p, i, j, c, s);
              GivensRotation(A, n, m, i, j, c, s);
            }
        }
    }

    for (int k = p - 1; k >= 0; k--) {
        for (int j = 0; j < m; j++) {
            sum = 0;
            for (int i = k ; i < p ; i++) {
                sum += W[k][i] * H[i][j];
            }
            H[k][j] = (A[k][j] - sum)/W[k][k];
        }
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------*/


/*------------------------------------------------- FUNÇÕES DESTINADAS À SEGUNDA TAREFA DO EXERCÍCIO PROGRAMA------------------------------------------*/



/* Função utilizada para a normalização de colunas de uma matriz. Esse método tem como parâmetros uma matriz (ponteiro duplo "W"),
seu número de linhas (inteiro "n") e de colunas (inteiro "p"). Nesse processo, dada uma coluna, soma-se os quadrados
de seus elementos e, por fim, os divide individualmente pela raiz da soma. Foram tomados os devidos cuidados para casos
onde a soma das colunas é nula. */
void NormalizeColumns(double** W, int n, int p) {

  double sum;

    for (int j = 0 ; j < p; j++) {
        sum = 0;
        for (int i = 0; i < n; i++) {
            sum += pow(W[i][j], 2);
        }

      if (sum > eps) {
        for (int k = 0; k < n; k++) {
            W[k][j] /= sqrt(sum);
        }
      }
  }
}


/* Função utilizada  para definir o max{hij, 0}. Ela tem como objetivo eliminar termos negativos de uma matriz. Esse método recebe como
parâmetros uma matriz (ponteiro duplo "H"), seu número de linhas (inteiro "p") e seu número de colunas (inteiro "m"). Nesse exercício programa,
serão considerados zero todos os números menores que a constante definida eps (1e-5) */
void DefineMax(double** H, int p, int m) {

  for (int i = 0; i < p; i++){
    for (int j = 0; j < m; j++){
      if (H[i][j] < eps)
        H[i][j] = 0;
    }
  }
}


/* Função utilizada para iniciar aleatoriamente uma matriz, ou seja, cada elemento da mesma é um número real gerado aleatoriamente. Nesse processo,
foi utilizada a função "srand()" que, dada uma seed, gera um número. Recebe como parâmetros uma matriz (ponteiro duplo "W"), seu número de
linhas (inteiro "lines") e de colunas (inteiro "columns"). */
void RandomInit(double** W, int lines, int columns) {

  srand(time(NULL)); //Incialização da função random através de uma semente aleatória (time()).

  for (int i = 0; i < lines; i++) {
      for (int j = 0; j < columns; j++) {
          W[i][j] = (double)rand();
      }
  }
}



/* Função destinada à aplicação da fatoração por matrizes não negativas por meio do MMQ Alternado. Ou seja, dada uma matriz A, a fatoramos
em um produto de outra duas, W e H. O procedimento sequencial aqui instruído segue o definido pelos professores da disciplina
no pdf do exercício programa. Recebe como parâmetros: uma matriz(ponteiro duplo "A"), suas dimensões (inteiros "n" e "m")
e um inteiro "p", que será uma das dimensões da matriz fatorada. */
double** AlternatingLeastSquares(double** A, int n, int p, int m){

    double** W = CreateMatrix(n, p); //Criação da matriz W que será utilizada no procedimento matemático.
    double** H = CreateMatrix(p, m); //Criação da matriz H que é solução do sistema WH = A.
    double** A_copy   = CreateMatrix(n, m); //Criação do recipiente B onde será guardado uma cópia da matriz original A.
    double** A_aprox = CreateMatrix(n, m); //Criação da matriz que receberá o resultado aproximado do produto W*H.
    double** A_t = CreateMatrix(m, n); //Criação de uma matriz recipiente que receberá a transposta da matriz A.
    double** H_t = CreateMatrix(m, p); //Criação de uma matriz recipiente que receberá a transposta da matriz H.
    double** W_t = CreateMatrix(p, n); //Criação de uma matriz recipiente que receberá a transposta da matriz W.

    int iterations = 0; //Número de iterações do processo.
    double Er = 1; //Variável referente ao erro quadrático relativo entre duas iterações.
    double E_atual; //Variável referente ao erro quadrático da iteração atual.
    double E_ant = 0; //Variável referente ao erro quadrático da iteração anterior.

    RandomInit(W, n, p); //Inicialização (de forma aleatória) dos elementos da matriz W criada.

    CreateCopy(A_copy, A, n, m); //Criação de uma cópia da matriz A original para um recipiente (no caso, "A_copy").

    //O processo de fatoração tem como condição de parada um número máximo de iterações (100) e o valor do erro relativo.
    while ((iterations < 100) & (Er > 0.0001)){

      E_atual = 0; //A cada iteração, o erro atual deve ser zerado.

      MakeElementsZero(H, p, m); //Zera todos os elementos da matriz H.

      MakeElementsZero(W_t, p, n); //Zera todos os elementos da matriz recipiente W_t.

      CreateCopy(A, A_copy, n, m); //Repassa os valores da matriz A original, armazenados na matriz recipiente A_copy, novamente para a mesma.

      NormalizeColumns(W, n, p); //Normalização das colunas da matriz W.

      SolveMultipleSystems(W, n, p, m, H, A); //Resolução do sistema WH = A.

      DefineMax(H, p, m); //Eliminação de termos negativos na matriz H.

      MatrixTransposte(H, H_t, p, m); //Cálculo da matriz transposta H_t, que será utilizada para resolver o sistema H_t * W_t = A_t.

      MatrixTransposte(A_copy, A_t, n, m); //Cálculo da matriz A_t, transposta da matriz A original.

      SolveMultipleSystems(H_t, m, p, n, W_t, A_t); //Resolução do sistema W_t * H_t = A_t.

      MatrixTransposte(W_t, W, p, n); //Após a resolução do sistema, obtém-se novamente a matriz W por meio da transposta de W_t.

      DefineMax(W, n, p); //Eliminação de termos negativos da matriz W.

      MatrixMultiplication(W, H, A_aprox, n, m, p); //Cáculo do A_aprox através da multiplicação das matrizes W e H calculadas.

      for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++) {
            E_atual += pow(A_copy[i][j] - A_aprox[i][j], 2); //Cáculo do erro quadrático.
        }
      }

      Er = fabs(E_atual-E_ant)/E_atual; //Cálculo do erro relativo, que será condição de parada para o processo global.

      E_ant = E_atual; //Como será iniciada uma nova iteração em breve, o erro do processo anterior deve ser atualizado.

      iterations += 1; //Atualização do número de iterações.

    }

    DeleteMatrix(A_t, m); //Desalocamento dos ponteiros alocados para o processo, a fim de diminuir o custo de memória.
    DeleteMatrix(H_t, m);
    DeleteMatrix(W_t, p);
    DeleteMatrix(A_copy, n);
    DeleteMatrix(A_aprox, n);
    DeleteMatrix(H, p);

    return W;
}



void Tasks() {

    printf("                                        ______________________________________                                       \n");
    printf("--------------------------------------- |PRIMEIRA TAREFA DO EXERCICIO PROGRAMA| ----------------------------------- \n \n \n ");
    printf("Os resultados da primeira tarefa serao demonstrados abaixo. Por questao de espaco, serao demonstrados apenas os vetores solucao e as matrizes de tamanho conveniente. \n  ");



    double** W = CreateMatrix(64, 64); //WX = b; Matrizes do item A.
    double** b = CreateMatrix(64, 1);
    double** X = CreateMatrix(64, 1);

    for (int i = 0; i < 64; i++) {  //Inicialização da matriz W do item A.
        for (int j = 0; j < 64; j++) {
            if (i == j)
                W[i][j] = 2;

            if (abs(i - j) == 1)
                W[i][j] = 1;

            if (abs(i - j) > 1)
                W[i][j] = 0;
        }
        b[i][0] = 1; //Inicialização da matriz b do item A.
    }


    QRDecomposition(W, 64, 64, X, b); //Aplicação da decomposição QR nas matrizes W e b e cálculo do vetor solução X.

    printf("A matriz 'W' do item A apos sucessivas aplicacoes da Rotacao de Givens, devido ao seu tamanho extenso, esta mostrada no relatorio. \n \n");

    printf(" A matriz 'b' do item A apos sucessivas aplicacoes da Rotacao de Givens sera: \n \n");
    PrintMatrix(b, 64, 1);

    printf("O vetor solucao do item A sera: \n \n");
    PrintMatrix(X, 64, 1);

    DeleteMatrix(W, 64);
    DeleteMatrix(b, 64);
    DeleteMatrix(X, 64);

    printf("----------------------------------------------------------------------------------PROXIMO ITEM---------------------------------------------------------------------------- \n \n");


    double** K = CreateMatrix(20, 17); //KY = c; Matrizes do item B.
    double** Y = CreateMatrix(17, 1);
    double** c = CreateMatrix(20, 1);

    for (int i = 0; i < 20; i++) { // Inicialização da matriz W do item B.
        for (int j = 0; j < 17; j++) {
            if (abs(i - j) <= 4)
                K[i][j] = (double)1/(i + j + 1);

            if (abs(i - j) > 4)
                K[i][j] = 0;
        }
        c[i][0] = 1; //Inicialização da matriz b do item C.
    }


    QRDecomposition(K, 20, 17, Y, c); //Aplicação da decomposição QR nas matrizes W e b e cálculo do vetor solução.

    printf("A matriz 'W' do item B apos sucessivas aplicacoes da Rotacao de Givens sera: \n \n");
    PrintMatrix(K, 20, 17);

    printf("A matriz 'b' do item B apos sucessivas aplicacoes da Rotacao de Givens sera: \n \n");
    PrintMatrix(c, 20, 1);

    printf("O vetor solucao do item B sera: \n \n");
    PrintMatrix(Y, 17, 1);

    DeleteMatrix(K, 20);
    DeleteMatrix(Y, 17);
    DeleteMatrix(c, 20);

    printf("----------------------------------------------------------------------------------PROXIMO ITEM---------------------------------------------------------------------------- \n \n");



    double** A = CreateMatrix(64, 3); //ZH = A; Matrizes do item C.
    double** H = CreateMatrix(64, 3);
    double** Z = CreateMatrix(64, 64);

    for (int i = 0; i < 64; i++) {  //Inicialização da matriz W do item C.
        for (int j = 0; j < 64; j++) {
            if (i == j)
                Z[i][j] = 2;

            if (abs(i - j) == 1)
                Z[i][j] = 1;

            if (abs(i - j) > 1)
                Z[i][j] = 0;
        }
        A[i][0] = 1;    //Inicializações pertinentes à matriz A do item C.
        A[i][1] = i;
        A[i][2] = 2*i - 1;
    }


    SolveMultipleSystems(Z, 64, 64, 3, H, A); //Resolução de múltiplos sistemas.


    printf("A matriz 'W' do item C apos a utilizacao da funcao de resolucao de sistemas multiplos, devido ao seu tamanho, sera mostrada no relatorio. \n");

    printf("A matriz 'b' do item C apos a utilizacao da funcao de resolucao de sistemas multiplos sera: \n \n");
    PrintMatrix(A, 64, 3);

    printf("O vetor solucao do item C sera: \n \n");
    PrintMatrix(H, 64, 3);

    DeleteMatrix(A, 64);
    DeleteMatrix(H, 64);
    DeleteMatrix(Z, 64);

    printf("----------------------------------------------------------------------------------PROXIMO ITEM---------------------------------------------------------------------------- \n \n");


    double** B = CreateMatrix(20, 3); //Nv = B; Matrizes item D.
    double** N = CreateMatrix(20, 17);
    double** v = CreateMatrix(17, 3);

    for (int i = 0; i < 20; i++) {  //Inicialização da matriz W do item D.
        for (int j = 0; j < 17; j++) {
            if (abs(i - j) <= 4)
                N[i][j] = (double)1/(i + j + 1);

            else
                N[i][j] = 0;
        }
        B[i][0] = 1;    //Inicializações pertinentes a matriz A do item D.
        B[i][1] = i;
        B[i][2] = 2*i - 1;
    }

    SolveMultipleSystems(N, 20, 17, 3, v, B); //Resolução de múltiplos sistemas.

    printf("A matriz 'W' do item D apos a utilizacao da funcao de resolucao de sistemas multiplos sera: \n \n");
    PrintMatrix(N, 20, 17);

    printf("A matriz 'b' do item D apos a utilizacao da funcao de resolucao de sistemas multiplos sera: \n \n");
    PrintMatrix(B, 20, 3);

    printf("O vetor solucao do item D sera: \n \n");
    PrintMatrix(v, 17, 3);

    DeleteMatrix(B, 20);
    DeleteMatrix(N, 20);
    DeleteMatrix(v, 17);


    printf("                                        ______________________________________                                       \n");
    printf("--------------------------------------- |SEGUNDA TAREFA DO EXERCICIO PROGRAMA| ----------------------------------- \n \n \n ");
    printf(" Os resultados da segunda tarefa serao demonstrados abaixo. \n");

    double** A_s = CreateMatrix(3, 3); //Matrizes A e W da segunda tarefa.
    double** W_s = CreateMatrix(3, 2);

    A_s[0][0] = 0.3; A_s[0][1] = 0.6; A_s[0][2] = 0.0;
    A_s[1][0] = 0.5; A_s[1][1] = 0.0; A_s[1][2] = 1.0;
    A_s[2][0] = 0.4; A_s[2][1] = 0.8; A_s[2][2] = 0.0;

    printf("A matriz A inserida no processo de fatoracao por MMQ alternados foi: \n \n");
    PrintMatrix(A_s, 3, 3);


    printf("As matrizes W e H resultantes do processo de fatoracao serao: \n \n");
    W_s = AlternatingLeastSquares(A_s, 3, 2, 3);
    PrintMatrix(W_s, 3, 2);

}


/*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*------------------------------ FUNÇÕES DESTINADAS À TAREFA PRINCIPAL DO EXERCÍCIO PROGRAMA E DO APRENDIZADO DA MÁQUINA--------------------------------------------------------*/



/* Função destinada a ler um arquivo, pegar seu conteúdo e colocá-lo em uma matriz, retornando-a para main. Esse método foi criado para
facilitar a captura de dados para treinar a máquina dos arquivos .txt fornecidos pelos professores da disciplina. Recebe como parâmetros:
as dimensões da matriz a ser criada (inteiros "n" e "m" que são as linhas e colunas, respec.) e um ponteiro para o arquivo que será lido.
Cada linha lida do arquivo torna-se uma coluna da matriz gerada. Além disso, o elemento a ser guardado na mesma já é normalizado
(divido por 255, conforme é sugerido no arquivo do exercício programa). */
double** ReadFile(int n, int m, char *arquivo){

  FILE *fpointer = fopen(arquivo, "r");
  double** W = CreateMatrix(n, m);
  char d = 'a';
  int q = 0;

  MakeElementsZero(W, n, m);

  for (int k = 0; k < n; k++){
      for (int j = 0; d != '\r' && d != '\n' && d != EOF; j++){ //Leitura do arquivo, sempre prestantado atenção ao final de uma linha.
        fscanf(fpointer,"%d%c", &q, &d);
        if (j < m)
          W[k][j] = (double)q/255.0;
      }
      d = '\0';
  }
  fclose(fpointer);
  return W;
}



/* Função destinada ao treino de um dígito, ou seja, em seu escopo há o recebimento de um arquivo de dados nos quais serão
aplicadas as metodologias definidas no pdf do exercício programa para o treino da máquina. Recebe como parâmetros: as
dimensões das matrizes envolvidas no processo (inteiros "n", "m" e "p") e um ponteiro para o arquivo .txt utilizado.
O funcionamento desse método baseia-se solenemente na utilização de outros criados. */
double** Train_Matrix(int n, int m, int p, char *arquivo){

    double** A = CreateMatrix(n, m);
    double** W = CreateMatrix(n, p);

    A = ReadFile(n, m, arquivo);
    W = AlternatingLeastSquares(A, n, p, m);

    DeleteMatrix(A, n);

    return W;
}

int main (){

    //Váriaveis para o treino da máquina.
    int ndig_treino;
    int n_test;
    int p;

    printf("Entre com o valor de ndig_treino:  ");
    scanf("%d", &ndig_treino);

    printf("Entre com o valor de n_test: ");
    scanf("%d", &n_test);

    printf("Entre com o valor de p:  ");
    scanf("%d", &p);


    int i; //Iterador para auxiliar na criação de matrizes nas quais serão armazenadas os dados para o treino dos dígitos.
    double** W[10]; //Alocação do espaço da memória para as matrizes, conforme o número de dígitos treinados (10).

    for (i = 0; i < 10; i++)
        W[i] = CreateMatrix(ndig_treino, p); //Criação das 10 matrizes W que terão dimensões "ndig_treino x p".

    double** A = CreateMatrix(784, n_test);
    double** B = CreateMatrix(784, n_test);

    W[0] = Train_Matrix(784, ndig_treino, p, "train_dig0.txt"); //Treinamento para o dígito 0.
    W[1] = Train_Matrix(784, ndig_treino, p, "train_dig1.txt"); //Treinamento para o dígito 1.
    W[2] = Train_Matrix(784, ndig_treino, p, "train_dig2.txt"); //Treinamento para o dígito 2.
    W[3] = Train_Matrix(784, ndig_treino, p, "train_dig3.txt"); //Treinamento para o dígito 3.
    W[4] = Train_Matrix(784, ndig_treino, p, "train_dig4.txt"); //Treinamento para o dígito 4.
    W[5] = Train_Matrix(784, ndig_treino, p, "train_dig5.txt"); //Treinamento para o dígito 5.
    W[6] = Train_Matrix(784, ndig_treino, p, "train_dig6.txt"); //Treinamento para o dígito 6.
    W[7] = Train_Matrix(784, ndig_treino, p, "train_dig7.txt"); //Treinamento para o dígito 7.
    W[8] = Train_Matrix(784, ndig_treino, p, "train_dig8.txt"); //Treinamento para o dígito 8.
    W[9] = Train_Matrix(784, ndig_treino, p, "train_dig9.txt"); //Treinamento para o dígito 9.

    int* answers = (int*) malloc(n_test * sizeof(int)); //Alocação dinâmica para as respostas dos testes serem salvas.
    int* results = (int*) malloc(n_test * sizeof(int)); //Alocação dinâmica para os resultados obtidos pela máquina.
    double* erros = (double*) malloc(n_test * sizeof(double)); //Alocação dinâmica para os erros cometidos pela máquina.

    FILE* fpointer = fopen("test_index.txt", "r"); //Como o treino dos dígitos já foi realizado, abre-se o arquivo de testes para leitura dos dados presentes no mesmo.

    for (i = 0; i < n_test; i++){ //Leitura do arquivo aberto, conforme o número de testes determinado pelo usuário.
        erros[i] = INFINITY; //Os erros inicialmente são tomados como infinito, apenas para ocupação de posição.
        fscanf(fpointer, "%d", &answers[i]);
    }
        fclose(fpointer); //Fechamento do arquivo.

        double** H[10]; //Alocação de memória para as matrizes H que será utilizadas.
        for (i = 0; i < 10; i++)
            H[i] = CreateMatrix(p, n_test); //Alocação das matrizes H que serão utilizadas, todas de dimensões "p x n_test".

        int k = 0, j = 0;

        A = ReadFile(784, n_test, "test_images.txt"); //Leitura do arquivo onde estão dispostas as imagens para testes.
        CreateCopy(B, A, 784, n_test); //Os dados do arquivo lido são guardados em um repositório (no caso a matriz B).

        double** WxH = CreateMatrix(784, n_test); //Criação do repositório que irá guardar o produto WxH.
        MakeElementsZero(WxH, 784, n_test); //Todos os elementos do repositório são zerados.

        double Erro = 0, Erro2 = 0, sum = 0; //Váriaveis relativas ao erro quadrático do processo de classificação.

        for (i = 0; i < 10; i++){ //Loop destinado a analisar os erros quadrádicos. Como são 10 dígitos, o mesmo varia de 0 a 10.

            CreateCopy(A, B, 784, n_test); //Uma cópia da matriz de dados de teste é colocada alocada na matriz B.

            SolveMultipleSystems(W[i], 784, p, n_test, H[i], A); //Resolução do sistema WH = A.


            //Cálculo da matriz A aproximada, da qual serão utilizados os elementos para o cálculo do erro quadrático a seguir.
            MatrixMultiplication(W[i], H[i], WxH, 784, n_test, p);

            for (k = 0; k < n_test; k++) { //Cálculo do erro quadrático.
                sum = 0;
                for (j = 0; j < 784; j++) {
                    Erro2 = A[j][k] - WxH[j][k];
                    sum += pow(Erro2, 2);
                }

                Erro = sqrt(sum); //Cálculo da norma do erro quadrático.

            /* O erro referente a cada dígito é alocado em um vetor de erros. Como anteriormente os elementos desse vetor foram considerados INFINITY,
            essa condição sempre será acessada.*/
                if (erros[k] > Erro) {
                    erros[k] = Erro;
                    results[k] = i;
                }
               }
        }

    int scores = 0; //Número de acertos da máquina em relação aos dígitos analisados.
    int correct_digit[10]; //O número de acerto relativo a cada dígito é salvo em uma posição do vetor correct_digit. A primeira posição faz referência do dígito "0".
    int total_digits[10]; //O número total de análises para cada dígito é salvo em cada posição do vetor total_digits. A primeira posição faz referência ao dígito "0".

    for (i = 0; i < 10; i++){ //Inicialização dos vetores de acerto e dígitos totais analisados.
        correct_digit[i] = 0;
        total_digits[i] = 0;
    }

    for (i = 0; i < n_test; i++) { //Loop destinado para contabilizar a quantidade de acertos e dígitos analisados.
        total_digits[answers[i]]++;
        if (answers[i] == results[i]) {
            scores++;
            correct_digit[answers[i]]++;
          }
    }

    for (i = 0; i < 10; i++) //Impressão da quantidade de acertos por quantidade de dígitos analisados.
        printf("Digito %d: %6d/%5d corretos - %lf%%\n", i, correct_digit[i], total_digits[i], (double)correct_digit[i]/total_digits[i] * 100);

    printf("%6d/%5d corretos - %lf%%\n\n", scores, n_test, (double) scores/n_test * 100); //Impressão da porcentagem geral total de acertos (precisão) da máquina.

    for (i = 0; i < 10; i++){ //Desalocamento das matrizes criadas para o treino e classificação dos dígitos.
        DeleteMatrix(W[i], 784);
        DeleteMatrix(H[i], p);
    }

        DeleteMatrix(A, 784); //Desalocamento de matrizes utilizadas para a classificação dos dígitos.
        DeleteMatrix(B, 784);
        DeleteMatrix(WxH, 784);


        free(answers); //Deslocamento dos vetores de acertos, erros e respostas.
        free(results);
        free(erros);

  return 0;

}
