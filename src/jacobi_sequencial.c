#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10 // O número de equações é 10

// Função para ler o sistema de equações do arquivo de entrada
void ler_entrada(double A[N][N], double b[N], const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Erro ao abrir o arquivo de entrada");
        exit(EXIT_FAILURE);
    }

    // Lê a matriz A (10 linhas, 10 valores por linha)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (fscanf(file, "%lf", &A[i][j]) != 1)
            {
                fprintf(stderr, "Erro ao ler a matriz A na linha %d, coluna %d\n", i + 1, j + 1);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Lê o vetor b (última linha do arquivo)
    for (int i = 0; i < N; i++)
    {
        if (fscanf(file, "%lf", &b[i]) != 1)
        {
            fprintf(stderr, "Erro ao ler o vetor b na coluna %d\n", i + 1);
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}

// Função para escrever a solução no arquivo de saída
void escrever_saida(double x[N], const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        perror("Erro ao criar o arquivo de saída");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; i++)
    {
        fprintf(file, "%.4f ", x[i]); // Saída com 4 casas decimais
    }
    fprintf(file, "\n");

    fclose(file);
}

// Implementação sequencial do método de Jacobi
void jacobi_sequencial(double A[N][N], double b[N], double x[N])
{
    double x_novo[N]; // Para armazenar a nova aproximação
    int iteracoes = 0;
    double erro;

    // Inicializa a solução com valores zero
    for (int i = 0; i < N; i++)
    {
        x[i] = 0.0;
    }

    do
    {
        erro = 0.0; // Reseta o erro para a nova iteração

        // Loop principal do método de Jacobi
        for (int i = 0; i < N; i++)
        {
            // Verifica se o elemento da diagonal é zero
            if (fabs(A[i][i]) < 1e-10)
            { // Usamos uma tolerância pequena para evitar erros de ponto flutuante
                fprintf(stderr, "Erro: Elemento da diagonal A[%d][%d] é zero, o que causa divisão por zero.\n", i, i);
                exit(EXIT_FAILURE);
            }
            double soma = 0.0;
            for (int j = 0; j < N; j++)
            {
                if (i != j)
                {
                    soma += A[i][j] * x[j];
                }
            }
            x_novo[i] = (b[i] - soma) / A[i][i];
        }

        // Calcula o erro e atualiza a solução
        for (int i = 0; i < N; i++)
        {
            double diff = fabs(x_novo[i] - x[i]);
            if (diff > erro)
            {
                erro = diff;
            }
            x[i] = x_novo[i];
        }

        iteracoes++;

    } while (erro > 1e-5); // Critério de parada: erro menor que 10e-5

    printf("Método sequencial de Jacobi convergiu em %d iterações.\n", iteracoes);
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "Uso: %s <arquivo_entrada> <arquivo_saida>\n", argv[0]);
        return EXIT_FAILURE;
    }

    double A[N][N];
    double b[N];
    double x[N];

    // Lê os dados
    ler_entrada(A, b, argv[1]);

    // Executa o método de Jacobi
    jacobi_sequencial(A, b, x);

    // Escreve a solução
    escrever_saida(x, argv[2]);

    return 0;
}