#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 2000

void ler_entrada(double **A, double *b, const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Erro ao abrir o arquivo de entrada");
        exit(EXIT_FAILURE);
    }
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

void escrever_saida(double *x, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        perror("Erro ao criar o arquivo de saída");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++)
    {
        fprintf(file, "%.4f ", x[i]);
    }
    fprintf(file, "\n");
    fclose(file);
}

void jacobi_sequencial(double **A, double *b, double *x)
{
    double *x_novo = (double *)malloc(N * sizeof(double));
    if (x_novo == NULL)
    {
        perror("Erro ao alocar memória para x_novo");
        exit(EXIT_FAILURE);
    }
    int iteracoes = 0;
    double erro;
    for (int i = 0; i < N; i++)
    {
        x[i] = 0.0;
    }
    do
    {
        erro = 0.0;
        for (int i = 0; i < N; i++)
        {
            if (fabs(A[i][i]) < 1e-10)
            {
                fprintf(stderr, "Erro: Elemento da diagonal A[%d][%d] é zero.\n", i, i);
                free(x_novo);
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
    } while (erro > 1e-5);
    printf("Método sequencial de Jacobi convergiu em %d iterações.\n", iteracoes);
    free(x_novo);
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "Uso: %s <arquivo_entrada> <arquivo_saida>\n", argv[0]);
        return EXIT_FAILURE;
    }
    double **A = (double **)malloc(N * sizeof(double *));
    if (A == NULL)
    {
        perror("Erro ao alocar memória para a matriz A");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < N; i++)
    {
        A[i] = (double *)malloc(N * sizeof(double));
        if (A[i] == NULL)
        {
            perror("Erro ao alocar memória para a linha da matriz A");
            for (int j = 0; j < i; j++)
            {
                free(A[j]);
            }
            free(A);
            return EXIT_FAILURE;
        }
    }
    double *b = (double *)malloc(N * sizeof(double));
    if (b == NULL)
    {
        perror("Erro ao alocar memória para o vetor b");
        for (int i = 0; i < N; i++)
            free(A[i]);
        free(A);
        return EXIT_FAILURE;
    }
    double *x = (double *)malloc(N * sizeof(double));
    if (x == NULL)
    {
        perror("Erro ao alocar memória para o vetor x");
        for (int i = 0; i < N; i++)
            free(A[i]);
        free(A);
        free(b);
        return EXIT_FAILURE;
    }
    ler_entrada(A, b, argv[1]);
    clock_t start_time = clock();
    jacobi_sequencial(A, b, x);
    clock_t end_time = clock();
    double cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("\nTempo de execução: %f segundos\n", cpu_time_used);
    escrever_saida(x, argv[2]);
    for (int i = 0; i < N; i++)
    {
        free(A[i]);
    }
    free(A);
    free(b);
    free(x);
    return 0;
}