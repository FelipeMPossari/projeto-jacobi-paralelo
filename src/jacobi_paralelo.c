#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define N 2000

void ler_entrada(double **A, double *b, const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        perror("Erro ao abrir o arquivo de entrada");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (fscanf(file, "%lf", &A[i][j]) != 1)
            {
                fprintf(stderr, "Erro ao ler a matriz A na linha %d, coluna %d\n", i + 1, j + 1);
                fclose(file);
                exit(EXIT_FAILURE);
            }

    for (int i = 0; i < N; i++)
        if (fscanf(file, "%lf", &b[i]) != 1)
        {
            fprintf(stderr, "Erro ao ler o vetor b na coluna %d\n", i + 1);
            fclose(file);
            exit(EXIT_FAILURE);
        }

    fclose(file);
}

void escrever_saida(double *x, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (!file)
    {
        perror("Erro ao criar o arquivo de saída");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++)
        fprintf(file, "%.4f ", x[i]);
    fprintf(file, "\n");
    fclose(file);
}

void jacobi_paralelo(double **A, double *b, double *x, int threads)
{
    double *x_novo = (double *)malloc(N * sizeof(double));
    if (!x_novo)
    {
        perror("Erro ao alocar memória para x_novo");
        exit(EXIT_FAILURE);
    }
    int iteracoes = 0;
    double erro;
    for (int i = 0; i < N; i++)
        x[i] = 0.0;

    omp_set_num_threads(threads);
    do
    {
        erro = 0.0;
#pragma omp parallel for schedule(runtime)
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
                if (i != j)
                    soma += A[i][j] * x[j];

            x_novo[i] = (b[i] - soma) / A[i][i];
        }

#pragma omp parallel for reduction(max : erro)
        for (int i = 0; i < N; i++)
        {
            double diff = fabs(x_novo[i] - x[i]);
            if (diff > erro)
                erro = diff;
            x[i] = x_novo[i];
        }

        iteracoes++;
    } while (erro > 1e-5);

    printf("Método paralelo de Jacobi convergiu em %d iterações.\n", iteracoes);
    free(x_novo);
}

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        fprintf(stderr, "Uso: %s <arq_entrada> <arq_saida> <num_threads> [<schedule_type>] [<chunk_size>]\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *valid_schedules[] = {"static", "dynamic", "guided"};
    const char *schedule_type = (argc > 4) ? argv[4] : "static";
    int is_valid = 0;
    for (int i = 0; i < 3; i++)
        if (strcmp(schedule_type, valid_schedules[i]) == 0)
            is_valid = 1;

    if (!is_valid)
    {
        fprintf(stderr, "Erro: Tipo de agendamento '%s' não é válido. Use: static, dynamic ou guided.\n", schedule_type);
        return EXIT_FAILURE;
    }

    double **A = (double **)malloc(N * sizeof(double *));
    double *b = (double *)malloc(N * sizeof(double));
    double *x = (double *)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++)
        A[i] = (double *)malloc(N * sizeof(double));

    int num_threads = atoi(argv[3]);
    int chunk_size = (argc > 5) ? atoi(argv[5]) : 0;

    char env_schedule[50];
    sprintf(env_schedule, "%s,%d", schedule_type, chunk_size);
    setenv("OMP_SCHEDULE", env_schedule, 1);

    ler_entrada(A, b, argv[1]);

    double tempos[3];
    double soma = 0.0;

    for (int k = 0; k < 3; k++)
    {
        double start_time = omp_get_wtime();
        jacobi_paralelo(A, b, x, num_threads);
        double end_time = omp_get_wtime();

        tempos[k] = end_time - start_time;
        soma += tempos[k];
        printf("Execução %d: tempo = %.6f s\n", k + 1, tempos[k]);
    }

    printf("Tempo médio: %.6f s\n", soma / 3.0);

    escrever_saida(x, argv[2]);

    for (int i = 0; i < N; i++)
        free(A[i]);
    free(A);
    free(b);
    free(x);

    return 0;
}
