Manual de Instalação e Execução

Dependências
Para compilar e executar os programas, você precisa ter instalado:

- GCC (compilador C) com suporte a OpenMP
- Make (para usar o Makefile)

No Ubuntu instale com:

    sudo apt update
    sudo apt install build-essential

------------------------------------------------------------

Compilação
Basta rodar o comando:

    make

Isso irá compilar e gerar dois executáveis:

    jacobi_sequencial
    jacobi_paralelo

------------------------------------------------------------

Execução

Programa Sequencial
Uso:
    ./jacobi_sequencial <arquivo_entrada> <arquivo_saida>

Exemplo:
    ./jacobi_sequencial entrada.txt saida_seq.txt


Programa Paralelo
Uso:
    ./jacobi_paralelo <arquivo_entrada> <arquivo_saida> <num_threads> [<schedule_type>] [<chunk_size>]

- <num_threads> → número de threads (ex: 2, 4, 8)
- <schedule_type> → tipo de escalonamento (static, dynamic, guided)
- <chunk_size> → tamanho do bloco para divisão das iterações (opcional, default = 0)

Exemplos:

    # 4 threads, escalonamento estático
    ./jacobi_paralelo entrada.txt saida_par.txt 4 static

    # 8 threads, escalonamento dinâmico com blocos de 25
    ./jacobi_paralelo entrada.txt saida_par.txt 8 dynamic 25

------------------------------------------------------------

Entrada e Saída
- Arquivo de entrada deve conter:
  - Uma matriz A de dimensão 2000x2000
  - Um vetor b de tamanho 2000 na ultima linha do arquivo
  (Todos os valores em formato numérico, separados por espaço ou quebra de linha)

- Arquivo de saída conterá o vetor solução x, com valores impressos em com 4 casas decimais

