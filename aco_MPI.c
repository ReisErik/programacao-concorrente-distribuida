#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <sys/time.h>
#include <windows.h>


#define MAX_IT 1000
#define MAX_VERTICES 200

typedef struct {
    int id;
    int *caminho;
    int *visitados;
    float custo;
} Formiga;

int carregarGrafo(const char *nomeArquivo, int grafo[MAX_VERTICES][MAX_VERTICES], int tamanho) {
    FILE *arquivo = fopen(nomeArquivo, "r");

    if (arquivo == NULL) {
        perror("Erro ao abrir arquivo");
        return 0;
    }

    for (int i = 0; i < tamanho; i++) {
        for (int j = 0; j < tamanho; j++) {
            if (fscanf(arquivo, "%d", &grafo[i][j]) != 1) {
                printf("Erro ao ler grafo.\n");
                fclose(arquivo);
                return 0;
            }
        }
    }

    fclose(arquivo);
    return 1;
}

void setFormigas(Formiga **formigas){
    for(int i = 0; i< MAX_VERTICES; i++){
        Formiga *f = (Formiga *)malloc(sizeof(Formiga));
        f->id = i;
        f->caminho = malloc(MAX_VERTICES * sizeof(int));
        f->visitados = calloc(MAX_VERTICES, sizeof(int));
        f->custo = 0;
        f->caminho[0] = i;
        f->visitados[i] = 1;

        formigas[i] = f;
    }
}

float probabilidade_caminho(int grafo[MAX_VERTICES][MAX_VERTICES], int origem, int destino, float feromonio[MAX_VERTICES][MAX_VERTICES], float alfa, float beta){
    float prob = (powf(feromonio[origem][destino], alfa) ) * (powf((1.0f/grafo[origem][destino]),beta));
    return prob;
}

double getTimestamp()
{
    struct timeval tv;

    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand((unsigned)time(NULL) ^ (rank * 7919));

    float alfa = 2.5; 
    float beta = 2.0;
    float decaimento = 0.6f;
    float deposito = 100.0f;

    int grafo[MAX_VERTICES][MAX_VERTICES];
    float feromonio[MAX_VERTICES][MAX_VERTICES];
    Formiga *formigas[MAX_VERTICES];

    // inicializa formigas
    setFormigas(formigas);
    
    FILE *arquivoSaida = NULL;
    if (rank == 0) {
        arquivoSaida = fopen("iteracoes_fast_mpi.csv", "w");
        if (arquivoSaida != NULL) {
            fprintf(arquivoSaida, "amostragem,processos,matriz,num_formigas,num_iteracoes,num_vertices_percorrido,tempo_total,inicio,fim\n");
        }
    }
    
    // tirar amostragens
    for(int amostragem = 0; amostragem < 5; amostragem++){
        int tamanhos[] = {150, 100, 50, 200};
        int qtd_tamanhos = 4;
        
        for (int g = 0; g < qtd_tamanhos; g++) {
            int NUM_VERTICES = tamanhos[g];

            char nomeArquivo[100];
            sprintf(nomeArquivo, "grafo_%d.txt", NUM_VERTICES);

            if (!carregarGrafo(nomeArquivo, grafo, tamanhos[g])) {
                MPI_Finalize();
                return 1;
            }

            Sleep(2000);
            double temp_inicio = getTimestamp();

            for (int i = 0; i < MAX_VERTICES; i++) {
                for (int j = 0; j < MAX_VERTICES; j++) {
                    feromonio[i][j] = 1.0f;
                }
            }

            // Divisão de formigas entre os processos MPI
            int NUM_FORMIGAS_POR_PROCESSO = NUM_VERTICES / size;
            int RESTO = NUM_VERTICES % size;

            int inicio = rank * NUM_FORMIGAS_POR_PROCESSO + (rank < RESTO ? rank : RESTO);
            int fim = inicio + NUM_FORMIGAS_POR_PROCESSO + (rank < RESTO ? 1 : 0);

            MPI_Barrier(MPI_COMM_WORLD);
            double tempo_inicio = MPI_Wtime();

            // inicia iterações
            for( int it = 0; it < MAX_IT ; it++ ){
                unsigned int seed = (unsigned int)time(NULL) ^ (rank * 104729 + it);
                
                for (int i = inicio; i < fim; i++){
                    Formiga *f = formigas[i];
                    // Loop caminho
                    for (int passo = 1; passo < NUM_VERTICES; passo++) {
                        float pesos[NUM_VERTICES];
                        float soma = 0;
                        int origem = f->caminho[passo - 1];
                        
                        // calcula peso (tij^alfa * 1/d^beta) e total para formula (somatorio dos pesos)
                        for (int j = 0; j < NUM_VERTICES; j++){
                            if(!f->visitados[j]){
                                pesos[j] = probabilidade_caminho(grafo, origem, j, feromonio, alfa, beta);
                            }
                            else{
                                pesos[j] = 0;
                            }
                            soma += pesos[j];
                        }
            
                        seed = (seed * 1103515245 + 12345) & 0x7fffffff;
                        float r = (float)seed / (float)0x7fffffff;
                        float acumulado = 0.0f;
                        int proximo_vertice = -1;
                        int ultimo_nao_visitado = -1;
            
                        // sorteia o proximo caminho
                        for (int j = 0; j < NUM_VERTICES; j++) {
                            if (!f->visitados[j]) {
                                ultimo_nao_visitado = j;
            
                                float probabilidade = pesos[j] / soma;
                                acumulado += probabilidade;
            
                                if (r <= acumulado) {
                                    proximo_vertice = j;
                                    break;
                                }
                            }
                        }
            
                        // caso o ultimo vertice seja sorteado, ou ocorra um erro, é atribuido
                        if (proximo_vertice == -1) {
                            proximo_vertice = ultimo_nao_visitado;
                        }
            
                        // salva informaçao
                        f->caminho[passo] = proximo_vertice;
                        f->visitados[proximo_vertice] = 1;
                        f->custo += grafo[origem][proximo_vertice];
                    }
            
                    // voltando a origem
                    int origem = f->caminho[NUM_VERTICES - 1];
                    int destino = f->caminho[0];
            
                    f->custo += grafo[origem][destino];
                }

                float pbest_local = 9999999.0f;
                float pbest_global;

                for (int i = inicio; i < fim; i++) {
                    if (formigas[i]->custo < pbest_local) {
                        pbest_local = formigas[i]->custo;
                    }
                }

                MPI_Allreduce(&pbest_local, &pbest_global, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

                // atualiza decaimento antes de somar (evitar sobrepor decaimento por formiga, apenas uma vez por iteraçao)
                for (int i = 0; i < NUM_VERTICES; i++) {
                    for (int j = 0; j < NUM_VERTICES; j++) {
                        feromonio[i][j] *= (1.0f - decaimento);
                    }
                }
                    
                // atualiza feromonio e reseta formiga
                float deposito_local[MAX_VERTICES][MAX_VERTICES] = {0};
                float deposito_global[MAX_VERTICES][MAX_VERTICES] = {0};

                for (int i = inicio; i < fim; i++){    
                    Formiga *f = formigas[i];    

                    // percorre caminho e calcula deposito local
                    for(int j = 0; j < (NUM_VERTICES - 1); j++){
                        int origem = f->caminho[j];
                        int destino = f->caminho[j+1];

                        deposito_local[origem][destino] += deposito / f->custo;
                    }
                        
                    // deposita a volta
                    int origem = f->caminho[NUM_VERTICES - 1];
                    int destino = f->caminho[0]; 

                    deposito_local[origem][destino] += deposito / f->custo;

                    // reset
                    for (int j = 0; j < NUM_VERTICES; j++) {
                        f->visitados[j] = 0;
                    }
                    
                    f->custo = 0.0f;
                    f->caminho[0] = f->id;
                    f->visitados[f->id] = 1;
                }

                MPI_Allreduce(deposito_local, deposito_global, MAX_VERTICES * MAX_VERTICES, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

                // Aplica feromonio global em todos os processos para manter matriz sincronizada
                for (int i = 0; i < NUM_VERTICES; i++) {
                    for (int j = 0; j < NUM_VERTICES; j++) {
                        feromonio[i][j] += deposito_global[i][j];
                    }
                }
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            double tempo_fim = MPI_Wtime();
            double temp_fim = getTimestamp();
            long long num_vertices_percorrido = (long long)tamanhos[g] * tamanhos[g] * MAX_IT;

            if (rank == 0) {
                if (arquivoSaida != NULL) {
                    fprintf(arquivoSaida, "%d,%d,%dx%d,%d,%d,%lld,%.6f,%.6f,%.6f\n",
                            amostragem, size, tamanhos[g], tamanhos[g], tamanhos[g],
                            MAX_IT, num_vertices_percorrido, tempo_fim - tempo_inicio, temp_inicio, temp_fim);
                    fflush(arquivoSaida);
                }
                printf("Matriz: %d Processos: %2d | Tempo: %.6f segundos | Inicio : %.6f | Amostragem %d\n",
                       tamanhos[g], size, tempo_fim - tempo_inicio, tempo_inicio, amostragem);
            }
        }
    }

    for (int i = 0; i < MAX_VERTICES; i++) {
        free(formigas[i]->caminho);
        free(formigas[i]->visitados);
        free(formigas[i]);
    }
    
    if (rank == 0 && arquivoSaida != NULL) {
        fclose(arquivoSaida);
    }

    MPI_Finalize();
    return 0;
}
}