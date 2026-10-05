#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <sys/time.h>

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

int main(){
    srand((unsigned)time(NULL));

    float alfa = 2.5; 
    float beta = 2.0;
    float decaimento = 0.6f;
    float deposito = 100.0f;

    int grafo[MAX_VERTICES][MAX_VERTICES];
    float feromonio[MAX_VERTICES][MAX_VERTICES];
    Formiga *formigas[MAX_VERTICES];

    // iniciaciza formiga
    
    FILE *arquivoSaida = fopen("iteracoes.csv", "w");
    fprintf(arquivoSaida, "amostragem,threads,matriz,num_formigas,num_iteracoes,num_vertices_percorrido,tempo_total,inicio,fim\n");
    
    // inicializa feromonio
    
    // tirar amostragens
    for(int amostragem = 0; amostragem < 5; amostragem++){
        int tamanhos[] = {50, 100, 150, 200};
        int qtd_tamanhos = 4;
        setFormigas(formigas);
        
        for (int g = 0; g < qtd_tamanhos; g++) {
            int NUM_VERTICES = tamanhos[g];

            char nomeArquivo[100];

            sprintf(
                nomeArquivo,
                "grafo_%d.txt",
                NUM_VERTICES
            );

            if (!carregarGrafo(nomeArquivo, grafo, tamanhos[g])) {
                return 1;
            }

        for (int num_threads = 1; num_threads < 13; num_threads++){
            omp_set_num_threads(num_threads);
            double temp_inicio = getTimestamp();

            for (int i = 0; i < MAX_VERTICES; i++) {
                for (int j = 0; j < MAX_VERTICES; j++) {
                    feromonio[i][j] = 1.0f;
                }
            }

            double inicio = omp_get_wtime();
            // inicia iterações
            for( int it = 0; it < MAX_IT ; it++ ){

                #pragma omp parallel
                {
                    unsigned int seed = (unsigned int)time(NULL) ^ omp_get_thread_num();
                    // loop formigas
                    #pragma omp for
                    for (int i = 0; i < NUM_VERTICES ; i++){
                        Formiga *f = formigas[i];
                        // Loop caminho
                        for (int passo = 1; passo < NUM_VERTICES; passo++) {
                            float pesos[NUM_VERTICES];
                            float soma = 0;
                            int origem = f->caminho[passo - 1];
                            
                            // calcula peso (tij^alfa * 1/d^beta) e total para formula (somatorio dos pesos)
                            #pragma omp simd reduction(+:soma)
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
            
                        //voltando a origem
                        int origem = f->caminho[NUM_VERTICES - 1];
                        int destino = f->caminho[0];
            
                        f->custo += grafo[origem][destino];
                    }

                }

                // calcula melhor custo da iteracao
                float pbest = formigas[0]->custo;

                #pragma omp parallel for reduction(min:pbest)
                for (int i = 1; i < NUM_VERTICES; i++) {
                    if (formigas[i]->custo < pbest) {
                        pbest = formigas[i]->custo;
                    }
                }

                // atualiza decaimento antes de somar (evitar sobrepor decaimento por formiga, apenas uma vez por iteraçao)
                #pragma omp parallel for collapse(2)
                for (int i = 0; i < NUM_VERTICES; i++) {
                    for (int j = 0; j < NUM_VERTICES; j++) {
                        feromonio[i][j] *= (1.0f - decaimento);
                    }
                }
                
                // atualiza feromonio e reseta formiga
                #pragma omp parallel for
                for (int i = 0; i < NUM_VERTICES ; i++){    
                    Formiga *f = formigas[i];    

                    //percorre caminho e atualiza feromonio
                    for(int j = 0; j < (NUM_VERTICES - 1); j++){
                        int origem = f->caminho[j];
                        int destino = f->caminho[j+1];

                        #pragma omp atomic
                        feromonio[origem][destino] += deposito/f->custo;
                    }
                    
                    //deposita a volta
                    int origem = f->caminho[NUM_VERTICES - 1];
                    int destino = f->caminho[0]; 

                    #pragma omp atomic
                    feromonio[origem][destino] += deposito / f->custo;

                    //reset
                    for (int j = 0; j < NUM_VERTICES; j++) {
                        f->visitados[j] = 0;
                    }
                    
                    f->custo = 0.0f;
                    f->caminho[0] = f->id;
                    f->visitados[f->id] = 1;
                }

            }
            

            double temp_fim = getTimestamp();
            double fim = omp_get_wtime();
            long long num_vertices_percorrido = (long long)tamanhos[g] * tamanhos[g] * MAX_IT;

            //"amostragem,threads,matriz,num_formigas,num_iteracoes,num_vertices_percorrido,tempo_total,inicio,fim
            fprintf(arquivoSaida,"%d,%d,%dx%d,%d,%d, %d,%.6f,%.6f,%.6f\n",amostragem,num_threads,tamanhos[g],tamanhos[g],tamanhos[g],MAX_IT,num_vertices_percorrido,fim-inicio,temp_inicio,temp_fim);
            printf("Matriz: %d Numero de Threads: %2d | Tempo: %.6f segundos | Inicio : %.6f | Amostragem %d\n", tamanhos[g] ,num_threads, fim - inicio, inicio, amostragem);
            }
        }

        for (int i = 0; i < MAX_VERTICES; i++) {
            free(formigas[i]->caminho);
            free(formigas[i]->visitados);
            free(formigas[i]);
        }
    }
    
    fclose(arquivoSaida);
    return 0;
}