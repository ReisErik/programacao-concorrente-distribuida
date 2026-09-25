#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NUM_VERTICES 10
#define MAX_IT 100

typedef struct {
    int id;
    int *caminho;
    int *visitados;
    float custo;
} Formiga;

int carregarGrafo(const char *nomeArquivo, int grafo[NUM_VERTICES][NUM_VERTICES]) {
    FILE *arquivo = fopen(nomeArquivo, "r");

    if (arquivo == NULL) {
        perror("Erro ao abrir arquivo");
        return 0;
    }

    for (int i = 0; i < NUM_VERTICES; i++) {
        for (int j = 0; j < NUM_VERTICES; j++) {
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
    for(int i = 0; i< NUM_VERTICES; i++){
        Formiga *f = (Formiga *)malloc(sizeof(Formiga));
        f->id = i;
        f->caminho = malloc(NUM_VERTICES * sizeof(int));
        f->visitados = calloc(NUM_VERTICES, sizeof(int));
        f->custo = 0;
        f->caminho[0] = i;
        f->visitados[i] = 1;

        formigas[i] = f;
    }
}

float probabilidade_caminho(int grafo[NUM_VERTICES][NUM_VERTICES], int origem, int destino, float feromonio[NUM_VERTICES][NUM_VERTICES], float alfa, float beta){
    float prob = (powf(feromonio[origem][destino], alfa) ) * (powf((1.0f/grafo[origem][destino]),beta));
    return prob;
}

int main(){
    srand((unsigned)time(NULL));

    float alfa = 2.5; 
    float beta = 2.0;
    float decaimento = 0.6f;
    float deposito = 100.0f;

    int grafo[NUM_VERTICES][NUM_VERTICES];
    float feromonio[NUM_VERTICES][NUM_VERTICES];
    Formiga *formigas[NUM_VERTICES];

    if (!carregarGrafo("grafo.txt", grafo)) {
        return 1;
    }

    printf("Grafo carregado:\n\n");

    for (int i = 0; i < NUM_VERTICES; i++) {
        for (int j = 0; j < NUM_VERTICES; j++) {
            printf("%d ", grafo[i][j]);
        }
        printf("\n");
    }

    // iniciaciza formiga
    setFormigas(formigas);

    // inicializa feromonio
    for (int i = 0; i < NUM_VERTICES; i++) {
        for (int j = 0; j < NUM_VERTICES; j++) {
            feromonio[i][j] = 1.0f;
        }
    }

    // inicia iterações
    for( int it = 0; it < MAX_IT ; it++ ){
        printf("---------------------- ITERACAO %d ----------------------\n", it);

        // loop formigas
        for (int i = 0; i < NUM_VERTICES ; i++){
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

                    float r = (float) rand() / RAND_MAX;
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

        // calcula melhor custo da iteracao
        float pbest = formigas[0]->custo;

        for (int i = 1; i < NUM_VERTICES; i++) {
            if (formigas[i]->custo < pbest) {
                pbest = formigas[i]->custo;
            }
        }

        // atualiza decaimento antes de somar (evitar sobrepor decaimento por formiga, apenas uma vez por iteraçao)
        for (int i = 0; i < NUM_VERTICES; i++) {
            for (int j = 0; j < NUM_VERTICES; j++) {
                feromonio[i][j] *= (1.0f - decaimento);
            }
        }
        
        // atualiza feromonio e reseta formiga
        for (int i = 0; i < NUM_VERTICES ; i++){    
            Formiga *f = formigas[i];    

            //percorre caminho e atualiza feromonio
            for(int j = 0; j < (NUM_VERTICES - 1); j++){
                int origem = f->caminho[j];
                int destino = f->caminho[j+1];
                feromonio[origem][destino] += deposito/f->custo;
            }
            
            //deposita a volta
            int origem = f->caminho[NUM_VERTICES - 1];
            int destino = f->caminho[0]; 
            feromonio[origem][destino] += deposito / f->custo;

            // impressao do resultado
            printf("Formiga %d: ", f->id);
            for (int j = 0; j < NUM_VERTICES; j++) {
                printf("%d ", f->caminho[j]);
            }
            printf("-> %d | custo = %.2f\n", f->caminho[0], f->custo);

            
            //reset
            for (int j = 0; j < NUM_VERTICES; j++) {
                f->visitados[j] = 0;
            }
            f->custo = 0.0f;
            f->caminho[0] = f->id;
            f->visitados[f->id] = 1;
        }

        printf("PBest = %.2f\n", pbest);

    }

    for (int i = 0; i < NUM_VERTICES; i++) {
        free(formigas[i]->caminho);
        free(formigas[i]->visitados);
        free(formigas[i]);
    }

    return 0;
}