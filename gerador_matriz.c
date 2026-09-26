#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){
    int n;
    int min = 1;
    int max = 100;
    char nome[50];

    srand(time(NULL));

    printf("Numero de vertices do grafo: ");
    scanf("%d", &n);

    sprintf(nome, "grafo_%d.txt", n);
    
    FILE *arquivo = fopen(nome, "w");

    if(arquivo == NULL){
        printf("Erro ao criar o arquivo");
        return 1;
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            int val;
            if(i == j){
                val = 0;
            }
            else{
                val = min + rand() % (max - min + 1);
            }
            fprintf(arquivo, "%d", val);
            if(j < n - 1) fprintf(arquivo, " ");
        }
        fprintf(arquivo, "\n");
    }
    fclose(arquivo);
    return 0;
}