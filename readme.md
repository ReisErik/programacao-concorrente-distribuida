# Compilar e executar ACO
Programa Serial com mudanças de matriz, iteraçoes e formigas manual

gcc aco.c -o aco -lm
./aco.exe

# Compilar e executar o pipeline do ACO
Programa paralelizado com pipeline de execucação rodando 5x com valores de 1 a 12 threads e matrizes de 50, 100, 150 e 200 com 1000 iterações

gcc aco_paralelizado.c -o aco_p.exe -fopenmp -lm
./aco_p.exe

