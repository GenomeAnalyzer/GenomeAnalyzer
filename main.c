#include "gene.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "mpi.h"
#include <sys/types.h>
#include <dirent.h>
#include <string.h>

int n = 1000;

unsigned int *alea_gen()
{
    unsigned int *gene = malloc(n * sizeof(int));

    srand(getpid());

    for (int i = 0; i < n; i++)
    {
        gene[i] = rand() % 2;
        // printf(" %d",gene[i]);
    }
    return gene;
}

int *test_gen()
{
    int *gene = calloc(sizeof(int *), MAX_GENES);

    // Example of a gene so the gene map register it
    gene[0] = 000011011100;

    return gene;
}

void print_test(gene_map_t *gene_map)
{
    printf("Début print de test\n");
    for (long long i = 0; i < gene_map->genes_counter; i++)
    {
        printf("Gene \x1b[31m%lld\x1b[0m: \x1b[32m[%lld:%lld]\x1b[0m\n", i, gene_map->gene_start[i], gene_map->gene_end[i]);
    }
}

/* Count files in a directory
 *
 */
int countfiles()
{
    int count = 0;
    struct dirent *entry;

    DIR *dir = opendir("fastas/");

    while ((entry = readdir(dir)) != NULL)
    {
        if (strcmp(entry->d_name, ".fasta") == 0)
            count++;
    }
    return count;
}

int readfiles()
{
    char **content = malloc(sizeof(char *) * countfiles());

    DIR *dir;

    FILE *input;

    // Open the directory which contain all the fastas files
    if ((dir = opendir("fastas")) == NULL)
        return printf("Error: Can't open fastas folder"), -1;

    struct dirent *file;

    // Iterate if a file exists in this directory
    while ((file == readdir(dir)))
    {
        // Skip parent directory ( linux)
        if ((!strcmp(file->d_name, ".")) && (!strcmp(file->d_name, "..")))
            continue;

        // Open fasta file
        if ((input = fopen(file->d_name, "r")) == NULL)
        {
            fclose(input);
            return printf("Error: Can't open fastas file %s", file->d_name), -1;
        }

        // first line , we need to skip it
        fscanf(input, "%s\n", content[0]);

        char *line;

        fgets( line, sizeof(line), input);

        printf("%s\n", line);

        fclose(input);
    }
    free(content);
}
int main(int argc, char **argv)
{
    int RANK_MASTER = 0;

    MPI_Init(&argc, &argv);

    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == RANK_MASTER)
    {
        readfiles();
    }

    MPI_Finalize();

    return 0;
}