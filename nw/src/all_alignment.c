#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>
#include "rdtsc.h"
#include <omp.h>

#include "gene_bin.h"
#include "gene_bin.c"

#define MAX(x, y) ((x) > (y) ? (x) : (y))

//////////////// Print Similarity Matrix
/**
 * Prints the similarity matrix
 *
 * in : F : The similarity matrix in COL MAJOR. Its dimensions are `m*n`.
 * in : A : Horizontal char array. Its size is `m-1`
 * in : B : Vertical char array. Its size is `n-1`
 * in : m : Number of columns of `F`. Also the `A` length + 1
 * in : n : Number of lines of `F`. Also the `B` length + 1
 * in : k : Colors background of selected case in `F`, which position is `k` and `l`. Defaults it to 0 if not intended to use.
 * in : l : Colors background of selected case in `F`, which position is `k` and `l`. Defaults it to 0 if not intended to use.
 */
void print_sim_mat(int* F, char A [], char B [], int m, int n, int k, int l) {
    printf("%3s  %3s  ", "*", "-");
    for (int j = 0; j < m; j++)
        printf("%3c  ", A[j]);
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("%3c  ", i == 0 ? '-' : B[i - 1]);
        for (int j = 0; j < m; j++) {
            if (i == k && j == l)
                printf("\033[101m% 3d  \033[0m", F[i * m + j]);
            else
                printf("% 3d  ", F[i * m + j]);
        }
        printf("\n");
    }
}


//////////////// Calculate Similarity Matrix
/**
 * Calculates the similarity matrix of two sequences.
 * Sequences are stored as pairs of two bits, for each nucleotide.
 * One extra space for the similarity matrix is required, its returned dimensions are `(m+1)*(n+1)`
 *
 * in : bin_A : Horizontal binary array. Its size is `m*2`
 * in : bin_B : Vertical binary array. Its size is `n*2`
 * in : m : Number of columns of `F`. Also the `A` length + 1
 * in : n : Number of lines of `F`. Also the `B` length + 1
 * in  : match : Score for a match.
 * in  : mismatch : Score for a mismatch.
 * in  : gap : Score for a gap.
 * out : F : The similarity matrix in COL MAJOR. <!> Its dimensions are `(m+1)*(n+1)`.
 */
int* bin_calculate_scoring_matrix(long int* bin_A, long int* bin_B, int m, int n, int match, int mismatch, int gap) {
    m++;
    n++;

    int* F = NULL;
    F = (int*)calloc(m * n, sizeof(int));

    int diag = 0;
    int left = 0;
    int up = 0;

    // Setting F borders
    for (int i = 0; i < n; i++)
        F[i] = gap * i;
    for (int j = 0; j < m; j++)
        F[j * n] = gap * j;

    // Relative position of iterator in bin_X
    int i_bin_A, j_bin_B;

    // Two right bits
    unsigned short int aa, bb;

    //    G C A T G C G 
    //<=> G C G T A C G 
    //  0b10011011000110

    // Calculate scoring matrix
    for (int i = 1; i < n; i++) {
        i_bin_A = i / int_SIZE + (i % int_SIZE != 0) - 1;
        aa = (bin_A[i_bin_A] >> ((i - 1) * 2) % int_SIZE) & 0b11;
        for (int j = 1; j < m; j++) {
            j_bin_B = j / int_SIZE + (j % int_SIZE != 0) - 1;
            bb = (bin_B[j_bin_B] >> ((j - 1) * 2) % int_SIZE) & 0b11;
            diag = F[(i - 1) * m + j - 1] + (aa == bb ? match : mismatch);
            up = F[(i - 1) * m + j] + gap;
            left = F[i * m + j - 1] + gap;
            int var = MAX(MAX(diag, left), up);
            F[i * m + j] = var;
        }
    }
    return F;
}

//////////////// String Insert
/**
 * Insert a char at beginning of a char array.
 *
 * in/out : dest   : Destination string, to which we insert the `source` char
 * in     : source : Char to be inerted at beginning of `dest`.
 */
void string_insert(char** dest, char* source) {
    if (*dest == NULL) {
        *dest = calloc(sizeof(char), 2);
        memcpy(*dest, source, 1);
        return;
    }
    char* temp = calloc(sizeof(char), strlen(*dest) + 2);
    memcpy(temp, source, 1);
    strcat(temp, *dest);
    free(*dest);
    *dest = temp;
}

//////////////// Align
/**
 * Backtrace the similarity matrix, and prints the two sequences aligned.
 * Returns the max similarity score of the two sequences.
 *
 * in  : F : The similarity matrix in COL MAJOR. Its dimensions are given by (`A` length + 1)*(`B` length + 1).
 * in  : A : Horizontal char array.
 * in  : B : Vertical char array.
 * in  : match : Score for a match.
 * in  : mismatch : Score for a mismatch.
 * in  : gap : Score for a gap.
 * in  : print : Enables (1) or Disables (0) alignment output. Returns max similarity score either ways.
 * out : Max similarity score of the `A` and `B` sequences.
 */
int align(int* F, char A [], char B [], int match, int mismatch, int gap, int print) {
    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    if (print) {
        int i = m - 1;
        int j = n - 1;

        char* A_aligned = NULL;
        char* B_aligned = NULL;

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && F[i * m + j] == F[(i - 1) * m + j - 1] + (A[i - 1] == B[j - 1] ? match : mismatch)) {
                string_insert(&A_aligned, A + (i - 1));
                string_insert(&B_aligned, B + (j - 1));
                i--;
                j--;
            }
            else if (i > 0 && F[i * m + j] == F[(i - 1) * m + j] + gap) {
                string_insert(&A_aligned, A + (i - 1));
                string_insert(&B_aligned, "-");
                i--;
            }
            else {
                // else if (j > 0 && F[i * m + j] == F[i * m + j - 1] + gap) {
                string_insert(&A_aligned, "-");
                string_insert(&B_aligned, B + (j - 1));
                j--;
            }
        }
        printf("%s\n", A_aligned);
        printf("%s\n", B_aligned);
    }
    return F[m * n - 1];
}

int* calculate_scoring_matrix(char A [], char B [], int match, int mismatch, int gap) {

    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    int* F = NULL;
    F = (int*)calloc(m * n, sizeof(int));

    int diag = 0;
    int left = 0;
    int up = 0;

    // Setting F borders
    for (int i = 0; i < n; i++)
        F[i] = gap * i;
    for (int j = 0; j < m; j++)
        F[j * n] = gap * j;

    // Calculate scoring matrix
    for (int i = 1; i < n; i++)
        for (int j = 1; j < m; j++) {
            diag = F[(i - 1) * m + j - 1] + (A[j - 1] == B[i - 1] ? match : mismatch);
            up = F[(i - 1) * m + j] + gap;
            left = F[i * m + j - 1] + gap;
            int var = MAX(MAX(diag, left), up);
            F[i * m + j] = var;
        }
    return F;
}

int* diag_calculate_scoring_matrix(char A [], char B [], int match, int mismatch, int gap) {
    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    int* F = NULL;
    F = (int*)calloc(m * n, sizeof(int));

    int diag = 0;
    int left = 0;
    int up = 0;

    // Setting F borders
    for (int i = 0; i < m; i++)
        F[i] = gap * i;
    for (int j = 0; j < n; j++)
        F[j * m] = gap * j;

    int max_l, i, j;
    for (int k = 1; k < m + n - 2; k++) {
        if (k > (m - 1)) {
            max_l = n - (k - (m - 1)) - 1;
            if (max_l > m - 1) max_l = m - 1;
        }
        else max_l = k;

#pragma omp parallel for schedule(guided)
        for (int l = 1; l <= max_l; l++) {
            if (k > (m - 1)) {
                i = k - (m - 1) + l;
                j = (m - 1) - l + 1;
            }
            else {
                i = l;
                j = k - i + 1;
            }
            diag = F[(i - 1) * m + j - 1] + (A[j - 1] == B[i - 1] ? match : mismatch);
            up = F[(i - 1) * m + j] + gap;
            left = F[i * m + j - 1] + gap;
            int var = MAX(MAX(diag, left), up);
            F[i * m + j] = var;
        }
    }
    return F;
}

int countfiles() {
    int count = 0;
    struct dirent* entry;

    DIR* dir = opendir("./fastas/");

    while ((entry = readdir(dir)) != NULL) {
        if (strstr(entry->d_name, ".fasta"))
            count++;
    }
    closedir(dir);
    return count;
}

typedef struct files {
    int files_count;
    char** content;
    char** filenames;
}files_t;

int readfiles(files_t* files, int max_files) {
    DIR* dir;
    FILE* input;
    struct dirent* file;
    char* filepath = NULL;

    int files_count = countfiles();
    files->files_count = files_count;
    files->content = malloc(sizeof(char*) * files_count);
    files->filenames = malloc(sizeof(char*) * files_count);

    char folder_name[50] = "./fastas2/";

    // Open the directory which contain all the fastas files
    if ((dir = opendir(folder_name)) == NULL)
        return printf("Error: Can't open fastas folder\n"), -1;

    int i = 0;

    // Iterate if a file exists in this directory
    while ((file = readdir(dir)) != NULL) {

        // Skip directories (linux)
        if (strstr(file->d_name, ".fasta") && i < max_files) {

            printf("Reading file : %s\n", file->d_name);

            // Get filepath
            free(filepath);
            filepath = malloc(sizeof(char) * (strlen(folder_name) + strlen(file->d_name)));
            strcpy(filepath, folder_name);
            strcat(filepath, file->d_name);

            files->filenames[i] = malloc(200 * sizeof(char));
            strcpy(files->filenames[i], file->d_name);

            // Get file size (will be used for content allocation)
            struct stat st;
            stat(filepath, &st);
            long filesize = st.st_size;

            // Open fasta file
            if ((input = fopen(filepath, "r")) == NULL)
                return printf("Error: Can't open fastas file %s\n", filepath), -1;

            char* line = NULL;
            size_t len = 0;
            ssize_t read;

            // Skipping first line, as it's only fasta's metadata.
            getline(&line, &len, input);

            // Init content[i] with second line (first data line)
            files->content[i] = malloc((filesize - read - 1) * sizeof(char));
            getline(&line, &len, input);
            line[strcspn(line, "\n") - 1] = '\0';
            strcpy(files->content[i], line);

            // Concat each lines in content[i], while toggling newline.
            while ((read = getline(&line, &len, input)) != -1) {
                line[strcspn(line, "\n") - 1] = '\0';
                strcat(files->content[i], line);
            }
            fclose(input);
            free(line);
            i++;
        }
    }
    files->files_count = i;
    free(filepath);
    closedir(dir);
    return files_count;
}



// Basic main test, using bin's functions.
int main() {
    int match = 1;
    int mismatch = -1;
    int gap = -1;

    unsigned long long before, after;
    double elapsed;

    files_t files;
    readfiles(&files, 500);
    fprintf(stdout, "# \"type\" \"memory usage\" \"Execution Time\" \"cycles\"\n");
    // printf("%d\n", files.files_count);
    for (int i = 0; i < files.files_count - 1; i++) {
        char* A = files.content[i];
        char* B = files.content[i];
        long int* bin_A = NULL;
        long int* bin_B = NULL;
        bin_A = set_binary_array(A, strlen(A));
        bin_B = set_binary_array(B, strlen(B));

        int m = strlen(A) + 1;
        int n = strlen(B) + 1;

        before = rdtsc();
        int* F = calculate_scoring_matrix(A, B, match, mismatch, gap);
        after = rdtsc();
        int score = F[m * n - 1];
        elapsed = (double)(after - before);
        long int mem = m * n * sizeof(int) + (strlen(A) + strlen(B)) * sizeof(char);
        printf("char\tMax score : %d\n", score);
        // printf("-F size : %.1lfe6\n", (double)m * n / 1e6);
        // printf("-Memory allocation : %ld\n", mem);
        // printf("-Cycles count : %f\n", elapsed);
        // printf("-Time : %.3f sec\n", elapsed/3.3e9);
        fprintf(stdout, "%s %ld %.0f %.3f %s\n", files.filenames[i], mem, elapsed, elapsed / 3.3e9, "char");
        free(F);

        before = rdtsc();
        F = bin_calculate_scoring_matrix(bin_A, bin_B, strlen(A), strlen(B), match, mismatch, gap);
        after = rdtsc();
        score = F[m * n - 1];
        elapsed = (double)(after - before);
        mem = m * n * sizeof(int) + (strlen(A) / 16 + strlen(B) / 16) * sizeof(long int);
        printf("bin\tMax score : %d\n", score);
        // printf("*F size : %.1lfe6\n", (double)m * n / 1e6);
        // printf("*Memory allocation : %ld\n", mem);
        // printf("*Cycles count : %f\n", elapsed);
        // printf("*Time : %.3f sec\n", elapsed/3.3e9);
        fprintf(stdout, "%s %ld %.0f %.3f %s\n", files.filenames[i], mem, elapsed, elapsed / 3.3e9, "bin");
        free(F);

        before = rdtsc();
        F = diag_calculate_scoring_matrix(A, B, match, mismatch, gap);
        after = rdtsc();

        elapsed = (double)(after - before);
        score = F[m * n - 1];
        mem = m * n * sizeof(int) + (strlen(A) + strlen(B)) * sizeof(char);
        printf("diag\tMax score : %d\n", score);
        // printf("=Cycles count : %.0f\n", elapsed);
        // printf("=Time : %.3f sec\n", elapsed / 3.3e9);
        fprintf(stdout, "%s %ld %.0f %.3f %s\n", files.filenames[i], mem, elapsed, elapsed / 3.3e9, "diag");

        free(F);
        free(files.content[i]);
        free(bin_A);
        free(bin_B);
    }
    free(files.content[files.files_count - 1]);
    free(files.content);
}