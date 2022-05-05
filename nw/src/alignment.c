#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX(x, y) ((x) > (y) ? (x) : (y))

void print_mat(int* F, char A [], char B [], int m, int n, int k, int l) {
    printf("%3s  %3s  ", "*", "-");
    for (int j = 0; j < m; j++)
        printf("%3c  ", A[j]);
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("%3c  ", i == 0 ? '-' : B[i - 1]);
        for (int j = 0; j < m; j++) {
            if (i == k && j == l)
                printf("\033[101m% 3d  \033[0m", F[i * n + j]);
            else
                printf("% 3d  ", F[i * n + j]);
        }
        printf("\n");
    }
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

void string_insert(char** dest, char* source) {
    if (*dest == NULL) {
        *dest = calloc(sizeof(char), 2);
        strncpy(*dest, source, 1);
        return;
    }
    char* temp = calloc(sizeof(char), strlen(*dest) + 2);
    strncpy(temp, source, 1);
    strcat(temp, *dest);
    free(*dest);
    *dest = temp;
}

void align(int* F, char A [], char B [], int match, int mismatch, int gap) {
    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    int i = m-1;
    int j = n-1;

    char* AlignmentA = NULL;
    char* AlignmentB = NULL;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && F[i * n + j] == F[(i - 1) * m + j - 1] + (A[i-1] == B[j-1] ? match : mismatch)) {
            string_insert(&AlignmentA, A + (i-1));
            string_insert(&AlignmentB, B + (j-1));
            i--;
            j--;
        }
        else if (i > 0 && F[i * m + j] == F[(i - 1) * m + j] + gap) {
            string_insert(&AlignmentA, A + (i-1));
            string_insert(&AlignmentB, "-");
            i--;
        }
        else {
        // else if (j > 0 && F[i * m + j] == F[i * m + j - 1] + gap) {
            string_insert(&AlignmentA, "-");
            string_insert(&AlignmentB, B + (j-1));
            j--;
        }
    }
    printf("%s\n", AlignmentA);
    printf("%s\n", AlignmentB);
}




int main() {
    char A [] = "GCATGCG";
    char B [] = "GATTACA";
    int m = strlen(A) + 1;
    int n = strlen(B) + 1;

    int match = 1;
    int mismatch = -1;
    int gap = -1;

    int* F = calculate_scoring_matrix(A, B, match, mismatch, gap);
    print_mat(F, A, B, m, n, n - 1, m - 1);

    align(F, A, B, match, mismatch, gap);

    free(F);
}