#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char *argv[]) {

  char csv_file[256];
  int N = 0;
  int COVNO = 0;
  char cov_line[512];

  FILE *config = fopen(argv[1], "r");

  fscanf(config, "file=%s\nn=%d\ncovno=%d\ncovariates=%s", csv_file, &N, &COVNO,
         cov_line);

  // not char arrays are pointer by themselves no need address of &

  printf("%s\n", csv_file);

  FILE *file = fopen(csv_file, "r");

  char line[1024 * 2];

  fgets(line, sizeof(line), file);

  printf("%s\n", line);

  char *token;

  char *tokens[10] = {NULL};

  int i = 0;

  token = strtok(line, ",");

  printf("%s\n", token);

  while (token != NULL) {
    int val = atoi(token);
    printf("Token (string): %s, Value (int): %d\n", token, val);

    token = strtok(NULL, ",");

    tokens[i] = token;
    i++;
    //    if (strcmp(token, "fin") == 0)
    //    printf("%s\n", "got it");
  }

  printf("%s\n", tokens[3]);

  for (int j = 0; j < i; j++) {
    if (tokens[j] != NULL && strcmp(tokens[j], "fin") == 0) {
      printf("got it at index %d\n", j);
    }
  }

  return 0;
}
