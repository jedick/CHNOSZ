#include <ctype.h>
#include <string.h>
#include <stdlib.h>

/* [20180217] adapted from https://stackoverflow.com/questions/13213422/count-the-number-of-occurrences-of-each-letter-in-string */

void count_letters(char **chararray, int *counts) {

  char *str = chararray[0];

  int i;
  size_t len = strlen(str);

  for (i = 0; i < len; i++) {
      char c = str[i];
      if (!isalpha(c)) continue;
      counts[(int)(tolower(c) - 'a')]++;
  }

}
