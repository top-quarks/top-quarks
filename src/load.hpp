// Minimal header library for progress bars
// Usage:
// for (int i = 0; i < 234; i++) {
//   load(234);
//   ...
// }

#include <stdio.h>

const int loading_len = 40;

void load_internal(long long a, long long b, const char*title = NULL, int keep_title = 0) {
  static long long prev = -1;
  long long pb = prev*b*2, x = a*2000+b;

  if (a == b) {
    if (title) {
      for (int i = 0; title[i]; i++)
	printf(" ");
    }
    for (int i = 0; i < loading_len+35+5; i++)
      printf(" ");
    printf("\r");
    fflush(stdout);
    if (title && keep_title) {
      printf("%s ", title);
      fflush(stdout);
    }
    return;
  }
  if (pb > x-b*2 && pb <= x) return;
  long long p = prev = x/(2*b);

  if (title) printf("%s ", title);
  printf("|");
  int dist = (a*2LL*loading_len+(b>>1))/b;
  for (int i = 0; i < dist>>1; i++)
    printf("=");
  if (dist&1)
    printf("-");
  for (int i = dist+1>>1; i < loading_len; i++)
    printf(" ");
  printf("| %lld / %lld (%lld.%lld%%)\r", a, b, p/10, p%10);
  fflush(stdout);
}

#define GET_FIRST_ARG(n,...) n
#define load(...) {static int cc = 0; load_internal(++cc, __VA_ARGS__); if (cc == GET_FIRST_ARG(__VA_ARGS__)) cc = 0;}
