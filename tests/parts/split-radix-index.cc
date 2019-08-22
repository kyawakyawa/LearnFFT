// http://www.kurims.kyoto-u.ac.jp/~ooura/fftman/ftmn1_24.html#sec1_2_4
#include <iostream>
using namespace std;
int main(void) {
  int n = 128;
  for (int m = n; m > 2; m >>= 1) {
    int mq = m >> 2;
    printf("m = %d\n", m);
    for (int i = 0; i < mq; i++) {
      for (int k = m; k <= n; k <<= 2) {
        for (int j = k - m + i; j < n; j += 2 * k) {
          printf("%d\n", j);
        }
      }
      printf("\n");
    }
    printf("\n\n");
  }

  printf("\n\n\n\n");

  for (int m = n; m > 2; m >>= 1) {
    int mq = m >> 2;
    printf("m = %d\n", m);
    for (int k = m; k <= n; k <<= 2) {
      for (int j = k - m; j < n; j += 2 * k) {
        for (int i = 0; i < mq; i++) {
          printf("%d\n", j + i);
        }
        printf("\n");
      }
    }
    printf("\n\n");
  }
  return 0;
}
