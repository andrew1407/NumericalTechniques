//steepest descent method
#include <iostream>
#include <math.h>

#define SIZE 4
#define ITERATIONS 20

const int X1 = 0, X2 = 4;
const float eps = pow(10, -6);
const float a = 0.3;

void F_x1(float* X) {
  const float x1 = X[0], x2 = X[1];
  const float res = sinf((x1 * 2) + 3) - 2 * x2 * cosf(x1 + 1.5) 
    - 5.8 * cosf(x1 + 1.5) + 2 * cosf(x2 - 2) + 2 * x1;
  X[0] = res;
}

void F_x2(float* X) {
  const float x1 = X[0], x2 = X[1];
  const float res = 2 * sinf(x1 + 1.5) + 2 * x2 + 2.9 
    - sinf(2 * x2 - 4) - 2 * x1 * sin(x2 - 2);
  X[1] = res;
}

void F(float* X) {
  F_x1(X);
  F_x2(X);
}

float maxValue(const float* x1, const float* x2) {
  const float res1 = abs(x1[0] - x2[0]);
  const float res2 = abs(x1[1] - x2[1]);
  return (res1 > res2) ? res1 : res2;
}

int main() {
  float X[2] = {X1, X2};
  float X_prev[2];

  bool breakpoint = true;
  float delta = 0;
  int iterationCounter = 0;

  printf("Interation 0:\n\tx1 = %f;\n\tx2 = %f;\n\tdelta = 0;\n", X[0], X[1]);

  while (breakpoint) {     
    X_prev[0] = X[0];
    X_prev[1] = X[1];
    F(X);
    delta = maxValue(X, X_prev);
    breakpoint = delta >= eps;
       
    printf("Iteration %d:\n", ++iterationCounter);
    printf("\tx1 = %.7f;\n", X[0]);
    printf("\tx2 = %.7f;\n", X[1]);
    printf("\tdelta = %.7f;\n", delta);

    if(iterationCounter > ITERATIONS) break;
  }

  return EXIT_SUCCESS;
}
