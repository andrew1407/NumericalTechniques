//Newton simplified method
#include <iostream>
#include <math.h>

#define SIZE 4
#define ITERATIONS 20

//input values
const int X1 = 0, X2 = 4;
const float eps = pow(10, -6);

void F(const float* X, float* F_x) {
  const float x1 = X[0], x2 = X[1];
  F_x[0] = sinf(x1 + 1.5) - x2 + 2.9;
  F_x[1] = cosf(x2 - 2) + x1;
}

void Jackobian(const float* X, float* J_x) {
  const float x1 = X[0], x2 = X[1];
  J_x[0] = cosf(x1 + 1.5);
  J_x[1] = 1;
  J_x[2] = -1;
  J_x[3] = -sinf(x2 - 2);
}

void Jackobian_min1(float* J) {
  const float det = J[0] * J[3] - J[1] * J[2];
  const float el = J[1];
  J[1] = J[2];
  J[2] = el;
  for (int i = 0; i < SIZE; J[i++] /= det);
}

float maxValue(const float* x1, const float* x2) {
  const float res1 = abs(x1[0] - x2[0]);
  const float res2 = abs(x1[1] - x2[1]);
  return (res1 > res2) ? res1 : res2;
}

int main() {
  float X[2] = {X1, X2};
  float X_prev[2];
  float F_x[2];
  float J[4];
  bool breakpoint = true;
  float delta = 0;
  int iterationCounter = 0;

  printf("Interation 0:\n\tx1 = %f;\n\tx2 = %f;\n\tdelta = 0;\n", X[0], X[1]);

  while (breakpoint) {     
    X_prev[0] = X[0];
    X_prev[1] = X[1];
    F(X, F_x);
    Jackobian(X, J);
    Jackobian_min1(J);
    X[0] = X_prev[0] - (J[0] * F_x[0] + J[1] * F_x[1]);
    X[1] = X_prev[1] - (J[2] * F_x[0] + J[3] * F_x[1]);
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
