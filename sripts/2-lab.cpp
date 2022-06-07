//Seidel iteration method
#include <stdio.h>
#include <iostream>

#define SIZE 4

// 4.855 * x1 + 1.239 * x2 + 0.272 * x3 + 0.258 * x4 = 1.192
// 1.491 * x1 + 4.954 * x2 + 0.124 * x3 + 0.236 * x4 = 0.256
// 0.456 * x1 + 0.285 * x2 + 4.354 * x3 + 0.254 * x4 = 0.852
// 0.412 * x1 + 0.335 * x2 + 0.158 * x3 + 2.874 * x4 = 0.862

//Main coefficients

float M[SIZE][SIZE + 1] = { /*{x1, x2, x3. x4, b}*/ 
  { 4.855, 1.239, 0.272, 0.258, 1.192 },
  { 1.491, 4.954, 0.124, 0.236, 0.256 },
  { 0.456, 0.285, 4.354, 0.254, 0.852 },
  { 0.412, 0.335, 0.158, 2.874, 0.862 }
};


//Functions

float maxNum(const float *x) {
  float m = x[0];
  for (int i = 1; i < SIZE; i++) m = (m > x[i]) ? m : x[i];
  return m;
}

float maxDelta(const float *x1, const float *x2) {
  float x[SIZE];
  float _x1, _x2;
  for (int i = 0; i < SIZE; i++) {
    _x1 = (x1[i] > 0) ? x1[i] : -x1[i];
    _x2 = (x2[i] > 0) ? x2[i] : -x2[i];
    float el = _x1 - _x2;
    x[i] = (el > 0) ? el : -el;
    printf("\tDelta x%d: %f\n", i+1, x[i]);
  }
  return maxNum(x);
}

//flag for main cycle
bool imparity(const float *x1, const float *x2, const float div) {
  const float delta = maxDelta(x1, x2);
  printf("Current delta (max): %f;\nPivotal misstep: %f\n", delta, div);
  return delta <= div;
}

void showX(char *str, const float *X) {
  printf(str);
  for (int i = 0; i < SIZE; i++) printf("x%d = %f;\t", i + 1, X[i]);
  printf("\n");
}

//[MAIN]
int main() {
  float X_prev[SIZE];
  float X[SIZE];
  float B[SIZE];
  float coefficients[SIZE][SIZE];

  const float eps = 0.001;
  int iterationsCounter = 1;
  bool breakpoint = true;
  float el;

  for (int i = 0; i < SIZE; i++) {
    el = M[i][i];
    for (int u = 0, k = 0; u < SIZE + 1; u++) {
      if (u == i) continue;
      coefficients[i][k] = M[i][u] / el;
      if (u != SIZE) coefficients[i][k++] *= -1;
    }
  }

  for (int i = 0; i < SIZE; i++)
  for (int u = 0; u < SIZE - 1; u++) {
    float el = coefficients[i][u];
    B[i] += (el > 0) ? el : -el;
  }

  for (int i = 0; i < SIZE; i++) {
    X[i] = coefficients[i][SIZE - 1];
    X_prev[i] = X[i];
  }

  const float bMax = maxNum(B);
  const float div = ((1.0 - bMax) / bMax) * eps;


//show the equations set
  printf("Equations set:\n");
  for (int i = 0; i < SIZE; i++) {
    for (int u = 0; u < SIZE; u++) printf("%f\t", coefficients[i][u]);
    printf("\n");
  }

//main iteration
  while (breakpoint) {
    printf("\nIteration %d:\n", iterationsCounter);
    showX("\tPrevious X\'es:\n\t\t", X_prev);
    showX("\tCurrent X\'es:\n\t\t", X);

    for (int i = 0; i < SIZE; i++) {
      X_prev[i] = X[i];
      X[i] = 0;
    }
        
    for (int i = 0; i < SIZE; i++) {
      for (int u = 0; u < SIZE - 1; u++)
        X[i] += X_prev[i] * coefficients[i][u];
      X[i] += coefficients[i][SIZE - 1];
    }

    breakpoint = !imparity(X_prev, X, div);
    iterationsCounter++;
  }

  printf("\n[FINAL RESULT]:\t");
  showX("", X);

  return EXIT_SUCCESS;
}
