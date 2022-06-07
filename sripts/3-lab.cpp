//chord method
#include <iostream>
#include <math.h>

//f(x) = 1 - 0.5*x^2*ln(x) + 0.3*sqrt(x)

const int A = 2, B = 3;
const float eps = pow(10, -6);

float f(const float& x) {
  return ( 1 - 0.5 * pow(x, 2) * log(x) + 0.3 * sqrt(x) );
}

float chordMethod(const float& x1, const float& x2) {
  return ( x1 - (f(x1) / (f(x2) - f(x1))) * (x2 - x1) );
}

void showResult(const int& n, const float& x, const float& f_x, const float& delta = 0.0) {
  printf("Iteration %d:\n", n);
  printf("\tx = %.7f;\n", x);
  printf("\tf(x) = %.8f;\n", f_x);
  printf("\t|x - x_prev| = %.8f;\n", delta);
}

int main() {
  float x = B, x_prev = A;
  bool breakpoint = true;
  int iterationsCounter = 0;
  showResult(iterationsCounter++, x_prev, f(x_prev));

  while (breakpoint) {
    float _x = x;
    x = chordMethod(x_prev, x);
    x_prev = _x;
    const float delta = abs(x - x_prev);
    showResult(iterationsCounter++, x, f(x), delta);
    breakpoint = eps < abs(x - x_prev);
  }

  printf("\n[FINAL RESULT]: %.7f\n",x);
    
  return EXIT_SUCCESS;
}
