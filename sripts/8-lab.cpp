//numerical integratiom
#include <iostream>
#include <math.h>

class Inequation {
private:
  //input values
  const int A = 0, B = 2;
  const float eps = pow(10, -4);
  const float h = 0.1;
  const int x0 = 0, y0 = 1;
  float x = x0, y = y0, h0 = h;

  //other values;
  float RK_Yvalues[20];

  //methods |PRIVATE|

  float _fn(const float& x, const float& y) {
    return 0.5 * (x + 1) * expf(x) * pow(y, 2) - x * y;
  }

  float _fnExact(const float& x, const float& y) {
    const float part1 = 2 * expf(-0.5 * pow(x, 2));
    const float part2 = 4.13273 * erff(0.707107 - 0.707107 * x);
    const float part3 = part2 + expf(x * (1 - 0.5 * x)) - 1.82137;
    return part1 / part3;
  }

  float _Y_k(const float& x, const float& h, const float& y) {
    const float F1 = _fn(x, y);
    const float F2 = _fn(x + h * 0.5, y + h * 0.5 * F1);
    const float F3 = _fn(x + h * 0.5, y + h * 0.5 * F2);
    const float F4 = _fn(x + h * 0.5, y + h * F3);
    float res = F1 + 2 * F2 + 2 * F3 + F4;
    res *= (h / 6);
    res += y;
    return res;
  }

//methods |PUBLIC|
public:

  void defineStep() {
    printf("[STEP SEARCHING]\n\n");
    const float greater = pow(10, -10);
    float y1, y2, _y2, _x = x;
    bool breakpoint = true;
    float dif;

    while (breakpoint) {
      y1 = _Y_k(_x + h0, h0, y0);
      y2 = _Y_k(_x + 2 * h0, h0, y0);
      _y2 = _Y_k(_x + 2 * h0, 2 * h0, y0);
      dif = abs(y2 - _y2);
      _x += h0;

      printf("y1 = %.6f\n", y1);
      printf("y2 = %.6f\n", y2);
      printf("~y2 = %.6f\n", _y2);
      printf("|y2 - ~y2| = %.6f\n", dif);

      breakpoint = (dif < eps);
      if (breakpoint) h0 *= 2;

      printf("h0 = %.6f\n\n", h0);
    }

    //after step searching
    int n = (B - A) / h0;
    if (n % 2) h0 = (B - A) / ++n;

    y1 = _Y_k(x + h0, h0, y0);
    y2 = _Y_k(x + 2 * h0, h0, y0);
    _y2 = _Y_k(x + 2 * h0, 2 * h0, y0);
    dif = abs(y2 - _y2);

    printf("[AFTER STEP SEARCHING]\n\n");
    printf("y1 = %.6f\n", y1);
    printf("y2 = %.6f\n", y2);
    printf("~y2 = %.6f\n", _y2);
    printf("|y2 - ~y2| = %.6f\n", dif);

    printf("\nn = %d\n", n);
    printf("h0 = %f\n", h0);
    printf("|y2 - ~y2| = %.6f\n", dif);
  }

  void MethodRungeKutta() {
    float deltaMax = 0, del;
    float yi = y0, _yi = y0;
    printf("\n[RUNGE-KUTTA(IV) METHOD]\n\n");

    for (float xi = A, i = 1; xi <= B; xi += h0, i++) {
      yi = _Y_k(xi, h0, yi);
      _yi = _Y_k(xi, 2 * h0, yi);
      del = abs(yi - _yi);
      if (del > deltaMax) deltaMax = del;
      RK_Yvalues[(int)i] = yi;

      printf("Iteration %d:\n", (int)i);
      printf("xi = %.6f\n", xi);
      printf("yi = %.6f\n", yi);
      printf("~yi = %.6f\n", _yi);
      printf("|yi - ~yi| = %.6f\n\n", del);
    }

    printf("Max delta |yi - ~yi| = %.6f\n", deltaMax);
  }

  void EulerMethod() {
    printf("\n[EULER METHOD]\n\n");
    float yi = y, _yi = y;
    float del, deltaMax = 0;

    for (float xi = A, i = 1; xi <= B; xi += h0, i++) {
      yi += h0 * _Y_k(xi, h0, yi);
      _yi += 2 * h0 * _Y_k(xi, h0, yi);
      del = abs(yi - _yi);
      if (del > deltaMax) deltaMax = del;

      printf("Iteration %d:\n", (int)i);
      printf("xi = %.6f\n", xi);
      printf("yi = %.6f\n", yi);
      printf("~yi = %.6f\n", _yi);
      printf("|yi - ~yi| = %.6f\n\n", del);
    }

    printf("Max delta |yi - ~yi| = %.6f\n", deltaMax);
  }

  void ExactInequation() {
    printf("\n[EXACT CALCULATING]\n\n");
    float yi = y, yiKG = y;
    float del, deltaMax = 0;

    for (float xi = A, i = 1; xi <= B; xi += h0, i++) {
      yi = _fnExact(xi, yi);
      yiKG = RK_Yvalues[(int)i];
      del = abs(yi - yiKG);
      if (del > deltaMax) deltaMax = del;

      printf("Iteration %d:\n", (int)i);
      printf("xi = %.6f\n", xi);
      printf("yi = %.6f\n", yi);
      printf("yi_RK = %.6f\n", yiKG);
      printf("|yi - ~yi| = %.6f\n\n", del);
    }

    printf("Max delta |yi - ~yi| = %.6f\n", deltaMax);
  }

};

int main() {
  Inequation inequation;
  inequation.defineStep();
  printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  inequation.MethodRungeKutta();
  printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  inequation.EulerMethod();
  printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  inequation.ExactInequation();

  return EXIT_SUCCESS;
}
