//numerical integratiom
#include <iostream>
#include <math.h>

//input values
const float A = 0.1, B = 1.1;       //interval
const float h = 0.125;              //step value
const int n = 8;                    //number of steps

float fn(const float &x) {
  return sqrtf(x + 1) * log10f(x + 1);
}

void TrapezeFn() {
  float h = ::h;
  float res1 = fn(A) / 2;
  for (int i = 1; i < n; i++) res1 += fn(A + i * h);
  res1 += fn(B) / 2;
  res1 *= h;

  h *= 2;
  float res2 = fn(A) / 2;
  const int m = n / 2;
  for (int i = 1; i < m; i++) res2 += fn(A + i * h);
  res2 += fn(B) / 2;
  res2 *= h;

  const float del = abs((res1 - res2) / 3);

  printf("Trapeze method:\n");
  printf("\tI(h) = %.6f:\n", res1);
  printf("\tI(2h) = %.6f:\n", res2);
  printf("\tI-I^tp(h) = %.6f:\n", del);

}

void SimpsonFn() {
  float h = ::h;
  float res1 = (fn(A) + fn(A + n * h)) / 2;
  for (int i = 1; i < n; i++)
    res1 += fn(A + i * h) * ((i % 2) ? 2 : 1);
  res1 *= h * 0.6667;
    
  h *= 2;
  int m = n / 2;
  float res2 = (fn(A) + fn(A + m * h)) / 2;
  for (int i = 1; i < m; i++)
    res2 += fn(A + i * h) * ((i % 2) ? 2 : 1);
  res2 *= h * 0.6667;

  float del = abs((res1 - res2) / 15);

  printf("Simpson\'s method:\n");
  printf("\tI(h) = %.8f:\n", res1);
  printf("\tI(2h) = %.8f:\n", res2);
  printf("\tI-I^tp(h) = %.10f:\n", del);
}

int main() {
  TrapezeFn();
  SimpsonFn();

  return EXIT_SUCCESS;
}
