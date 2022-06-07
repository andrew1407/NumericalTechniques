//table function approximation (least square technique)
#include <iostream>
#include <math.h>

#define SIZE 9

const float X[SIZE] =  { 0.235, 0.569, 0.903, 1.237, 1.571, 1.906, 2.240, 2.574, 2.908 };
const float Y[SIZE] =  { 1.082, 1.563, 2.182, 3.078, 3.816, 4.500, 5.443, 6.209, 7.082 };

void func1Sys() {
  float sumX4 = 0, sumX3 = 0, sumX2 = 0, sumX = 0;
  float sumYX2 = 0, sumYX = 0, sumY = 0; 
  for (int i = 0; i < SIZE; i++) {
    sumX4 += powf(X[i], -4);
    sumX3 += powf(X[i], -3);
    sumX2 += powf(X[i], -2);
    sumX += powf(X[i], -1);
    sumYX2 += Y[i] / powf(X[i], 2);
    sumYX += Y[i] /(X[i]);
    sumY += Y[i];
  }

  printf("Function 1 (normal system):\n");
  printf("a * %.4f + b * %.4f + c * %.4f = %.4f\n", sumX4, sumX3, sumX2, sumYX2);
  printf("a * %.4f + b * %.4f + c * %.4f = %.4f\n", sumX3, sumX2, sumX, sumYX);
  printf("a * %.4f + b * %.4f + c * %.4f = %.4f\n", sumX2, sumX, SIZE + 0.0, sumY);

}

void discrepasny1(const float& a, const float& b, const float& c) {
  printf("Function 1:\n\ty = %.4f * (1 / x^2) + %.4f * (1 / x) + %.4f\n", a, b, c);
  float S = 0;
  for (int i = 0; i < SIZE; i++) {
    const float f = a * powf(X[i], -2) + b * (powf(X[i], -1)) + c;
    S += powf((f - Y[i]), 2);
  }
  printf("\tS = %.4f\n\tdel = %.4f\n", S, sqrtf(S));
}

void func2Sys() {
  float sumX2 = 0, sumX = 0, sumX2exp = 0, sumXexp = 0, sumExp = 0, sumExp2 = 0;
  float sumYX = 0, sumYExp = 0, sumY = 0; 
  for (int i = 0; i < SIZE; i++) {
    sumX2 += powf(X[i], -2);
    sumX += powf(X[i], -1);
    sumYX += Y[i] / powf(X[i], -1);
    sumY += Y[i];
    sumXexp += expf(X[i]) / X[i];
    sumX2exp += expf(X[i]) / powf(X[i], 2);
    sumExp += expf(X[i]);
    sumExp2 += powf(expf(X[i]), 2);
    sumY += Y[i];
    sumYExp += Y[i] * expf(X[i]);
    sumYX += Y[i] / X[i];
  }

  printf("Function 2 (normal system):\n");
  printf("a * %.4f + b * %.4f + c * %.4f = %.4f\n", sumX2, sumXexp, sumX, sumYX);
  printf("a * %.4f + b * %.4f + c * %.4f = %.4f\n", sumXexp, sumExp2, sumExp, sumYExp);
  printf("a * %.4f + b * %.4f + c * %.4f = %.4f\n", sumX, sumExp, SIZE + 0.0, sumY);
}

void discrepasny2(const float& a, const float& b, const float& c) {
  printf("Function 2:\n\ty = %.4f * (1 / x) + %.4f * exp(x) + %.4f\n", a, b, c);
  float S = 0;
  for (int i = 0; i < SIZE; i++) {
    const float f = a * powf(X[i], -1) + b * expf(X[i]) + c;
    S += powf((f - Y[i]), 2);
  }
  printf("\tS = %.4f\n\tdel = %.4f\n", S, sqrtf(S));
}

int  main()
{
  const float a1 = 0.1177, b1 = 1.7591, c1 = 1.681;
  func1Sys();
  printf("\n");
  discrepasny1(a1, b1, c1);
  printf("\n");

  const float a2 = -0.4091, b2 = -0.5904, c2 = 12.2568;
  func2Sys();
  printf("\n");
  discrepasny2(a2, b2, c2);

  return EXIT_SUCCESS;
}
