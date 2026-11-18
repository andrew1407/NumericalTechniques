//table function approximation (least square technique)
#include <iostream>
#include <math.h>

#define SIZE 9

const float X[SIZE] =  { 0.235, 0.569, 0.903, 1.237, 1.571, 1.906, 2.240, 2.574, 2.908 };
const float Y[SIZE] =  { 1.082, 1.563, 2.182, 3.078, 3.816, 4.500, 5.443, 6.209, 7.082 };

float det3(const float m[3][3]) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

//Cramer's rule for the 3x3 normal system {a, b, c | rhs}
void solve3(const float S[3][4], float& a, float& b, float& c) {
  float m[3][3], res[3];
  for (int i = 0; i < 3; i++)
    for (int u = 0; u < 3; u++) m[i][u] = S[i][u];
  const float det = det3(m);

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++)
      for (int u = 0; u < 3; u++) m[i][u] = (u == k) ? S[i][3] : S[i][u];
    res[k] = det3(m) / det;
  }
  a = res[0]; b = res[1]; c = res[2];
}

void func1Sys(float& a, float& b, float& c) {
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

  const float S[3][4] = {
    { sumX4, sumX3, sumX2, sumYX2 },
    { sumX3, sumX2, sumX, sumYX },
    { sumX2, sumX, SIZE + 0.0f, sumY }
  };
  solve3(S, a, b, c);
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

void func2Sys(float& a, float& b, float& c) {
  float sumX2 = 0, sumX = 0, sumX2exp = 0, sumXexp = 0, sumExp = 0, sumExp2 = 0;
  float sumYX = 0, sumYExp = 0, sumY = 0; 
  for (int i = 0; i < SIZE; i++) {
    sumX2 += powf(X[i], -2);
    sumX += powf(X[i], -1);
    sumY += Y[i];
    sumXexp += expf(X[i]) / X[i];
    sumX2exp += expf(X[i]) / powf(X[i], 2);
    sumExp += expf(X[i]);
    sumExp2 += powf(expf(X[i]), 2);
    sumYExp += Y[i] * expf(X[i]);
    sumYX += Y[i] / X[i];
  }

  printf("Function 2 (normal system):\n");
  printf("a * %.4f + b * %.4f + c * %.4f = %.4f\n", sumX2, sumXexp, sumX, sumYX);
  printf("a * %.4f + b * %.4f + c * %.4f = %.4f\n", sumXexp, sumExp2, sumExp, sumYExp);
  printf("a * %.4f + b * %.4f + c * %.4f = %.4f\n", sumX, sumExp, SIZE + 0.0, sumY);

  const float S[3][4] = {
    { sumX2, sumXexp, sumX, sumYX },
    { sumXexp, sumExp2, sumExp, sumYExp },
    { sumX, sumExp, SIZE + 0.0f, sumY }
  };
  solve3(S, a, b, c);
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
  float a1, b1, c1;
  func1Sys(a1, b1, c1);
  printf("\n");
  discrepasny1(a1, b1, c1);
  printf("\n");

  float a2, b2, c2;
  func2Sys(a2, b2, c2);
  printf("\n");
  discrepasny2(a2, b2, c2);

  return EXIT_SUCCESS;
}
