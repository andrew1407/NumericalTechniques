//table function interpolation
#include <iostream>
#include <math.h>
#include <string>

#define SIZE 5
#define String std::string
#define toString std::to_string

float X[SIZE] = { 0.092, 0.772, 1.385, 2.108, 2.938 };
float Y[SIZE] = { 3.161, 1.357, -0.158, -0.129, -4.438 };

class PolynumFn {
private:
  float *X, *Y;
  float DivDif[SIZE];

  void GenerateDivDif() {
    float Y_arr1[SIZE], Y_arr2[SIZE];
    for (int i = 0; i < SIZE; i++) 
      Y_arr1[i] = Y_arr2[i] = Y[i];

    //first-state output
    printf("Divided differences:\n");
    for (int u = 0; u < SIZE; u++) printf("%.4f  ", Y_arr2[u]);
    printf("\n");
    int shift = 1;
    for (int i = 1; i < SIZE; i++) {
      for (int u = i; u < SIZE; u++)
        Y_arr1[u] = (Y_arr2[u] - Y_arr2[u - 1]) / (X[u] - X[u - shift]);
      shift++;
      for (int u = i; u < SIZE; u++) Y_arr2[u] = Y_arr1[u];
      //next-state output
      for (int u = i; u < SIZE; u++) printf("%.4f  ", Y_arr2[u]);
      printf("\n");
    }

    for (int i = 0; i < SIZE; i++) DivDif[i] = Y_arr1[i];
  }

public:
  PolynumFn(float* X_args, float* Y_args): X(X_args), Y(Y_args) { }

  void GenerateFinDif() {
    const int len = SIZE - 1;
    const float h = pow(SIZE, -1);
    float X_dif[SIZE];

    printf("Finite differences:\n");
    for (int i = len; i >= 0; i--) {
      for (int k = 0; k <= i; k++) {
        X_dif[k] = X[len - i] + (h * k);
        printf("%.4f  ", X_dif[k]);
      }
      printf("\n");
    }
  }

  float LxFunction(const float& x) {
    float dif_PX[SIZE];
    float koefs_Px[SIZE];
    float koef, res;
    String Px_strList[SIZE];
    String Lx_strList[SIZE];
    String p_str;

    for (int i = 0; i < SIZE; i++) {
      dif_PX[i] = koef = 1;
      p_str = "";
      
      for (int u = 0; u < SIZE; u++) {
        if (u == i) continue;
        koef *= (X[i] - X[u]);
        dif_PX[i] *= (x - X[u]);
        p_str += "(x - " + toString(X[u]) + ")";
      }

      koef = powf(koef, -1);
      float mpl = koef * Y[i];
      dif_PX[i] *= mpl;
      Px_strList[i] = toString(koef) + p_str;
      Lx_strList[i] = toString(mpl) + p_str;
      res += dif_PX[i];
    }

    //output
    printf("P(x):\n");
    for (int i = 0; i < SIZE; i++) {
      printf("x%d:\t", i);
      printf((Px_strList[i] + "\n").c_str());
    }

    printf("\nL(x):\n");

    for (int i = 0; i < SIZE; i++)
      switch (i) {
        case 0:
          printf((Lx_strList[i] + " +\n").c_str());
          break;
        case SIZE - 1:
          printf(("+ " + Lx_strList[i]).c_str());
          break;
        default:
          printf(("+ " + Lx_strList[i] + " +\n").c_str());
          break;
      }

    printf("\n");

    //RETURN L(X) VALUE
    return res;
  }

  float NxFunction(const float& x) {
    GenerateDivDif();
    printf("\n");
    float res = DivDif[0];
    float mlp;
    String Nx_strList[SIZE];
    String Nx_str;
    Nx_strList[0] = toString(res);

    for (int i = 1; i < SIZE; i++) {
      mlp = DivDif[i];
      Nx_str = "";
      for (int u = 0; u < i; u++) {
        mlp *= (x - X[u]);
        Nx_str += "(x - " + toString(X[u]) +")";
      }

      Nx_strList[i] = toString(DivDif[i]) + Nx_str;
      res += mlp;
    }
        
    //output
    Nx_strList[1] = Nx_strList[0] + " + " + Nx_strList[1];
    printf("N(x):\n");
    for (int i = 1; i < SIZE; i++)
      switch (i) {
        case 1:
          printf((Nx_strList[i] + " + \n").c_str());
          break;
        case SIZE - 1:
          printf(("+ " + Nx_strList[i] + "\n").c_str());
          break;
        default:
          printf(("+ " + Nx_strList[i] + " +\n").c_str());
          break;
      }

    //RETURN N(X) RESULT
    return res;        
  }
};

int main()
{
  PolynumFn fn(X, Y);
  float res;
  res = fn.LxFunction(X[1] + X[2]);
  printf("\nL(x1 + x2) = %.4f\n\n", res);
  fn.GenerateFinDif();
  printf("\n");
  res = fn.NxFunction(X[1] + X[2]);
  printf("\nN(x1 + x2) = %.4f\n", res);

  return EXIT_SUCCESS;
}