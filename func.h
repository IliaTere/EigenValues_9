#include <iostream>
#include <cstring>
#include <cmath>
#include <fstream>
#include <sstream>

#define EPS 1e-16

using namespace std;

class Args
{
public:
    int n = 0;
    double *A = nullptr;
    double *x = nullptr;
    double eps = 0;

    double *a = nullptr;
    double *b = nullptr;
    double *c = nullptr;
    double *xk = nullptr;
    double *y = nullptr;
    double *z = nullptr;

    double norm = 0;

    int its = 0;
    double three_diagonal_time = 0;
    double eigenvalues_time = 0;

    ~Args()
    {
        delete[] a;
        delete[] b;
        delete[] c;
        delete[] xk;
        delete[] y;
        delete[] z;
    }
};

double get_full_time();
double get_CPU_time();
bool isNumber(std::string& str);
bool is_double(const std::string& s);
int toDouble(const char* str, double* ptr);
int read_ff(const std::string& filename, double* result, int n);
void PrintDouble(double* matrix, int n, int r);
void input_formula(double *A, int s, int n);
bool is_symmetric(double *A, int n);
int solve(Args *arg);
double m_norm(double *A, int n);
void get_column(double *A, double *a1, int n, int k);
double vector_norm(double *x, int n);
int vector_division(double *a1, int n, double alpha);
void UAUt(Args *arg);
void matrix_product(double *A, double* B, double* C, int n, int s, int m);
int matrixSubtraction(double* A1, double *A2, double *C, int n, int m);
void vector_multiplying(double *a1, int n, double alpha);
double r1(double *A, int n, double *x);
double r2(double *A, int n, double *x);
