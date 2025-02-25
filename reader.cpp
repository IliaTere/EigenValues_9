#include "func.h"
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/resource.h>

double get_full_time()
{
    struct timeval buf;
    gettimeofday(&buf, NULL);
    return buf.tv_sec + buf.tv_usec / 1.e6;
}

double get_CPU_time()
{
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1.e6;
}

bool isNumber(std::string& str)
{
    std::string::iterator it = std::begin(str);
    while (it != str.end() && std::isdigit(*it)) {
        it++;
    }
    return !str.empty() && it == str.end();
}

bool is_double(const std::string& s) {
    if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+') && (s[0] != '.'))) return false;
    char* p;
    strtod(s.c_str(), &p);
    return (*p == 0);
}

int toDouble(const char* str, double* ptr)
{
    char* e;

    errno = 0;
    *ptr = strtod(str, &e);

    if (!errno && *e == '\0')
        return 0;
    else
        return -1;
}

int read_ff(const std::string& filename, double* result, int n) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Error opening file: " << filename << std::endl;
        return -1;
    }
    std::string line;
    int i = 0;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
            if (is_double(token)) {
                if (i > n-1) {
                    printf("Incorrect size\n");
                    return -3;
                }
                result[i] = std::stod(token);
                i++;
            } else {
                std::cout << "Invalid double value: " << token << std::endl;
                return -2;
            }
        }
    }
    if (i != n) {
        std::cout << "Invalid length" << std::endl;
        return -3;
    }
    return 0;
}

void PrintDouble(double* matrix, int n, int r) {
	for(int i = 0; i < std::min(n, r); i++) {
        for(int j = 0; j < std::min(n, r); j++) {
            printf("%10.3e ", matrix[i*n+j]);
        }
        std::cout << std::endl;
    }
	printf("\n");
}

void input_formula(double *A, int s, int n)
{
    for(int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= n; j++)
        {
            switch(s) {
                case 1:
                    A[(i - 1) * n + j - 1] = n - max(i, j) + 1;
                    break;

                case 2:
                    if(i == j)
                        A[(i - 1) * n + j - 1] = 2.;
                    else if(abs(i - j) == 1)
                        A[(i - 1) * n + j - 1] = -1;
                    else
                        A[(i - 1) * n + j - 1] = 0;
                    break;

                case 3:
                    if(i == j and i != n - 1)
                        A[(i - 1) * n + j - 1] = 1;
                    else if(j == n - 1)
                        A[(i - 1) * n + j - 1] = i;
                    else if(i == n - 1)
                        A[(i - 1) * n + j - 1] = j;
                    else
                        A[(i - 1) * n + j - 1] = 0;
                    break;

                case 4:
                    A[(i - 1) * n + j - 1] = 1. / (i + j - 1);
                    break;
            }
        }
    }
}
bool is_symmetric(double *A, int n)
{
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(fabs(A[i * n + j] - A[j * n + i]) > EPS) {
                return false;
            }
        }
    }
    return true;
}
void get_column(double *A, double *a1, int n, int k)
{
    memset(a1, 0, sizeof(double) * n);

    for(int i = k + 1; i < n; i++)
        a1[i] = A[i * n + k];
}
double vector_norm(double *x, int n)
{
    double sum = 0;
    for(int i = 0; i < n; i++)
        sum += x[i] * x[i];

    return sqrt(sum);
}
int vector_division(double *a1, int n, double alpha)
{
    if(abs(alpha) < 1e-100) return -1;

    for(int i = 0; i < n; i++)
        a1[i] = a1[i] / alpha;

    return 0;
}

void UAUt(Args *arg)
{
    double *a = arg->a, *b = arg->b, *c = arg->c, *xk = arg->xk, *y = arg->y, *z = arg->z;
    int n = arg->n;

    matrix_product(a, xk, y, n, n, 1);
    matrix_product(xk, y, z, 1, n, 1);
    double xy = z[0];

    for(int i = 0; i < n; i++)
        z[i] = xk[i] * xy;

    matrixSubtraction(y, z, z, 1, n); 
    vector_multiplying(z, n, 2);
    matrix_product(z, xk, c, n, 1, n); 
    matrix_product(xk, z, b, n, 1, n); 
    matrixSubtraction(a, c, a, n, n);
    matrixSubtraction(a, b, a, n, n);
}