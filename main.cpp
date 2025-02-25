#include "func.h"
#include <fenv.h>

int main(int argc, char **argv)
{
    feenableexcept(FE_DIVBYZERO | FE_INVALID
        | FE_OVERFLOW | FE_UNDERFLOW);
    if (argc < 5 || argc > 6) {
        printf("error: ./a.out n m eps k\n");
        return -1;
    }
    for (int i = 1; i < 5; i++) {
        std::string str(argv[i]);
        if (isNumber(str) == false && i != 3) {
            std::cout << "error: Invalid argument: " << str << std::endl;
            return -4;
        }
    }
    std::string s1(argv[1]);
    std::string s2(argv[2]);
    std::string s3(argv[3]);
    int s;
    int n = std::stoi(s1); // Размерность матрицы
    int m = std::stoi(s2); // Размерность блока
    double eps;
    if(toDouble(argv[3], &eps) == -1)
        {
            cout << "Can't read parameter eps." << endl;
            return -1;
        }
    // if ( m == 0  || m > n) {
    //     printf("invalid block\n");
    //     return -8;
    // }
    double* matr = new double[n*n];
    if (strcmp(argv[4],"0") == 0) {
        if (argc != 6) {
            std::cout << "error: File not found" << std::endl;
            delete[] matr;
            return -2;
        }
        std::string name(argv[5]);
		int t = read_ff(name , matr , n*n);
        if(t != 0) {
            
            delete[] matr;
            return -5;
        }
        std::string tmp(argv[4]);
        s = stoi(tmp);
    }
    if ((strcmp(argv[4],"0") != 0)) 
    {
        if (argc > 5) {
            std::cout << "error: To many(few) arguments" << std::endl;
            delete[] matr;
            return -6;
        }
        std::string tmp(argv[4]);
        s = stoi(tmp);
        if (s < 1 || s > 4) {
            std::cerr << "error: Parametr s is not a valid number" << std::endl;
            delete[] matr;
            return -7;
        }
        input_formula(matr, s, n);
    }
    double* x_vals = new double[n];
    if(m > n ) {
        m = n;
    }
    PrintDouble(matr, n, m);
    
    if (is_symmetric(matr, n) != 1) {
        cout << "Matrix is not symmetric." << endl;
        delete[] matr;
        delete[] x_vals;
        return EXIT_FAILURE;
    }
    
    Args arguments;
    arguments.n = n;
    arguments.A = matr;
    arguments.x = x_vals;
    arguments.eps = eps;
    
    arguments.a = new double[n*n];
    arguments.b = new double[n*n];
    arguments.c = new double[n*n];
    arguments.xk = new double[n];
    arguments.y = new double[n];
    arguments.z = new double[n];
    
    std::copy(matr, matr + n*n, arguments.a);
    
    int calc_result = solve(&arguments);
    double residual1 = 0.0, residual2 = 0.0;
    bool success = (calc_result != -1);
    
    if (success) {
        if (m > 0) {
            cout << "Result:" << endl;
            for (int idx = 0; idx < m; ++idx)
                printf("%10.3e ", x_vals[idx]);
            putchar('\n');
        }
        residual1 = r1(matr, n, x_vals);
        residual2 = r2(matr, n, x_vals);
    } else {
        cerr << "Unsupported matrix type for this method." << endl;
    }
    
    int total_iter = arguments.its;
    double timing1 = arguments.three_diagonal_time, timing2 = arguments.eigenvalues_time;
    
    printf("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Timing1 = %.2f Timing2 = %.2f\n",
           argv[0], residual1, residual2, total_iter, total_iter / n, timing1, timing2);
    

    delete[] matr;
    delete[] x_vals;
    
    return EXIT_SUCCESS;
}