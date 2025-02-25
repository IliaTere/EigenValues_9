#include "func.h"

void matrix_product(double *A, double* B, double* C, int n, int s, int m)
{
    int k1, l1, k2, l2, k3, l3, i, j, t;
    double t00, t01, t02, t10, t11, t12, t20, t21, t22;
    k1 = n/3;
    k2 = s/3;
    k3 = m/3;
    l1 = n - k1*3;
    l2 = s - k2*3;
    l3 = m - k3*3;

    for(i = 0; i < k1; i++){
        for(j = 0; j < k3; j++){
            t00 = 0;
            t01 = 0;
            t02 = 0;
            t10 = 0;
            t11 = 0;
            t12 = 0;
            t20 = 0;
            t21 = 0;
            t22 = 0;
            for(t = 0; t < k2; t++){
                t00 += A[3*i*s+3*t]*B[3*t*m+3*j] + A[3*i*s+3*t+1]*B[(3*t+1)*m+3*j] + A[3*i*s+3*t+2]*B[(3*t+2)*m+3*j];
                t01 += A[3*i*s+3*t]*B[3*t*m+3*j+1] + A[3*i*s+3*t+1]*B[(3*t+1)*m+3*j+1] + A[3*i*s+3*t+2]*B[(3*t+2)*m+3*j+1];
                t02 += A[3*i*s+3*t]*B[3*t*m+3*j+2] + A[3*i*s+3*t+1]*B[(3*t+1)*m+3*j+2] + A[3*i*s+3*t+2]*B[(3*t+2)*m+3*j+2];
                t10 += A[(3*i+1)*s+3*t]*B[3*t*m+3*j] + A[(3*i+1)*s+3*t+1]*B[(3*t+1)*m+3*j] + A[(3*i+1)*s+3*t+2]*B[(3*t+2)*m+3*j];
                t11 += A[(3*i+1)*s+3*t]*B[3*t*m+3*j+1] + A[(3*i+1)*s+3*t+1]*B[(3*t+1)*m+3*j+1] + A[(3*i+1)*s+3*t+2]*B[(3*t+2)*m+3*j+1];
                t12 += A[(3*i+1)*s+3*t]*B[3*t*m+3*j+2] + A[(3*i+1)*s+3*t+1]*B[(3*t+1)*m+3*j+2] + A[(3*i+1)*s+3*t+2]*B[(3*t+2)*m+3*j+2];
                t20 += A[(3*i+2)*s+3*t]*B[3*t*m+3*j] + A[(3*i+2)*s+3*t+1]*B[(3*t+1)*m+3*j] + A[(3*i+2)*s+3*t+2]*B[(3*t+2)*m+3*j];
                t21 += A[(3*i+2)*s+3*t]*B[3*t*m+3*j+1] + A[(3*i+2)*s+3*t+1]*B[(3*t+1)*m+3*j+1] + A[(3*i+2)*s+3*t+2]*B[(3*t+2)*m+3*j+1];
                t22 += A[(3*i+2)*s+3*t]*B[3*t*m+3*j+2] + A[(3*i+2)*s+3*t+1]*B[(3*t+1)*m+3*j+2] + A[(3*i+2)*s+3*t+2]*B[(3*t+2)*m+3*j+2];
            }
            for(t = 0; t < l2; t++){
                t00 += A[3*i*s+t+3*k2]*B[(t+3*k2)*m+3*j];
                t01 += A[3*i*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+1];
                t02 += A[3*i*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+2];
                t10 += A[(3*i+1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j];
                t11 += A[(3*i+1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+1];
                t12 += A[(3*i+1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+2];
                t20 += A[(3*i+2)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j];
                t21 += A[(3*i+2)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+1];
                t22 += A[(3*i+2)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+2];
            }
            C[3*i*m+3*j] = t00;
            C[3*i*m+3*j+1] = t01;
            C[3*i*m+3*j+2] = t02;
            C[(3*i+1)*m+3*j] = t10;
            C[(3*i+1)*m+3*j+1] = t11;
            C[(3*i+1)*m+3*j+2] = t12;
            C[(3*i+2)*m+3*j] = t20;
            C[(3*i+2)*m+3*j+1] = t21;
            C[(3*i+2)*m+3*j+2] = t22;
        }

        for(j = 0; j < l3; j++){
            t00 = 0;
            t10 = 0;
            t20 = 0;
            for(t = 0; t < k2; t++){

                t00 += A[3*i*s+3*t]*B[3*t*m+j+3*k3] + A[3*i*s+3*t+1]*B[(3*t+1)*m+j+3*k3] + A[3*i*s+3*t+2]*B[(3*t+2)*m+j+3*k3];
                t10 += A[(3*i+1)*s+3*t]*B[3*t*m+j+3*k3] + A[(3*i+1)*s+3*t+1]*B[(3*t+1)*m+j+3*k3] + A[(3*i+1)*s+3*t+2]*B[(3*t+2)*m+j+3*k3];
                t20 += A[(3*i+2)*s+3*t]*B[3*t*m+j+3*k3] + A[(3*i+2)*s+3*t+1]*B[(3*t+1)*m+j+3*k3] + A[(3*i+2)*s+3*t+2]*B[(3*t+2)*m+j+3*k3];

            }
            for(t = 0; t < l2; t++){
                t00 += A[3*i*s+t+3*k2]*B[(t+3*k2)*m+j+3*k3];
                t10 += A[(3*i+1)*s+(t+3*k2)]*B[(t+3*k2)*m+j+3*k3];
                t20 += A[(3*i+2)*s+(t+3*k2)]*B[(t+3*k2)*m+j+3*k3];
            }
            C[3*i*m+j+3*k3] = t00;
            C[(3*i+1)*m+j+3*k3] = t10;
            C[(3*i+2)*m+j+3*k3] = t20;
        }
    }

    for(i = 0; i < l1; i++){
        for(j = 0; j < k3; j++){
            t00 = 0;
            t01 = 0;
            t02 = 0;
            for(t = 0; t < k2; t++){
                t00 += A[(i+3*k1)*s+3*t]*B[3*t*m+3*j] + A[(i+3*k1)*s+3*t+1]*B[(3*t+1)*m+3*j] + A[(i+3*k1)*s+3*t+2]*B[(3*t+2)*m+3*j];
                t01 += A[(i+3*k1)*s+3*t]*B[3*t*m+3*j+1] + A[(i+3*k1)*s+3*t+1]*B[(3*t+1)*m+3*j+1] + A[(i+3*k1)*s+3*t+2]*B[(3*t+2)*m+3*j+1];
                t02 += A[(i+3*k1)*s+3*t]*B[3*t*m+3*j+2] + A[(i+3*k1)*s+3*t+1]*B[(3*t+1)*m+3*j+2] + A[(i+3*k1)*s+3*t+2]*B[(3*t+2)*m+3*j+2];

            }
            for(t = 0; t < l2; t++){
                t00 += A[(i+3*k1)*s+t+3*k2]*B[(t+3*k2)*m+3*j];
                t01 += A[(i+3*k1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+1];
                t02 += A[(i+3*k1)*s+(t+3*k2)]*B[(t+3*k2)*m+3*j+2];
            }
            C[(i+3*k1)*m+3*j] = t00;
            C[(i+3*k1)*m+3*j+1] = t01;
            C[(i+3*k1)*m+3*j+2] = t02;
        }
        for(j = 0; j < l3; j++){
            t00 = 0;
            for(t = 0; t < k2; t++){
                t00 += A[(i+3*k1)*s+3*t]*B[3*t*m+j+3*k3] + A[(i+3*k1)*s+3*t+1]*B[(3*t+1)*m+j+3*k3] + A[(i+3*k1)*s+3*t+2]*B[(3*t+2)*m+j+3*k3];
            }
            for(t = 0; t < l2; t++){
                t00 += A[(i+3*k1)*s+t+3*k2]*B[(t+3*k2)*m+j+3*k3];
            }
            C[(i+3*k1)*m+j+3*k3] = t00;
        }
    }
}
int matrixSubtraction(double* A1, double *A2, double *C, int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
            C[m * i + j] = A1[m * i + j] - A2[m * i + j];
    }

    return 0;
}
void vector_multiplying(double *a1, int n, double alpha)
{
    for(int i = 0; i < n; i++)
        a1[i] = a1[i] * alpha;
}

double r1(double *A, int n, double *x)
{
    double trace = 0, lambda = 0, norm = m_norm(A, n);

    for(int i = 0; i < n; i++)
    {
        trace += A[i * n + i];
        lambda += x[i];
    }

    return abs(trace - lambda) / norm;
}

double r2(double *A, int n, double *x)
{
    double A_lenght = 0, x_lenght = 0, norm = m_norm(A, n);

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            A_lenght += A[i * n + j] * A[j * n + i];

    for(int i = 0; i < n; i++)
        x_lenght += x[i] * x[i];

    return abs(sqrt(A_lenght) - sqrt(x_lenght)) / norm;
}