#include "func.h"

double m_norm(double *A, int n)
{
    double norm = 0, helper = 0;

    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            helper += fabs(A[i * n + j]);

        if(helper > norm)
            norm = helper;

        helper = 0;
    }

    return norm;
}
bool is_three_diagonal(Args *arg)
{
    int n = arg->n;
    double *a = arg->a;
    
    for (int j = 0; j < n - 2; j++)
        for(int i = j + 2; i < n; i++)
            if(fabs(a[i * n + j]) > 1e-100)
                return false;
            
    return true;
}
int improved_three_diagonal(Args *arg)
{
    double *a = arg->a;
    int n = arg->n;
    
    // Allocate temporary storage
    double *d = new double[n];   // Diagonal elements
    double *e = new double[n-1]; // Off-diagonal elements
    
    // Householder reduction to tridiagonal form
    for(int i = n-1; i > 0; i--) {
        // Extract the column vector
        double scale = 0.0;
        for(int k = 0; k < i; k++) {
            scale += fabs(a[k*n+i]);
        }
        
        if(scale == 0.0) {
            e[i-1] = a[(i-1)*n+i];
            continue;
        }
        
        double h = 0.0;
        for(int k = 0; k < i; k++) {
            a[k*n+i] /= scale;
            h += a[k*n+i] * a[k*n+i];
        }
        
        double f = a[(i-1)*n+i];
        double g = (f >= 0.0) ? -sqrt(h) : sqrt(h);
        e[i-1] = scale * g;
        h -= f * g;
        a[(i-1)*n+i] = f - g;
        f = 0.0;
        
        // Apply Householder transformation
        for(int j = 0; j < i; j++) {
            a[i*n+j] = a[j*n+i] / h;
            g = 0.0;
            
            for(int k = 0; k <= j; k++) {
                g += a[j*n+k] * a[k*n+i];
            }
            for(int k = j+1; k < i; k++) {
                g += a[k*n+j] * a[k*n+i];
            }
            
            e[j] = g / h;
            f += e[j] * a[j*n+i];
        }
        
        double hh = f / (h + h);
        for(int j = 0; j < i; j++) {
            f = a[j*n+i];
            g = e[j] - hh * f;
            e[j] = g;
            
            for(int k = 0; k <= j; k++) {
                a[j*n+k] -= f * e[k] + g * a[k*n+i];
            }
        }
    }
    
    // Copy tridiagonal elements back to the matrix
    for(int i = 0; i < n; i++) {
        d[i] = a[i*n+i];
        if(i < n-1) {
            a[i*n+i+1] = e[i];
            a[(i+1)*n+i] = e[i]; // Symmetry
        }
    }
    
    // Zero out the rest of the matrix
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(abs(i-j) > 1) {
                a[i*n+j] = 0.0;
            }
        }
        a[i*n+i] = d[i];
    }
    
    delete[] d;
    delete[] e;
    return 0;
}
int three_diagonal(Args *arg)
{
    double *a = arg->a, *xk = arg->xk;
    int n = arg->n;
    double norm;

    for(int k = 0; k < n - 2; k++)
    {
        //Считаем x(k)
        get_column(a, xk, n, k);

        norm = vector_norm(xk, n);
        xk[k + 1] -= norm;
        norm = vector_norm(xk, n);
        if(vector_division(xk, n, norm) == -1)
            continue;
        
        UAUt(arg);
    }

    return 0;
}
int sign(double l)
{
    if(l < 0)
        return -1.0;
    else
        return 1.0;
}
int sign_changes(Args* arg, double b)
{
    const int n = arg->n;
    double* matrix = arg->a; // Full matrix storage
    const double eps = arg->eps;
    const double norm = arg->norm;
    
    if(n <= 0) return 0;
    
    // Extract diagonal and off-diagonal elements from the full matrix
    double* diag = new double[n];       // Main diagonal
    double* offdiag = new double[n-1];  // Off-diagonal elements
    
    for(int i = 0; i < n; i++) {
        diag[i] = matrix[i*n+i];
        if(i < n-1) {
            offdiag[i] = matrix[i*n+i+1]; // Upper diagonal (equals lower diagonal for symmetric matrix)
        }
    }
    
    // Sturm sequence calculation (more numerically stable for symmetric tridiagonal matrices)
    int count = 0;
    double d_prev = diag[0] - b;
    
    // Count sign changes in the Sturm sequence
    if(d_prev < 0) count++;
    
    for(int i = 1; i < n; i++) {
        // Calculate next term in the sequence
        double d_curr;
        if(fabs(d_prev) < eps * norm) {
            // Avoid division by very small numbers
            d_curr = diag[i] - b - fabs(offdiag[i-1]*offdiag[i-1]) / (eps * norm);
        } else {
            d_curr = diag[i] - b - (offdiag[i-1]*offdiag[i-1]) / d_prev;
        }
        
        // Count sign changes
        if(d_curr * d_prev < 0) count++;
        
        d_prev = d_curr;
    }
    
    delete[] diag;
    delete[] offdiag;
    
    return count;
}
int solve(Args *arg)
{
    int n = arg->n;
    double *a = arg->a, *x = arg->x;
    double eps = arg->eps;

    arg->norm = m_norm(a, n);
    double norm = arg->norm;
    
    if (n == 1)
    {
        x[0] = a[0];
        return 0;
    }

    // Reduction to tridiagonal form
    double time = get_full_time();
    if(!is_three_diagonal(arg))
    {
        // Use the improved Householder reduction instead of rotation method
        if(improved_three_diagonal(arg) == -1)
        {
            cout << "Problem with tridiagonal reduction" << endl;
            return -1;
        }
        arg->three_diagonal_time = get_full_time() - time;
    }

    // Calculate eigenvalue bounds for the bisection method
    double b0 = 0;
    for(int i = 0; i < n; i++) {
        double row_sum = 0;
        for(int j = 0; j < n; j++) {
            row_sum += fabs(a[i*n+j]);
        }
        b0 = std::max(b0, row_sum);
    }
    
    double ai = -b0, bi = b0;

    // Find eigenvalues using bisection method
    time = get_full_time();
    for(int k = 0; k < n; k++)
    {
        bool is_end = false;
        while(bi - ai > eps * norm && !is_end)
        {
            double c = (ai + bi) / 2;

            if(bi - c < EPS * norm || c - ai < EPS * norm)
                is_end = true;

            int n_minus = sign_changes(arg, c);
            if(n_minus < k + 1)
                ai = c;
            else
                bi = c;

            arg->its++;
        }

        double lambda = (ai + bi) / 2;
        
        // Determine multiplicity
        int s1 = sign_changes(arg, bi), s2 = sign_changes(arg, ai);
        int multiplicity = s1 - s2;

        if(multiplicity == 0)
        {
            // Try a slightly wider interval
            double delta = eps * norm * 10;
            s1 = sign_changes(arg, bi + delta);
            s2 = sign_changes(arg, ai - delta);
            multiplicity = s1 - s2;
            if (multiplicity == 0) multiplicity = 1;
        }
        
        // Store eigenvalues
        for(int i = k; i < k + multiplicity && i < n; i++)
        {
            x[i] = lambda;
        }

        k += multiplicity - 1;
        bi = b0;
        ai = lambda + eps * norm;
    }
    arg->eigenvalues_time = get_full_time() - time;

    return 0;
}