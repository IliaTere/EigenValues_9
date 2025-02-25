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
    const double* t = arg->a; // Трехдиагональная матрица в упакованном виде
    const double eps = arg->eps;
    const double norm = arg->norm;

    if(n <= 0) return 0;
    
    const double abs_eps = std::max(eps * norm, 1e-300);
    int count = 0;
    int current_sign = 1;
    
    // Первый элемент главной диагонали
    double diag = t[0] - b;
    
    // Регуляризация с сохранением знака
    if(fabs(diag) < abs_eps) {
        diag = (diag >= 0) ? abs_eps : -abs_eps;
    }
    current_sign = (diag > 0) ? 1 : -1;

    // Обработка трехдиагональной структуры
    for(int i = 1; i < n; ++i) {
        const int offset = 2*i - 1;
        const double lower = t[offset];     // Элемент нижней диагонали
        const double upper = t[offset + 1]; // Элемент верхней диагонали
        
        // Стабильное вычисление элемента разложения
        diag = t[2*i] - b - (lower * upper) / diag;
        
        // Адаптивная регуляризация
        const double adaptive_eps = std::max(abs_eps, 1e-15 * fabs(t[2*i]));
        if(fabs(diag) < adaptive_eps) {
            diag = (current_sign > 0) ? adaptive_eps : -adaptive_eps;
        }
        
        // Проверка изменения знака
        const int new_sign = (diag > 0) ? 1 : -1;
        if(new_sign != current_sign) {
            count++;
            current_sign = new_sign;
        }
    }
    
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

    //Приведение к трехдиагональному виду
    double time = get_full_time();
    if(!is_three_diagonal(arg))
    {
        if(three_diagonal(arg) == -1)
        {
            cout << "Pr with three_diagonal" << endl;
            return -1;
        }
        for (int j = 0; j < n - 2; j++)
            for(int i = j + 2; i < n; i++)
            {
                a[i * n + j] = 0;
                a[j * n + i] = 0;
            }
        arg->three_diagonal_time = get_full_time() - time;
    }

    //Вычисляем b0
    double b0 = arg->norm;
    double ai = -b0, bi = b0, c = 0;

    //Вычисляем собственные значения
    time = get_full_time();
    for(int k = 0; k < n; k++)
    {
        bool is_end = false;
        while(bi - ai > eps * norm and !is_end)
        {
            c = (ai + bi) / 2;

            if(bi - c < EPS * norm or c - ai < EPS * norm)
                is_end = true;

            int n_minus = sign_changes(arg, c);
            if(n_minus < k + 1)
                ai = c;
            else
                bi = c;

            arg->its++;
        }

        double lambda = (ai + bi) / 2;
        int s1 = sign_changes(arg, bi), s2 = sign_changes(arg, ai);
        int multiplicity = s1 - s2;

        if(multiplicity == 0)
        {
            s1 = sign_changes(arg, bi + eps * norm);
            s2 = sign_changes(arg, ai - eps * norm);
            multiplicity = s1 - s2;
            if (multiplicity == 0) multiplicity++;
        }
        
        for(int i = k; i < k + multiplicity; i++)
        {
           if(i > n - 1) break;
            x[i] = lambda;
        }

        k += multiplicity - 1;
        bi = b0;
        ai = lambda + eps * norm;
    }
    arg->eigenvalues_time = get_full_time() - time;

    return 0;
}