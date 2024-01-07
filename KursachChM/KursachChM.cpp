#include <iostream>
#include <fstream>

#define SIMPSON_SECTIONS 20

using namespace std;


// нужна функция фи от икс
double Phi(double x);
// нужна функция пси от икс. Для неё нужно посчитать интеграл от x-at до x+at конечно-разностным методом
double Psi(double x);
double integrate_psi(double lower, double upper);
double integrate_Simpson(double lower, double upper, double (*f)(double));
double u(double x, double t);
void write_u_to_file();

double const a = 1.0; // коэффициент а в формуле Даламбера
const int n_phi = 5;
double restr_phi[n_phi] = { 1, 2, 3, 5, 7 }; // restrictions (ограничения в кусочно заданной функции)
double (*phi[n_phi + 1])(double) = {    [](double x) { return 0.0; },
                                        [](double x) { return x - 1; },
                                        [](double x) { return -x + 3; },
                                        [](double x) { return 1.5*x - 4.5; },
                                        [](double x) { return -1.5*x + 10.5; },
                                        [](double x) { return 0.0; } };

const int n_psi = 2;
double restr_psi[n_psi] = { 1, 3 };
double (*psi[n_psi + 1])(double) = {    [](double x) { return 0.0; },
                                        [](double x) { return 1.0; },
                                        [](double x) { return 0.0; } };
int main()
{
    try
    {
        write_u_to_file();
    }
    catch (invalid_argument& e)
    {
        cout << e.what();
    }
}

void write_u_to_file()
{
    double x0 = -3;
    double xn = 11;
    int num_of_x = 560;
    double delta_x = (xn - x0) / num_of_x;
    double tn = 5;
    int num_of_t = 80;
    double delta_t = tn / num_of_t;
    ofstream out;
    out.open("points.txt");
    if (out.is_open())
    {
        cout << "File opened" << endl;
        out << num_of_x + 1 << endl; // кол-во иксов
        for (int i = 0; i <= num_of_x; i++)
            out << x0 + i * delta_x << " ";
        out << endl;

        out << num_of_t + 1 << endl; // кол-во t
        for (int i = 0; i <= num_of_t; i++)
            out << i * delta_t << " ";
        out << endl;
        
        for (double t = 0; t <= tn; t += delta_t)
        {
            for (double x = x0; x <= xn; x += delta_x)
            {
                out << u(x, t) << " ";
            }
            out << endl;
        }
    }
    else cout << "File not opened" << endl;
    
}

double u(double x, double t)
{
    if (t < 0) throw invalid_argument("Time variable t can't be negative");
    double xt1 = x - a * t;
    double xt2 = x + a * t;
    double phi_part = (Phi(xt1) + Phi(xt2)) / 2;
    double psi_part = integrate_psi(xt1, xt2) / (2 * a);
    return phi_part + psi_part;
}

// задаём кусочную функцию на бесконечности
// точек ограничений n => функций n+1
double Phi(double x)
{
    int i = 0;
    for (i = 0; i < n_phi; i++)
    {
        if (x < restr_phi[i])
            break;
    }
    return phi[i](x);
}

double Psi(double x)
{
    int i = 0;
    for (i = 0; i < n_psi; i++)
    {
        if (x < restr_psi[i])
            break;
    }
    return psi[i](x);
}

// проинтегрировать отдельно каждую функцию из ограничений
double integrate_psi(double lower, double upper)
{
    double res = 0;
    int cur_func = 0;
    for (cur_func = 0; cur_func < n_psi; cur_func++)
        if (restr_psi[cur_func] >= lower)
            break;

    while (cur_func < n_psi && restr_psi[cur_func] < upper)
    {
        res += integrate_Simpson(lower, restr_psi[cur_func], psi[cur_func]);
        lower = restr_psi[cur_func++];
    }
    res += integrate_Simpson(lower, upper, psi[cur_func]);
    

    return res;
}

double integrate_Simpson(double lower, double upper, double (*f)(double))
{
    double h = (upper - lower) / SIMPSON_SECTIONS;
    double res = 0;

    res += f(lower); // первый член формулы Симпсона
    for (int i = 1; i < SIMPSON_SECTIONS; i++)
    {
        double xi = lower + i * h;
        res += (2 * (i % 2) + 2) * f(xi); // коэффициент в формуле (4 при нечётном i и 2 при чётном)
    }
    res += f(upper);

    res *= h / 3.0;
    return res;
}