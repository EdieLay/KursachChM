#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#define _USE_MATH_DEFINES
#include <math.h>

#define SIMPSON_SECTIONS 2 * 10

using namespace std;


// нужна функция фи от икс
double Phi(double x);
// нужна функция пси от икс. Для неё нужно посчитать интеграл от x-at до x+at конечно-разностным методом
double Psi(double x);
double integrate_psi(double lower, double upper);
double integrate_Simpson(double lower, double upper, double (*f)(double));
double u(double x, double t);
double u_half_finite(double x, double t, bool left_odd);
double u_finite(double x, double t, double l, bool left_odd, bool right_odd);
void write_u_to_file();

double const a = 1.0; // коэффициент а в формуле Даламбера
const int n_phi = 4;
double restr_phi[n_phi] = { 1, 3, 6, 9 }; // restrictions (ограничения в кусочно заданной функции)
double (*phi[n_phi + 1])(double) = { [](double x) { return 0.0; },
                                        [](double x) { return sin(M_PI / 2 * (x - 1)); }, // sin(M_PI / 2 * (x - 1))
                                        [](double x) { return 0.0; },
                                        [](double x) { return -2 * sin(M_PI / 3 * (x - 6)); },
                                        [](double x) { return 0.0; } };


const int n_psi = 2;
double restr_psi[n_psi] = { 1, 3 };
double (*psi[n_psi + 1])(double) = {    [](double x) { return 0.0; },
                                        [](double x) { return 0.0; },
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
    double x0 = 0;
    double xn = 5;
    int num_of_x = 100 * (int)xn;
    double delta_x = (xn - x0) / num_of_x;
    double tn = 6;
    int num_of_t = 10 * (int)tn;
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
                out << u_finite(x, t, xn, false, true) << " ";
            }
            out << endl;
        }
    }
    else cout << "File not opened" << endl;
    
}

// колебание струны на отрезке
// l - права граница отрезка (и его длина)
// left_odd - зафиксирована ли левая точка в нуле (т.е. слева нечетное отражение)
// right_odd - зафикирована ли правая точка в нуле (т.е. справа нечетное отражение)
double u_finite(double x, double t, double l, bool left_odd, bool right_odd)
{
    double xt1 = x - a * t;
    int part1 = 0; // часть, в которую попала точка xt1 (рисунок в папке "Курсач ЧМы")
    if (xt1 < 0)
    {
        while (xt1 < 0) // приводим точку к начальному отрезку
        {
            xt1 += l;
            part1++;
        }

        if (part1 % 2 == 1) xt1 = l - xt1; // в зависимости от кол-ва отражений точку надо будет ещё раз отразить
    }

    double xt2 = x + a * t;
    int part2 = 0;
    if (xt2 > l)
    {
        while (xt2 > l)
        {
            xt2 -= l;
            part2++;
        }

        if (part2 % 2 == 1) xt2 = l - xt2;
    }

    double phi1 = Phi(xt1);
    switch(part1 % 4) // отражения повторяются с циклом 4 (в 4 части в любом случае все точки принимают изначальный вид)
    { // отражаем относительно Ox в зависимости от части и четности отражений
    case 0:
        break;
    case 1:
        if (left_odd) phi1 *= -1;
        break;
    case 2:
        if (left_odd != right_odd) phi1 *= -1;
        break;
    case 3:
        if (right_odd) phi1 *= -1;
        break;
    default:
        break;
    }

    double phi2 = Phi(xt2);
    switch (part2 % 4)
    {
    case 0:
        break;
    case 1:
        if (right_odd) phi2 *= -1;
        break;
    case 2:
        if (left_odd != right_odd) phi2 *= -1;
        break;
    case 3:
        if (left_odd) phi2 *= -1;
        break;
    default:
        break;
    }

    // мы берем интеграл от отрицательной точки до нуля
    // проходимся отдельно по каждой части
    // в первую итерацию цикла мы берём интеграл от xt1 до l, где 0 <= xt1 < l
    // поэтому мы должны предусмотреть, в какой части это происходит, чтобы не взять больше или меньше, чем нужно
    // в последующие итерации мы всегда берём от 0 до l, 
    // поэтому присваивания xt1 и становление его то нижним, то верхним пределом - не имеет значения
    // просто это написано в общем виде
    double psi1 = 0;
    while (part1 > 0)
    {
        switch (part1 % 4)
        {
        case 0:
            psi1 += integrate_psi(xt1, l);
            xt1 = l;
            break;
        case 1:
            if (left_odd) psi1 -= integrate_psi(0, xt1);
            else psi1 += integrate_psi(0, xt1);
            xt1 = 0;
            break;
        case 2:
            if (left_odd != right_odd) psi1 -= integrate_psi(xt1, l);
            else psi1 += integrate_psi(xt1, l);
            xt1 = l;
            break;
        case 3:
            if (right_odd) psi1 -= integrate_psi(0, xt1);
            else psi1 += integrate_psi(0, xt1);
            xt1 = 0;
            break;
        default:
            break;
        }
        part1--;
    }

    double psi2 = 0;
    // то же самое, что и для xt1
    // только теперь мы берем интеграл от l до большего числа
    while (part2 > 0)
    {
        switch (part2 % 4)
        {
        case 0:
            psi2 += integrate_psi(0, xt2);
            xt2 = 0;
            break;
        case 1:
            if (right_odd) psi2 -= integrate_psi(xt2, l);
            else psi2 += integrate_psi(xt2, l);
            xt2 = l;
            break;
        case 2:
            if (left_odd != right_odd) psi2 -= integrate_psi(0, xt2);
            else psi2 += integrate_psi(0, xt2);
            xt2 = 0;
            break;
        case 3:
            if (left_odd) psi2 -= integrate_psi(xt2, l);
            else psi2 += integrate_psi(xt2, l);
            xt2 = l;
            break;
        default:
            break;
        }
        part2--;
    }

    // берем интеграл от начального отрезка, учитывая что xt1 мог стать нулём, если изначально он был меньше нуля
    // а xt2 мог стать l, если изначально он был больше l
    double psi_part = (psi1 + integrate_psi(xt1, xt2) + psi2) / (2 * a);

    double phi_part = (phi1 + phi2) / 2;
    return phi_part + psi_part;
}

// колебание струны на луче
double u_half_finite(double x, double t, bool left_odd)
{
    double xt1 = x - a * t;
    double xt2 = x + a * t;

    double phi1 = Phi(xt1);
    double phi2 = Phi(xt2);
    double psi_part = 0;
    // по сути то же самое, что и для отрезка, только теперь есть всего одно отражение
    if (xt1 < 0)
    {
        phi1 = Phi(-xt1);
        if (left_odd)
        {
            phi1 *= -1;
            psi_part -= integrate_psi(0, -xt1);
        }
        else
        {
            psi_part += integrate_psi(0, -xt1);
        }
        xt1 = 0;
    }
    psi_part += integrate_psi(xt1, xt2);
    psi_part /= (2 * a);

    double phi_part = (phi1 + phi2) / 2;

    return phi_part + psi_part;
}

// колебание струны на прямой
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
        if (x < restr_phi[i])
            break;
    return phi[i](x);
}

double Psi(double x)
{
    int i = 0;
    for (i = 0; i < n_psi; i++)
        if (x < restr_psi[i])
            break;
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
    res += f(upper); // последний член формулы Симпсона

    res *= h / 3.0;
    return res;
}