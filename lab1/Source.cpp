#pragma comment(linker, "/STACK:67108864")

#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

const double Z = 1.202056903159594;

double fix(double x, int d) {
    ostringstream oss;
    oss << fixed << setprecision(14) << x;
    string s = oss.str();
    int pos = 0;
    while (s[pos] == '0' || s[pos] == '.')
        ++pos;
    for (int st = 0; pos < s.length(); ++pos)
        if (isdigit(s[pos]) && (st++ > d))
            s[pos] = '0';
    return stod(s);
}

double sum(double left, double right, int d) {
    return fix(left + right, d);
}

double mul(double left, double right, int d) {
    return fix(left * right, d);
}

double div(double left, double right, int d) {
    return fix(left / right, d);
}

double elem1(double n, int d) {
    double f1 = mul(mul(0.2, n, d), n, d);
    double f2 = mul(1.3, n, d);
    return div(1, sum(sum(f1, f2, d), 1.7, d), d);
}

double elem2(double n, int d) {
    double num = sum(mul(13, n, d), 17, d);
    double sq = mul(n, n, d);
    double f1 = mul(sq, sq, d);
    double f2 = mul(sq, mul(6.5, n, d), d);
    double f3 = mul(sq, 8.5, d);
    double den = sum(sum(f1, f2, d), f3, d);
    return div(num, den, d);
}

double elem3(double n, int d) {
    double num = sum(mul(135, n, d), 221, d);
    double sq = mul(n, n, d);
    double tr = mul(n, mul(n, n, d), d);
    double f1 = mul(tr, sq, d);
    double f2 = mul(tr, mul(6.5, n, d), d);
    double f3 = mul(tr, 8.5, d);
    double den = sum(sum(f1, f2, d), f3, d);
    return div(num, den, d);
}

double elem4(double n, int d) {
    double num = sum(mul(1313, n, d), 2295, d);
    double sq = mul(n, n, d);
    double qd = mul(sq, sq, d);
    double f1 = mul(qd, sq, d);
    double f2 = mul(qd, mul(6.5, n, d), d);
    double f3 = mul(qd, 8.5, d);
    double den = sum(sum(f1, f2, d), f3, d);
    return div(num, den, d);
}

double convert1(double res, int d) {
    return res;
}

double convert2(double res, int d) {
    return sum(5 * M_PI * M_PI / 6, mul(-2.5, res, d), d);
}

double convert3(double res, int d) {
    res = sum(65 * Z / 2, mul(-1.25, res, d), d);
    return sum(5 * M_PI * M_PI / 6, -res, d);
}

double convert4(double res, int d) {
    res = sum(15 * pow(M_PI, 4) / 8, mul(-0.625, res, d), d);
    res = sum(65 * Z / 2, -res, d);
    return sum(5 * M_PI * M_PI / 6, -res, d);
}

double exec(int m, int d, double elem(double, int), double convert(double, int)) {
    double res = 0;
    for (int i = m; i > 0; --i) {
        res = sum(res, elem(static_cast<double>(i), d), d);
    }
    return convert(res, d);
}

typedef double(*func)(double, int);

int main()
{
    ios_base::sync_with_stdio(false);

    int sizes[] = { 20, 100, 1000 };
    func elem_funcs[] = { elem1, elem2, elem3, elem4 };
    func convert_funcs[] = { convert1, convert2, convert3, convert4 };

    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 4; ++i) {
            cout << fixed << setprecision(10) << exec(sizes[j], 15, elem_funcs[i], convert_funcs[i]) << '\t';
        }
        cout << endl;
    }

    cout << fixed << setprecision(10) << exec(2e5, 10, elem_funcs[0], convert_funcs[0]) << endl;
    cout << fixed << setprecision(10) << exec(807, 8, elem_funcs[1], convert_funcs[1]) << endl;
    cout << fixed << setprecision(10) << exec(132, 7, elem_funcs[2], convert_funcs[2]) << endl;
    cout << fixed << setprecision(10) << exec(54, 7, elem_funcs[3], convert_funcs[3]) << endl;

    return 0;
}