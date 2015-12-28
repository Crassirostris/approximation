#pragma comment(linker, "/STACK:67108864")

#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

double f(double x) {
    return cos(3 * x) - pow(x, 3);
}

double phi(double x, double q) {
    return x + f(x) * q;
}

pair<double, int> fpi(double x, double q) {
    int it = 0;
    while (abs((phi(x, q) - x) / 3 / q) >= 5e-6) {
        x = phi(x, q);
        ++it;
    }
    return make_pair(x, it);
}

pair<double, int> secant(double x, double y) {
    int it = 0;
    while (abs(f(x) / 3) >= 5e-6) {
        double t = x - f(x) * (x - y) / (f(x) - f(y));
        y = x;
        x = t;
        ++it;
    }
    return make_pair(x, it);
}

void print(pair<double, int> ans) {
    cout << fixed << setprecision(10) << ans.first << ' ' << ans.second << endl;
}

int main()
{
    ios_base::sync_with_stdio(false);

    double x = M_PI / 6, y = 5 * M_PI / 36;
    cout << "FPI" << endl;
    print(fpi(x, 0.25));
    print(fpi(x, 1 / 3.5));
    print(fpi(x, 0.01));
    print(fpi(6, 0.3));
    cout << "Secant" << endl;
    print(secant(x, y));
    print(secant(1, 0));
    print(secant(6, 0));
    print(secant(6, 0.5));

    return 0;
}