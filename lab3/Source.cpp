#pragma comment(linker, "/STACK:67108864")

#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

const double Z = 1.202056903159594;
const double EPS = 1e-9;

#pragma region Arithmetics

double fix(double x, int d) {
    double p = pow(10, d);
    double m = x + EPS > 0 ? 1 : -1;
    return m * floor(m * x * p + 0.5 + EPS) / p;
}

double sum(double left, double right, int d) {
    return fix(left + right, d);
}

double sub(double left, double right, int d) {
    return fix(left - right, d);
}

double mul(double left, double right, int d) {
    return fix(left * right, d);
}

double div(double left, double right, int d) {
    return fix(left / right, d);
}

#pragma endregion

#pragma region Printing

const char *LatexDelimiter = " ";

void print(const vector<vector<double>> &m, int d) {
    for (int i = 0; i < m.size(); ++i) {
        for (int j = 0; j < m.size(); ++j)
            cout << fixed << setprecision(d) << m[i][j] << ' ';
        cout << endl;
    }
    cout << endl;
}

void print(const vector<double> &v, int d) {
    for (int i = 0; i < v.size(); ++i)
        cout << fixed << setprecision(d) << v[i] << endl;
    cout << endl;
}

void print(const pair<vector<vector<double>>, vector<vector<double>>> &lu, int d) {
    cout << "L:" << endl;
    print(lu.first, d);
    cout << "R:" << endl;
    print(lu.second, d);
}

void print_latex(const vector<double> &v, int d, const char *name) {
    cout << "$$" << name << "=\\left(\\begin{matrix}" << LatexDelimiter;
    for (int i = 0; i < v.size(); ++i)
        cout << fixed << setprecision(d) << v[i] << "\\\\" << LatexDelimiter;
    cout << "\\end{matrix}\\right)$$" << endl;
}

void print_latex(vector<vector<double>> &matrix, int d, const char *name) {
    cout << "$$" << name << "=\\left(\\begin{matrix}" << LatexDelimiter;
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            if (j > 0)
                cout << " & ";
            cout << fixed << setprecision(d) << matrix[i][j];
        }
        cout << "\\\\" << LatexDelimiter;
    }
    cout << "\\end{matrix}\\right)$$" << endl;
}

void print_system(const vector<vector<double>> &a, const vector<double> &b, const vector<int> &trans, int d) {
    int n = a.size();
    cout << "$$ \\left(\\begin{matrix}" << LatexDelimiter;
    for (int i = 0; i < n; ++i) {
        cout << '\t';
        for (int j = 0; j < n; ++j)
            cout << fixed << setprecision(d) << a[i][trans[j]] << (j == n - 1 ? " \\\\" : " & ");
        cout << LatexDelimiter;
    }
    cout << "\\end{matrix}\\right)" << LatexDelimiter << "\\left(\\begin{matrix}" << LatexDelimiter;
    for (int i = 0; i < n; ++i)
        cout << "x_" << trans[i] + 1 << " \\\\" << LatexDelimiter;
    cout << "\\end{matrix}\\right)" << LatexDelimiter << "=" << LatexDelimiter << "\\left(\\begin{matrix}" << LatexDelimiter;
    for (int i = 0; i < n; ++i)
        cout << fixed << setprecision(d) << b[i] << " \\\\" << LatexDelimiter;
    cout << "\\end{matrix}\\right) $$" << LatexDelimiter << endl;
}

#pragma endregion

#pragma region Support

vector<double> solve_triangular(vector<vector<double>> &coeffs, vector<double> &values, int d, bool up) {
    int n = values.size();
    vector<double> result(n, 0);

    for (int i = (up ? n - 1 : 0); up ? (i >= 0) : (i < n); up ? (--i) : (++i)) {
        double s = 0;
        for (int j = 0; j < n; ++j)
            s = sum(s, mul(result[j], coeffs[i][j], d), d);
        s = sub(values[i], s, d);
        result[i] = div(s, coeffs[i][i], d);
    }

    return result;
}

vector<double> prepare(const vector<double> &a, int d) {
    vector<double> res(a.size());
    for (int i = 0; i < a.size(); ++i)
        res[i] = fix(a[i], d);
    return res;
}

vector<vector<double>> prepare(const vector<vector<double>> &a, int d) {
    vector<vector<double>> res(a.size());
    for (int i = 0; i < a.size(); ++i)
        res[i] = prepare(a[i], d);
    return res;
}

#pragma endregion

#pragma region LU

pair<vector<vector<double>>, vector<vector<double>>> calc_lu(const vector<vector<double>> &a, int d) {
    vector<vector<double>> l(a.size(), vector<double>(a.size(), 0));
    vector<vector<double>> u(a.size(), vector<double>(a.size(), 0));
    for (int i = 0; i < a.size(); ++i)
        u[i][i] = 1;

    for (int k = 0; k < a.size(); ++k) {
        for (int i = k; i < a.size(); ++i) {
            l[i][k] = a[i][k];
            for (int m = 0; m <= k - 1; ++m)
                l[i][k] = sub(l[i][k], mul(l[i][m], u[m][k], d), d);
        }
        for (int j = k + 1; j < a.size(); ++j) {
            u[k][j] = a[k][j];
            for (int m = 0; m <= k - 1; ++m)
                u[k][j] = sub(u[k][j], mul(l[k][m], u[m][j], d), d);
            u[k][j] = div(u[k][j], l[k][k], d);
        }
    }

    return make_pair(l, u);
}

void solve_lu(const vector<vector<double>> &a, const vector<double> &b, int d) {
    cout << "LU with " << d << " digits" << endl << endl;

    auto a_fixed = prepare(a, d);
    auto b_fixed = prepare(b, d);

    auto lu = calc_lu(a_fixed, d);
    auto y = solve_triangular(lu.first, b_fixed, d, false);
    auto x = solve_triangular(lu.second, y, d, true);

    print_latex(lu.first, d, "L");
    print_latex(lu.second, d, "U");
    print_latex(y, d, "y");
    print_latex(x, d, "x");
    cout << endl;
    cout << "L" << endl;
    print(lu.first, d);
    cout << "R" << endl;
    print(lu.second, d);
    cout << "y" << endl;
    print(y, d);
    cout << "x" << endl;
    print(x, d);
}

#pragma endregion

#pragma region Gauss

void solve_gauss(const vector<vector<double>> &a, const vector<double> &b, int d) {
    cout << "Gauss with " << d << " digits" << endl << endl;
    auto a_fixed = prepare(a, d);
    auto b_fixed = prepare(b, d);

    int n = a.size();
    vector<int> trans(n);
    for (int i = 0; i < n; ++i)
        trans[i] = i;

    print_system(a_fixed, b_fixed, trans, d);
    for (int i = 0; i < n; ++i) {
        double max_val = 0, max_i = i;
        for (int j = i; j < n; ++j)
            if (max_val < abs(a_fixed[i][trans[j]])) {
                max_val = abs(a_fixed[i][trans[j]]);
                max_i = j;
            }
        swap(trans[i], trans[max_i]);
        if (i != max_i)
            print_system(a_fixed, b_fixed, trans, d);

        for (int j = i + 1; j < n; ++j) {
            auto c = div(a_fixed[j][trans[i]], a_fixed[i][trans[i]], d);
            b_fixed[j] = sub(b_fixed[j], mul(b_fixed[i], c, d), d);

            for (int k = i; k < n; ++k)
                a_fixed[j][trans[k]] = sub(a_fixed[j][trans[k]], mul(a_fixed[i][trans[k]], c, d), d);

            print_system(a_fixed, b_fixed, trans, d);
        }
    }

    vector<vector<double>> t(n, vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            t[i][j] = a_fixed[i][trans[j]];

    cout << endl;
    cout << "A:" << endl;
    print(t, d);
    cout << "b:" << endl;
    print(b_fixed, d);

    auto y = solve_triangular(t, b_fixed, d, true);
    cout << "y:" << endl;
    print(y, d);

    vector<double> solution(n);
    for (int i = 0; i < n; ++i)
        solution[trans[i]] = y[i];
    cout << "x:" << endl;
    print(solution, d);
}

#pragma endregion

int main()
{
    ios_base::sync_with_stdio(false);

    vector<vector<double>> a = {
        { 1.2345, 3.1415, 1 },
        { 2.3456, 5.9690, 0 },
        { 3.4567, 2.1828, 2.3 },
    };

    vector<double> b = { 7.6475, 14.2836, 8.1213 };

    for (int i = 1; i < 4; ++i)
        solve_lu(a, b, 2 * i);

    solve_gauss(a, b, 2);
    solve_gauss(a, b, 4);

    return 0;
}