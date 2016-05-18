#include <iostream>
#include <cmath>
#include <functional>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

const double SEGMENT_LEFT_END = 0;
const double SEGMENT_RIGHT_END = 5;

const double LEFT_END_VALUE = -6.4;
const double RIGHT_END_VALUE = exp(5) + 2 * exp(-5) - 0.9 * 71;

const double INITIAL_STEP = 1;
const int ITERATIONS_COUNT = 5;

const double COMPRESSING_FACTOR = 0.8;

const double EPS = 1e-9;

struct Point {
    Point(double x, double y) : x(x), y(y) { }

    Point operator+(const Point &other) const {
        return Point(x + other.x, y + other.y);
    }

    Point operator-(const Point &other) const {
        return Point(x - other.x, y - other.y);
    }

    Point operator*(double v) const {
        return Point(x * v, y * v);
    }

    Point operator/(double v) const {
        return Point(x / v, y / v);
    }

    double x;
    double y;
};

double given_poly(double x) {
    return 0.9 * (x * x * x - x * x + 2);
}

Point fun(double x, Point y) {
    return Point(y.y, y.x + given_poly(x));
}

Point calc_next_y_rk3(double x, Point y, double h) {
    auto k1 = fun(x, y) * h;
    auto k2 = fun(x + h, y + k1) * h;
    auto k3 = fun(x + h / 2, y + k1 / 4.0 + k2 / 4.0) * h;
    return y + k1 / 6 + k2 / 6 + k3 * 4 / 6;
}

std::vector<double> get_solution(double param, double step) {
    Point current = Point(param, LEFT_END_VALUE);

    std::vector<double> result = { param };
    for (double x = SEGMENT_LEFT_END; x < SEGMENT_RIGHT_END - EPS; x += step) {
        current = calc_next_y_rk3(x, current, step);
        result.push_back(current.x);
    }

    return result;
}

double find_cauchy_last(double param, double step) {
    Point current = Point(param, LEFT_END_VALUE);

    for (double x = SEGMENT_LEFT_END; x < SEGMENT_RIGHT_END - EPS; x += step) {
        current = calc_next_y_rk3(x, current, step);
    }

    return current.y;
}

double find_parameter(double step) {
    double l = -100, r = 100;
    while (abs(l - r) >= EPS) {
        double kauchy_last = find_cauchy_last((l + r) / 2, step);

        if (abs(kauchy_last - RIGHT_END_VALUE) < EPS) {
            return (l + r) / 2;
        }

        if ((kauchy_last - RIGHT_END_VALUE) * (find_cauchy_last(l, step) - RIGHT_END_VALUE) > 0) {
            l = l * COMPRESSING_FACTOR + r * (1 - COMPRESSING_FACTOR);
        } else {
            r = l * (1 - COMPRESSING_FACTOR) + r * COMPRESSING_FACTOR;
        }
    }
    return (l + r) / 2;
}

void print_solutions(std::vector<std::vector<double>> solutions) {
    int max_len = 0;
    for (auto &sol : solutions) {
        max_len = std::max(max_len, (int) sol.size());
    }

    for (int j = 0; j < max_len; ++j) {
        for (auto &solution : solutions) {
            if (j < solution.size()) {
                std::cout << solution[j];
            }
            std::cout << '\t';
        }
        std::cout << std::endl;
    }
}

std::vector<double> get_ticks(double step) {
    std::vector<double> result;
    for (double x = SEGMENT_LEFT_END; x < SEGMENT_RIGHT_END + EPS; x += step) {
        result.push_back(x);
    }
    return result;
}

double true_fun(double x) {
    return 1.00018161 * exp(x) + 2.00018161 * exp(-x) - 0.9 * (x * x * x - x * x + 6 * x);
}

double get_error(std::vector<double> &solution, std::vector<double> &ticks) {
    double error = 0;
    for (int i = 0; i < solution.size(); ++i) {
        error = std::max(error, abs(solution[i] - true_fun(ticks[i])));
    }
    return error;
}

std::vector<double> get_true_solution(const std::vector<double> &ticks) {
    std::vector<double> result;
    for (auto x : ticks) {
        result.push_back(true_fun(x));
    }
    return result;
}

std::vector<std::vector<double>> perform_shooting() {
    double step = INITIAL_STEP;
    std::vector<std::vector<double>> solutions;
    for (int i = 0; i < ITERATIONS_COUNT; ++i, step /= 2) {
        std::cout << std::fixed << std::setprecision(10) << "Shooting with step " << step << std::endl;
        double parameter = find_parameter(step);
        std::cout << std::fixed << std::setprecision(10) << "Found parameter: " << parameter << std::endl;
        std::vector<double> solution = get_solution(parameter, step);
        auto ticks = get_ticks(step);
        std::cout << std::fixed << std::setprecision(10) << "Error: " << get_error(solution, ticks) << std::endl;
        solutions.push_back(solution);
    }
    return solutions;
}

std::vector<double> solve_triagonal(double step) {
    int n = static_cast<int>((SEGMENT_RIGHT_END - SEGMENT_LEFT_END + EPS) / step) + 1;
    std::vector<double> second_diagonal(n);
    std::vector<double> values(n);

    for (int i = 1; i < n - 1; ++i) {
        second_diagonal[i] = -(2 + step * step);
        values[i] = given_poly(SEGMENT_LEFT_END + step * i) * step * step;
    }
    second_diagonal[0] = -1;
    values[0] = step * LEFT_END_VALUE;
    second_diagonal[n - 1] = -1;
    values[n - 1] = -step * RIGHT_END_VALUE;

    // Fake points method
    //second_diagonal[0] = -(2 + step * step) / 2;
    //values[0] = step * LEFT_END_VALUE + given_poly(SEGMENT_LEFT_END) * step * step / 2;
    //second_diagonal[n - 1] = -(2 + step * step) / 2;
    //values[n - 1] = - step * RIGHT_END_VALUE - given_poly(SEGMENT_RIGHT_END) * step * step / 2;

    for (int i = 0; i < n - 1; ++i) {
        double c = -1 / second_diagonal[i];
        second_diagonal[i + 1] += 1 * c;
        values[i + 1] += values[i] * c;
    }

    for (int i = n - 1; i > 0; --i) {
        double c = -1 / second_diagonal[i];
        values[i - 1] += values[i] * c;
    }

    std::vector<double> solution;
    for (int i = 0; i < n; ++i) {
        solution.push_back(values[i] / second_diagonal[i]);
    }

    std::cout << std::fixed << std::setprecision(10) << "3diag error: " << get_error(solution, get_ticks(step)) << std::endl;
    return solution;
}

std::vector<std::vector<double>> get_tridiagonal_soltuions() {
    double step = INITIAL_STEP;
    std::vector<std::vector<double>> solutions;

    for (int i = 0; i < ITERATIONS_COUNT; ++i, step /= 2) {
        std::vector<double> solution = solve_triagonal(step);
        solutions.push_back(solution);
    }

    return solutions;
}

std::vector<std::vector<double>> make_solutions(std::vector<std::vector<double>> shooting_solutions, std::vector<std::vector<double>> tridiagonal_solutions) {
    double step = INITIAL_STEP;
    std::vector<std::vector<double>> solutions;
    for (int i = 0; i < ITERATIONS_COUNT; ++i, step /= 2) {
        solutions.push_back(get_ticks(step));
        solutions.push_back(get_true_solution(get_ticks(step)));
        solutions.push_back(shooting_solutions[i]);
        solutions.push_back(tridiagonal_solutions[i]);
    }
    return solutions;
}

int main() {
    FILE *f;
    freopen_s(&f, "output.txt", "wt", stdout);

    auto shooting_solutions = perform_shooting();
    auto tridiagonal_solutions = get_tridiagonal_soltuions();

    auto solutions = make_solutions(shooting_solutions, tridiagonal_solutions);

    print_solutions(solutions);
}
