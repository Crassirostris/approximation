#include <iostream>
#include <cmath>
#include <functional>
#include <iomanip>
#include <string>

const double SEGMENT_LEFT_END = 0;
const double SEGMENT_RIGHT_END = 1;
const double LEFT_END_VALUE = 0;
const double INITIAL_STEP = 0.1;
const int ITERATIONS_COUNT = 4;

const double EPS = 1e-14;

double fun(double x, double y) {
    return -5 * y * cos(5 * x) + 25 * sin(10 * x);
}

double fun_der_y(double x, double y) {
    return -5 * cos(5 * x);
}

double true_y(double x) {
    return 10 * (sin(5 * x) + exp(-sin(5 * x)) - 1);
}

void print_hline() {
    std::cout << "------------------------------------------------------------" << std::endl;
}

void analyze_method(std::string method_name, double initial_step, int iterations_count, std::function<double(double, double, double)> calc_next) {
    std::cout << method_name << std::endl;
    print_hline();

    double step = initial_step;
    for (int i = 0; i < iterations_count; ++i, step /= 2) {
        std::cout << "Step: " << step << std::endl;
        double cur_value = LEFT_END_VALUE;
        for (double x = SEGMENT_LEFT_END; x < SEGMENT_RIGHT_END - EPS; x += step) {
            std::cout << std::fixed << std::setprecision(10) << cur_value << std::endl;
            cur_value = calc_next(step, x, cur_value);
        }
        std::cout << std::fixed << std::setprecision(10) << cur_value << std::endl;
    }

    print_hline();
}

double calc_true(double step, double x, double y) {
    return true_y(x + step);
}

double calc_next_euler_1(double step, double x, double y) {
    return y + step / 2 * (fun(x, y) + fun(x + step, y + step * fun(x, y)));
}

double calc_next_euler_2(double step, double x, double y) {
    double y0 = y;
    double y1 = y + step / 2 * (fun(x, y) + fun(x + step, y + step * fun(x, y)));

    std::function<double(double)> g = [step, y, x](double z) { return z - y - step * fun(x + step, z); };
    std::function<double(double)> g_der = [step, y, x](double z) { return 1 - step * fun_der_y(x + step, z); };

    while (abs(y1 - y0) >= EPS) {
        double y2 = y1 - g(y1) / g_der(y1);
        std::swap(y0, y1);
        std::swap(y1, y2);
    }

    return y1;
}

double calc_next_adams_1(double step, double x, double y) {
    double y0 = y;
    double y1 = y + step / 2 * (fun(x, y) + fun(x + step, y + step * fun(x, y)));

    std::function<double(double)> g = [step, y, x](double z) { return z - y - step / 2 * (fun(x + step, z) + fun(x, y)); };
    std::function<double(double)> g_der = [step, y, x](double z) { return 1 - step / 2 * fun_der_y(x + step, z); };

    while (abs(y1 - y0) >= EPS) {
        double y2 = y1 - g(y1) / g_der(y1);
        std::swap(y0, y1);
        std::swap(y1, y2);
    }

    return y1;
}

double calc_next_kr_3(double step, double x, double y) {
    double k1 = step * fun(x, y);
    double k2 = step * fun(x + step, y + k1);
    double k3 = step * fun(x + step / 2, y + k1 / 4 + k2 / 4);
    return y + k1 / 6 + k2 / 6 + k3 * 4 / 6;
}


int main() {
    analyze_method("True", INITIAL_STEP, ITERATIONS_COUNT, calc_true);
    analyze_method("Euler 1", INITIAL_STEP, ITERATIONS_COUNT, calc_next_euler_1);
    analyze_method("Euler 2", INITIAL_STEP, ITERATIONS_COUNT, calc_next_euler_2);
    analyze_method("KR3", INITIAL_STEP, ITERATIONS_COUNT, calc_next_kr_3);
}
