#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
#include <string>

const double SEGMENT_LEFT_END = 0.3;
const double SEGMENT_RIGHT_END = 1.1;
const double INITIAL_STEP = 0.1;
const double ACTUAL_INT_VALUE = -0.447545;
const int ITERATIONS_COUNT = 5;

const double EPS = 1e-9;

double fun(double x) {
    return log(atan(x));
}

double fun_der(double x) {
    return 1 / (1 + x * x) / atan(x);
}

double fun_corrected(double x) {
    return 0.4 * fun(0.4 * x + 0.7);
}

double runge(double int_part, double int_double_part, int order) {
    return exp2(order) * (int_double_part - int_part) / (exp2(order) - 1);
}

void analyze_method(std::string method_name, std::function<double(double)> calc, int method_order, int iterations, double initial_step, double actual) {
    std::vector<double> values(iterations + 1);
    double current_step = initial_step;
    for (int i = 0; i <= iterations; ++i, current_step /= 2) {
        values[i] = calc(current_step);
    }

    std::vector<double> runge_errors(iterations);
    for (int i = 0; i < iterations; ++i) {
        runge_errors[i] = runge(values[i], values[i + 1], method_order);
    }

    current_step = initial_step;
    for (int i = 0; i < iterations; ++i, current_step /= 2) {
        std::cout << std::fixed << std::setw(12) << std::setprecision(10)
            << method_name << '\t'
            << current_step << '\t'
            << values[i] << '\t'
            << abs(runge_errors[i]) << '\t'
            << abs(actual - values[i]) << std::endl;
    }
}

double gauss(double step) {
    return fun_corrected(-sqrt(0.6)) * 5 / 9 + fun_corrected(0) * 8 / 9 + fun_corrected(sqrt(0.6)) * 5 / 9;
}

double euler(double step) {
    double sum = 0;
    for (double point = SEGMENT_LEFT_END; point < SEGMENT_RIGHT_END + EPS; point += step) {
        sum += fun(point);
    }
    return step * (2 * sum - fun(SEGMENT_LEFT_END) - fun(SEGMENT_RIGHT_END)) / 2 + step * step / 12 * (fun_der(SEGMENT_LEFT_END) - fun_der(SEGMENT_RIGHT_END));
}

double center_rect(double step) {
    double sum = 0;
    for (double point = SEGMENT_LEFT_END; point < SEGMENT_RIGHT_END + EPS; point += step) {
        sum += fun(point + step / 2);
    }
    return sum * step;
}

double gauss_comp(double step) {
    double sum = 0;
    for (double point = SEGMENT_LEFT_END; point < SEGMENT_RIGHT_END + EPS; point += step) {
        double center = point + step / 2;
        double d = step / 2 * sqrt(0.6);
        sum += fun(center - d) * 5 / 9 + fun(center) * 8 / 9 + fun(center + d) * 5 / 9;
    }
    return sum * step / 2;
}

int main() {
    analyze_method("C Rectangles", center_rect, 2, ITERATIONS_COUNT, INITIAL_STEP, ACTUAL_INT_VALUE);
    analyze_method("Euler", euler, 4, ITERATIONS_COUNT, INITIAL_STEP, ACTUAL_INT_VALUE);
    analyze_method("Gauss", gauss, 2, 2, INITIAL_STEP, ACTUAL_INT_VALUE);
    analyze_method("Gauss Comp", gauss_comp, 2, ITERATIONS_COUNT, INITIAL_STEP, ACTUAL_INT_VALUE);
}