#include <iostream>
#include <random>
#include <cmath>

double definite_integral(double a, double b, double (*function)(double x), int N) {
    double res = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    double rank = (b - a) / N;
    std::uniform_real_distribution<double> dist(a, b);
    for (int i = 0; i < N; ++i) {
        double x = dist(gen);
        res += function(x) * rank;
    }
    return res;
}

double calc_pi(int N) {
    int in_circle = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
        double x = dist(gen), y = dist(gen);
        if (std::pow(x, 2) + std::pow(y, 2) <= 1)
            in_circle += 1;
    }
    return (double)(4.0 * in_circle) / N;
}

int main() {
    double a = 2.312, b = 11.7;
    int N = 1000000;
    auto function { [](double x) -> double { return (std::pow(2, x) * std::cos(x)); } };
    double integral = definite_integral(a, b, function, N);
    double pi = calc_pi(N);
    std::cout << "definite integral of function (2 ^ x) * cos(x) from " << a << " to " << b << " is " << integral << "\n";
    std::cout << "Число пи: " << pi << "\n";
    return 0;
}