#include <iostream>
#include <immintrin.h>
#include <random>
#include <cmath>
#include <chrono>

float definite_integral_128(float a, float b, float (*function)(float x), int N, std::mt19937 gen) {
    int size = 4;
    float res = 0;
    float rank = (b - a) / N;
    std::uniform_real_distribution<float> dist(a, b);
    __m128 intrinsic_res = _mm_setzero_ps();
    __m128 intrinsic_rank = _mm_set1_ps(rank);
    for (int i = 0; i < N; i += size) {
        __m128 group_x = _mm_set_ps(function(dist(gen)), function(dist(gen)),
        function(dist(gen)), function(dist(gen)));
        intrinsic_res = _mm_add_ps(intrinsic_res, _mm_mul_ps(group_x, intrinsic_rank));
    }
    float tmp_sum[size];
    _mm_store_ps(&tmp_sum[0], intrinsic_res);
    for (int i = 0; i < size; ++i)
        res += tmp_sum[i];
    return res;
}

float definite_integral_256(float a, float b, float (*function)(float x), int N, std::mt19937 gen) {
    int size = 8;
    float res = 0;
    float rank = (b - a) / N;
    std::uniform_real_distribution<float> dist(a, b);
    __m256 intrinsic_res = _mm256_setzero_ps();
    __m256 intrinsic_rank = _mm256_set1_ps(rank);
    for (int i = 0; i < N; i += size) {
        __m256 group_x = _mm256_set_ps(function(dist(gen)), function(dist(gen)),
        function(dist(gen)), function(dist(gen)), function(dist(gen)), function(dist(gen)),
        function(dist(gen)), function(dist(gen)));

        intrinsic_res = _mm256_add_ps(intrinsic_res, _mm256_mul_ps(group_x, intrinsic_rank));
    }
    float tmp_sum[size];
    _mm256_store_ps(&tmp_sum[0], intrinsic_res);
    for (int i = 0; i < size; ++i)
        res += tmp_sum[i];
    return res;
}

float definite_integral(float a, float b, float (*function)(float x), int N, std::mt19937 gen) {
    float res = 0;
    float rank = (b - a) / N;
    std::uniform_real_distribution<float> dist(a, b);
    for (int i = 0; i < N; ++i) {
        float x = dist(gen);
        res += function(x) * rank;
    }
    return res;
}

float calc_pi(int N, std::mt19937 gen) {
    int in_circle = 0;
    std::uniform_real_distribution<float> dist(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
        float x = dist(gen), y = dist(gen);
        if (std::pow(x, 2) + std::pow(y, 2) <= 1)
            in_circle += 1;
    }
    return (float)(4.0 * in_circle) / N;
}

int main() {
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    std::random_device rd;
    std::mt19937 gen(rd());
    int N = 1000000;
    float pi = calc_pi(N, gen);
    float a1 = 2.312, b1 = 11.7; auto function1 { [](float x) -> float { return (std::pow(2, x) * std::cos(x)); } };
    float a2 = pi / 2, b2 = pi; auto function2 { [](float x) -> float { return std::sin(x) / (std::pow(std::cos(x), 2) + 1); } };
    float a3 = -0.5, b3 = 0.5; auto function3 { [](float x) -> float { return std::sin(x) / (std::acos(2 * x)); } };

    float integral_1 = definite_integral_128(a1, b1, function1, N, gen); // choose function what you want to execute
    float integral_2 = definite_integral_128(pi / 2, pi, function2, N, gen);
    float integral_3 = definite_integral_128(a3, b3, function2, N, gen);
    std::cout << "Pi number: " << pi << "\n";
    std::cout << "definite integral of function (2 ^ x) * cos(x) from " << a1 << " to " << b1 << " is " << integral_1 << "\n";
    std::cout << "definite integral of function sin(x) / (cos(x)^2 + 1) from " << a2 << " to " << b2 << " is " << integral_2 << "\n";
    std::cout << "definite integral of function arccos(2x) from " << a3 << " to " << b3 << " is " << integral_3 << "\n";
    
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Execution time: " << elapsed_seconds.count() << " seconds" << std::endl;
    return 0;
}