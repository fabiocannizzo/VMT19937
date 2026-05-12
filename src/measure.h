/**
 * @file AdaptiveBenchmark.hpp
 * @brief Statistical Adaptive Code Benchmarking with Outlier Rejection.
 *
 * THE CORE LOGIC:
 * This tool measures the execution time of a function 'f' repeatedly until the
 * sample average is guaranteed to be within a specific error margin (alpha%)
 * of the true population mean, with a given confidence level (p).
 *
 * STATISTICAL PRINCIPLES:
 * 1. Central Limit Theorem: As n increases, the distribution of the sample mean
 *    becomes normal, even if the timing data is skewed.
 * 2. Student's t-Distribution: Used because the population standard deviation
 *    (sigma) is unknown and must be estimated from the samples.
 * 3. One-Sided Outlier Rejection: Timing jitter (OS context switches, interrupts)
 *    only pushes execution time UP. We filter samples > Mean + N*Sigma to
 *    capture the "true" code performance under ideal conditions.
 * 4. Adaptive Stopping: The loop continues until:
 *    (t_critical * (StDev / sqrt(n))) / Mean <= Alpha
 */

#ifndef ADAPTIVE_BENCHMARK_HPP
#define ADAPTIVE_BENCHMARK_HPP

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <numeric>
#include <algorithm>

struct BenchmarkParams {
    double alpha = 0.01;         // Relative accuracy (0.01 = 1% error margin)
    double confidence = 0.95;    // Confidence level (0.95 = 95%)
    double outlier_sigma = 3.0;  // Threshold for one-sided outlier removal
    int min_reps = 30;           // Min samples to satisfy Central Limit Theorem
    int max_reps = 10000;        // Safety break to prevent infinite loops
    int warmup = 10;             // Initial runs to "heat" the CPU cache/branch predictor
};

struct BenchmarkResult {
    double mean_us;              // Final average time in microseconds
    double stdev_us;             // Standard deviation of filtered samples
    double rel_error;            // The actual relative error achieved
    size_t total_samples;        // Total iterations including outliers
    size_t filtered_samples;     // Iterations used in the final calculation
};

namespace StatsDetail {
    /**
     * Approximates the inverse CDF of the Student's t-distribution.
     * Uses the Abramowitz & Stegun approximation (26.7.10) to transform
     * a Normal Z-score into a T-statistic based on degrees of freedom.
     */
    inline double get_t_critical(double confidence, int df) {
        if (df < 1) return 12.706;

        // Probability for two-tailed test
        double p = 1.0 - ((1.0 - confidence) / 2.0);

        // Hastings approximation for Z-score
        double t_z = std::sqrt(-2.0 * std::log(1.0 - p));
        double z = t_z - (2.515517 + 0.802853 * t_z + 0.010328 * t_z * t_z) /
                       (1.0 + 1.432788 * t_z + 0.189269 * t_z * t_z + 0.001308 * t_z * t_z * t_z);

        // Transform Z to T
        double v = static_cast<double>(df);
        return z + (std::pow(z, 3) + z) / (4.0 * v) +
               (5.0 * std::pow(z, 5) + 16.0 * std::pow(z, 3) + 3.0 * z) / (96.0 * v * v);
    }
}

/**
 * @brief Executes the adaptive benchmark.
 * @param func The function to measure.
 * @param p The configuration parameters.
 */
template <typename Func>
inline BenchmarkResult run_adaptive_benchmark(Func&& func, BenchmarkParams p) {
    std::vector<double> raw_timings;
    raw_timings.reserve(p.min_reps * 2);

    // 1. Warm-up Phase: Discard results to stabilize hardware state
    for (int i = 0; i < p.warmup; ++i) {
        func();
    }

    // 2. Adaptive Measurement Loop
    while (raw_timings.size() < (size_t)p.max_reps) {
        // Collect a new sample
        auto start = std::chrono::high_resolution_clock::now();
        func();
        auto end = std::chrono::high_resolution_clock::now();
        raw_timings.push_back(std::chrono::duration<double, std::micro>(end - start).count());

        // We need at least min_reps to perform meaningful statistics
        if (raw_timings.size() < (size_t)p.min_reps) continue;

        // 3. Statistical Analysis
        double sum = std::accumulate(raw_timings.begin(), raw_timings.end(), 0.0);
        double raw_mean = sum / raw_timings.size();

        double sq_sum = std::inner_product(raw_timings.begin(), raw_timings.end(), raw_timings.begin(), 0.0);
        double raw_stdev = std::sqrt(std::max(0.0, (sq_sum / raw_timings.size()) - (raw_mean * raw_mean)));

        // 4. One-Sided Outlier Filtering (Remove high latency spikes)
        std::vector<double> filtered;
        double upper_bound = raw_mean + (p.outlier_sigma * raw_stdev);
        for (double val : raw_timings) {
            if (val <= upper_bound) filtered.push_back(val);
        }

        // 5. Check Stopping Condition on Filtered Data
        size_t n = filtered.size();
        if (n < 2) continue; // Safety check

        double f_sum = std::accumulate(filtered.begin(), filtered.end(), 0.0);
        double f_mean = f_sum / n;
        double f_sq_sum = std::inner_product(filtered.begin(), filtered.end(), filtered.begin(), 0.0);

        // Bessel's correction for sample variance: s^2 = (sum(x-avg)^2) / (n-1)
        double f_var = (f_sq_sum - (f_sum * f_sum) / n) / (n - 1);
        double f_stdev = std::sqrt(std::max(0.0, f_var));

        double t_stat = StatsDetail::get_t_critical(p.confidence, (int)n - 1);
        double margin_of_error = t_stat * (f_stdev / std::sqrt(static_cast<double>(n)));
        double current_rel_error = margin_of_error / f_mean;

        // Exit condition
        if (current_rel_error <= p.alpha && raw_timings.size() >= (size_t)p.min_reps) {
            return { f_mean, f_stdev, current_rel_error, raw_timings.size(), n };
        }
    }

    // If we hit max_reps, return the best estimate we have
    return { 0.0, 0.0, 0.0, raw_timings.size(), 0 };
}

#endif // ADAPTIVE_BENCHMARK_HPP

