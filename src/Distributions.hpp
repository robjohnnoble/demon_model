#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP

#include <random>

class RandomNumberGenerator {
public:
    RandomNumberGenerator() : rng(std::random_device{}()), dist(0.0, 1.0) {}

    // Allow manual setting of the seed
    void set_seed(unsigned int seed) {
        rng.seed(seed);
    }

    // U(0,1)
    double unif_ran() {
        return dist(rng);
    }

    // Poisson(lambda)
    int poisson(double lambda) {
        std::poisson_distribution<int> dist(lambda);
        return dist(rng);
    }

    // Exp(lambda)
    double exp(double lambda) {
        std::exponential_distribution<double> dist(lambda);
        return dist(rng);
    }

    // stochastic rounding
    int stochastic_round(double a) {
        float rnd = dist(rng);
        float fractional_part = a - static_cast<int>(a);
        if (rnd < fractional_part) {
            return static_cast<int>(a) + 1;
        } else {
            return static_cast<int>(a);
        }
    }

    // Hypergeometric(n, N, K)
    unsigned int hypergeometric(unsigned int n1, unsigned int n2, unsigned int t) {
        const unsigned int n = n1 + n2; // Total population size
        unsigned int a = n1;
        unsigned int b = n1 + n2;
        unsigned int k = 0; // Number of samples found so far

        if (t > n) t = n; // Sample size can't exceed population size

        if (t < n / 2) {
            for (unsigned int i = 0; i < t; ++i) {
                double u = unif_ran();
                if (b * u < a) {
                    ++k;
                    if (k == n1) return k;
                    --a;
                }
                --b;
            }
        } else {
            for (unsigned int i = 0; i < n - t; ++i) {
                double u = unif_ran();
                if (b * u < a) {
                    ++k;
                    if (k == n1) return n1 - k;
                    --a;
                }
                --b;
            }
            return n1 - k;
        }

        return k;
    }

    // get engine
    std::mt19937 get_engine() {
        return rng;
    }

private:
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist;
};

#endif // DISTRIBUTIONS_HPP