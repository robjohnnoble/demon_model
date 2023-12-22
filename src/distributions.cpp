#include "distributions.hpp"

// constructor
RandomNumberGenerator::RandomNumberGenerator() : rng(std::random_device()()), dist(0.0, 1.0) {}

// get instance
RandomNumberGenerator& RandomNumberGenerator::getInstance() {
    static RandomNumberGenerator instance;
    return instance;
}

// manual seed set
void RandomNumberGenerator::setSeed(unsigned int seed) {
    rng.seed(seed);
}

// U(0,1)
double RandomNumberGenerator::unitUnifDist() {
    return dist(rng);
}

// Poisson(lambda)
int RandomNumberGenerator::poissonDist(double lambda) {
    std::poisson_distribution<int> dist(lambda);
    return dist(rng);
}

// Exp(lambda)
double RandomNumberGenerator::expDist(double lambda) {
    std::exponential_distribution<double> dist(lambda);
    return dist(rng);
}

// stochastic rounding
int RandomNumberGenerator::stochasticRound(double a) {
    float rnd = dist(rng);
    float fractional_part = a - static_cast<int>(a);
    if (rnd < fractional_part) {
        return static_cast<int>(a) + 1;
    } else {
        return static_cast<int>(a);
    }
}

// Hypergeometric(n, N, K)
unsigned int RandomNumberGenerator::hypergeometricDist(unsigned int n1, unsigned int n2, unsigned int t) {
    const unsigned int n = n1 + n2; // Total population size
    unsigned int a = n1;
    unsigned int b = n1 + n2;
    unsigned int k = 0; // Number of samples found so far

    if (t > n) t = n; // Sample size can't exceed population size

    if (t < n / 2) {
        for (unsigned int i = 0; i < t; ++i) {
            double u = unitUnifDist();
            if (b * u < a) {
                ++k;
                if (k == n1) return k;
                --a;
            }
            --b;
        }
    } else {
        for (unsigned int i = 0; i < n - t; ++i) {
            double u = unitUnifDist();
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
std::mt19937 RandomNumberGenerator::getEngine() {
    return rng;
}