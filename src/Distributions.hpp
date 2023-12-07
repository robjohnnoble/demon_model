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

private:
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist;
};

#endif // DISTRIBUTIONS_HPP