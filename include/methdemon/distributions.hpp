#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP

#include <random>

class RandomNumberGenerator {
public:
    static RandomNumberGenerator& getInstance();
    void setSeed(unsigned int seed);
    double unitUnifDist();
    int poissonDist(double lambda);
    double expDist(double lambda);
    int stochasticRound(double a);
    unsigned int hypergeometricDist(unsigned int n1, unsigned int n2, unsigned int t);
    std::mt19937 getEngine();
private:
    RandomNumberGenerator();
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist;
};

#endif // DISTRIBUTIONS_HPP