#include "Objects.hpp"

// Constructor
Deme::Deme(int K, std::string side, int identity, int population, int fissions,
    float deathRate, float sumBirthRates, float sumMigrationRates)
    : K(K), side(side), identity(identity), population(population), fissions(fissions),
        death_rate(deathRate),
        sum_birth_rates(sumBirthRates), sum_migration_rates(sumMigrationRates) {}

void Deme::calculate_sum_of_rates() {
    sum_rates = sum_birth_rates + sum_migration_rates + death_rate;
}

// rates handling
void Deme::set_death_rate(const InputParameters& params) {
    if (population <= K) {
        death_rate = params.baseline_death_rate;
    }
    else {
        death_rate = params.baseline_death_rate + 100;
    }
}

void Deme::increment(int increment, const InputParameters& params) {
    population += increment;
    set_death_rate(params);
}

