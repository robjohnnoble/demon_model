#include "Objects.hpp"

// Constructor
Deme::Deme(int K, std::string side, int identity, int population, int fissions,
    float deathRate, float sumBirthRates, float sumMigrationRates)
    : K(K), side(side), identity(identity), population(population), fissions(fissions),
        death_rate(deathRate),
        sum_birth_rates(sumBirthRates), sum_migration_rates(sumMigrationRates) {}

void Deme::calculate_sum_of_rates() {
    sum_rates = sum_birth_rates + sum_migration_rates + population * death_rate;
}

// rates handling
void Deme::set_death_rate(const InputParameters& params) {
    if (population <= K) {
        death_rate = params.baseline_death_rate;
    }
    else {
        death_rate = params.baseline_death_rate + 1000;
    }
}

void Deme::increment(int increment, const InputParameters& params, std::string origin) {
    population += increment;
    set_death_rate(params);
    // check population sum
    int num_clones_in_deme = clones_list.size();
    if (population != num_clones_in_deme) {
        std::cout << "Population does not equal number of clones in deme " << identity 
        << " population: " << population << " num_clones_in_deme: " << num_clones_in_deme << std::endl
        << " increment: " << increment 
        << " origin: " << origin << std::endl;
        exit(1);
    }
}

