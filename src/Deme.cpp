#include "Objects.hpp"

// Constructor
Deme::Deme(int K, std::string side, int identity, int population, int fissions,
    float deathRate, float migrationModifier, float sumBirthRates, float sumMigrationRates)
    : K(K), side(side), identity(identity), population(population), fissions(fissions),
        death_rate(deathRate), migration_modifier(migrationModifier),
        sum_birth_rates(sumBirthRates), sum_migration_rates(sumMigrationRates) {}

void Deme::calculate_sum_of_rates() {
    sum_rates = sum_birth_rates + sum_migration_rates + death_rate;
}

// Methods
void Deme::calculate_average_array(const DerivedParameters& d_params, std::vector<Clone>& clones) {
    std::vector<float> res(d_params.fcpgs, 0);

    for (int i = 0; i < clones_list.size(); i++) {
        for (int j = 0; i < d_params.fcpgs; j++) {
            res[j] += clones_list[i]->meth_array[j] * clones_list[i]->population;
        }
    }

    for (int i = 0; i < d_params.fcpgs; i++) {
        avg_meth_array[i] = res[i] / static_cast<float>(population);
    }
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

