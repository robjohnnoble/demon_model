#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

struct InputParameters {
    // capacity
    int log2_deme_carrying_capacity;

    // dispersal
    float init_migration_rate;
    int migration_rate_scales_with_K;

    // fission times (in numbers of fissions)
    int t0;
    int tL1;
    int tL2;
    int tL3;
    int tR1;
    int tR2;
    int tR3;

    // fitness
    float normal_birth_rate;
    float baseline_death_rate;
    float s_driver_birth;
    float s_driver_migration;
    float max_relative_birth_rate;
    float max_relative_migration_rate;

    // mutation
    float mu_driver_birth;
    float mu_driver_migration;

    // methylation
    float meth_rate;
    float demeth_rate;
    int fCpG_loci_per_cell;
    float manual_array;

    // seed
    int seed;

    // stopping conditions
    int max_time;
    int max_generations;
    int max_fissions;

    // initial conditions
    int init_pop;
    
    // output indicators
    int write_demes_file;
    int write_clones_file;
};

struct DerivedParameters {
    int K; // deme carrying capacity
    int fcpgs; // total number of fCpG sites per cell
    int density_dept_death_rate=100; // should be much larger than birth rates
    float dmax; // any genotypes that differ by at least this many mutations are considered totally unrelated (i.e. the distance between them is 1)
    int max_bintree_elements=0; // maximum number of bintree elements = max_demes * (1/SET_SIZE + 1/SET_SIZE^2 + ...) < max_demes/(SET_SIZE-1); allow an extra 10% due to rounding up
    int max_bintree_clone_elements_per_deme=0;
    bool use_clone_bintrees=false; // whether to use clone bintrees: faster for large K, but a waste of memory for small K
    int max_clones_per_deme;
    int max_genotypes;
    int max_driver_genotypes;
    int max_clones;
    int max_demes = 8;
};

struct EventCounter {
    int birth=0;
    int death=0;
    int mutation=0;
    int fission=0;
    int methylation=0;
    int demethylation=0;
};

#endif