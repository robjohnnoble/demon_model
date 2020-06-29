#ifndef CONFIGURATION_STRUCTURES_H_INCLUDED
#define CONFIGURATION_STRUCTURES_H_INCLUDED

//---------------conf file---------------------------begin
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
//---------------conf file-----------------------------end



struct SpatialStructure{
	int log2_deme_carrying_capacity;  // log2 of deme carrying capacity
															 // examples: 0, 3, 10.
};

struct Dispersal{
	float migration_type; 	// Type of cell migration: 
													// 0=new cell always tries to migrate after birth event;
													// 1=migration rate independent of birth events;
													// 2=deme fission when population size reaches carrying capacity;
													// 3=deme fission rate independent of population size.
	float init_migration_rate; 	// Initial cell migration rate or deme fission rate
												  // examples: 0.1 or 1.
    int migration_edge_only;  	// If 1 then cells cannot migrate into demes that are already full
    											  // (i.e. demes that have cancer cell population >= K);
												  // examples: 0 or 1.
    int migration_rate_scales_with_K; 	// If 1 then migraton rates are divided by sqrt(K)
												  // examples: 0 or 1.

};

struct MutationRates{
	float mu_driver_birth;         		 // Driver mutation rate per cancer cell division;
														 // examples: 0 or 0.00001.
	float mu_passenger;              // Passenger mutation rate per cancer cell division;
														 // examples: 0 or 0.00001.
	float mu_driver_migration;              // Migration rate mutation rate per cancer cell division;
														 // examples: 0 or 0.00001.
	int passenger_pop_threshold;              // population size at which passenger mutations stop occurring;
														 // examples: -1 (no threshold) or 100000
};

struct FitnessEffects{
	float normal_birth_rate; // Normal cell birth rate, relative to cancer cell birth rate,
													 // except that -ve values imply total absence of normal cells;
													 // examples: 0.9 or -1.
	float baseline_death_rate; // Baseline death rate regardless of deme population size
													 // examples: 0 or 0.9
	float s_passenger;			         // Negative effect on birth rate per passenger mutation;
													 // examples: 1E-3.
	float s_driver_birth;              // Positive effect on birth rate per driver mutation;
													 // examples: 1E-1.
	float s_driver_migration;							 // Positive effect on migration rate per migration rate mutation;
													 // examples: 1E-1.
	float max_relative_birth_rate;              // maximum birth rate, relative to initial birth rate
													 // examples: 10.
	float max_relative_migration_rate;              // maximum migration rate, relative to initial migration rate
													 // examples: 10.
};



struct NonBiologicalParameters{
	int seed; 						// Seed for random number generator.
	int init_pop;       	// Initial cancer cell population;
												// examples: 10.
	int max_time;    			// Max elapsed time before program stops, in seconds;
												// examples: 86400.
	int max_pop;          // Max cancer cell population before program stops;
												// examples: 20000.
	int max_generations;          // Max generations before program stops;
												// examples: 500.
	int matrix_max;       // Max number of genotypes before program stops
												// (larger value => more memory demanded);
												// examples: 5000.
	int write_grid; // Whether to make plots of grid states
																// (can significantly increase running time when grid size is large);
																// examples: 1.
	int write_clones_file; // Whether to write contents of all clones
																// (can be a very large file but is essential for analysis);
																// examples: 1.
	int write_demes_file; // Whether to write contents of all demes
																// (can be a very large file and not essential);
																// examples: 0.
	int record_matrix; // Whether to record distance matrix for all genotypes (not just driver genotypes)
								// automatically set to 0 if mu_passenger = 0
																// (can be a very large file and not essential);
																// examples: 0.
	int write_phylo; // Whether to write phylogeny for all genotypes (not just driver genotypes)
																// (can be a very large file);
																// examples: 
																		// 0 (never write);
																		// 0.5 (write data for whole tumour but not for samples);
																		// 1 (write both).
	int calculate_total_diversity; // Whether to calculate diversity across all genotypes
																// examples: 0 (don't calculate)
																// examples: 1 (calculate)
    int biopsy_size_per_sample; // max number of cells per biopsy sample
																// examples: 100.
};

struct Parameters{
public:
	SpatialStructure 				 spatial_structure;
	Dispersal 			 			   dispersal;
	MutationRates    				 mutation_rates;
	FitnessEffects   				 fitness_effects;
	NonBiologicalParameters  non_biological_parameters;

	void GetFromFile(std::string & file_name);
	void SetToFile(std::string & file_name);
	void PrintToConsole();
	boost::property_tree::ptree getPTree();
	void GetFromPTree(boost::property_tree::ptree pt);

	private:
	boost::property_tree::ptree p_tree;
};

#endif // CONFIGURATION_STRUCTURES_H_INCLUDED
