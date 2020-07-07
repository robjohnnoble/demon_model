#include "demon.h"

// two-dim arrays:
int **clone_ints;
double **bintree_clone_doubles;
int **bintree_clone_ints;
int **genotype_ints, **driver_genotype_ints;
float **genotype_floats, **driver_genotype_floats;
int **deme_ints;
float **deme_floats;
double **deme_doubles;
int **clones_list_in_deme, **clones_pops_in_deme;
double **clones_rates_in_deme;
int **matrix, **driver_matrix; // distances between genotypes, driver genotypes
int **grid; // maps from grid coordinates to deme number
int **genotype_relatives;
int **freq_table;
float **depth_diversity;
float **depth_diversity_bigsample;

// one-dim arrays:
double *bintree_deme_doubles;
int *empty_cols, *empty_driver_cols; // list of columns of distance matrix that can be replaced (because counts[] = 0)
int *position_in_edge_list, *demes_at_edge, *sides_at_edge, *genotype_edge_pop, *at_edge;
float *buff_array;
int *allele_count;
float *within_deme_diversity, *within_deme_driver_diversity;
int *max_layer_needed_array;
int *min_deme_bintree_index_by_layer, *min_clone_bintree_index_by_layer;
float *buff_small_float;
double *buff_small_double;
int *buff_array_int;

// parameters read from config file:
int log2_deme_carrying_capacity; // log2 of deme carrying capacity
int migration_type; // type of cell migration / deme fission
float init_migration_rate; // initial migration / deme fission rate
int migration_edge_only; // whether migration / deme fission occurs only at the tumour edge
int migration_rate_scales_with_K; // whether migration / deme fission rate should be divided by sqrt(K)
float mu_driver_birth; // rate of driver mutations affecting birth rate
float mu_driver_migration; // rate of driver mutations affecting migration rate
float mu_passenger; // passenger mutation rate
int passenger_pop_threshold; // population size at which passenger mutations stop occurring
float normal_birth_rate; // birth rate of normal cells, relative to initial cancer cell birth rate
float baseline_death_rate; // baseline death rate regardless of deme population size
float s_driver_birth; // mean positive effect on birth rate per driver mutation
float s_driver_migration; // mean positive effect on migration rate per driver mutation
float s_passenger; // mean cost to birth rate per passenger mutation
float max_relative_birth_rate; // maximum birth rate, relative to initial birth rate
float max_relative_migration_rate; // maximum migration rate, relative to initial birth rate
int init_pop; // initial population size
int matrix_max; // maximum number of matrix columns (for assigning memory)
long max_time; // max execution time (seconds) before the program halts
int max_generations; // max generations before the program halts
int seed; // seed for random number generator
int write_grid; // whether to plot grids using gnuplot
int write_clones_file; // whether to write clones file
int write_demes_file; // whether to write demes file
int record_matrix; // whether to record distance matrix for all genotypes (not just driver genotypes)
int write_phylo; // whether to write phylo file for all genotypes (not just drivers)
int calculate_total_diversity; // whether to calculate diversity for all genotypes (not just drivers)
int biopsy_size_per_sample; // max number of cells per biopsy sample

// derived parameters:
int K; // deme carrying capacity
int max_pop, max_demes, dim_grid, max_clones_per_deme, max_genotypes, max_driver_genotypes, max_clones, max_distinct_allele_freqs, length_of_max_layer_needed_array;
int well_mixed; // whether the model is being run within a single deme (depends on other parameter values)
int filled_grid; // whether the grid is initially fully occupied by cancer cells
int density_dept_death_rate = 100; // should be much larger than birth rates
float dmax; // any genotypes that differ by at least this many mutations are considered totally unrelated (i.e. the distance between them is 1)
int max_bintree_deme_elements; // maximum number of bintree elements = max_demes * (1/SET_SIZE + 1/SET_SIZE^2 + ...) < max_demes/(SET_SIZE-1); allow an extra 10% due to rounding up
int max_bintree_clone_elements_per_deme;
int use_clone_bintrees; // whether to use clone bintrees: faster for large K, but a waste of memory (and possibly slower) for small K
int max_gens = 1000000; // needed only for creating image file names; should be larger than the generations count when the program halts (but will be increased if needed)

// output files:
FILE *output_parameters, *output_pops, *output_diversities, *output_demes, *output_genotype_counts, *output_driver_genotype_counts, *output_phylo, *output_driver_phylo, *output_clones;
FILE *output_matrix, *output_driver_matrix, *output_popgrid, *output_passengersgrid, *output_normalcellsgrid, *output_deathrates_grid, *output_birthratesgrid, *output_migrationratesgrid, *output_driversgrid;
FILE *sample_size_log, *output_allele_counts, *output_driver_allele_counts, *output_genotype_properties, *output_driver_genotype_properties;
FILE *gp; // gnuplot pipe

int iterations; // number of iterations (attempted cell events)
int total_events[7]; // count events throughout simulation
long int powers_list[21]; // powers of SET_SIZE

// error tracking:
FILE *error_log;
int exit_code;

//////////////

int main(int argc, char *argv[])
{
	char *input_and_output_path = get_input_path(argc, argv);

	// config_file_with_path contains the name of the config file with its path:
	char *config_file_with_path = (char *) malloc(strlen(input_and_output_path) + strlen(argv[2]) + 1);
	strcpy(config_file_with_path, input_and_output_path);
	strcat(config_file_with_path, argv[2]);

	boost::property_tree::ptree pt; // create pt, which will hold all the parameter values
	boost::property_tree::info_parser::read_info(config_file_with_path, pt); // read from the config file to pt

	read_parameters(pt); // read parameter values from file

	assign_memory();

	run_sim(input_and_output_path, config_file_with_path); // run the simulation

	free(input_and_output_path);
	free(config_file_with_path);
	free_memory();

	return exit_code;
}

/////////////////////// set up:

// get input path:
char* get_input_path(int argc, char *argv[])
{
	if(argc < 3) { // if config file name is missing
		printf("Missing argument providing path to the config file and/or the config file name.\n");
		exit(0);
	}

	unsigned len = strlen(argv[1]); // argv[1] is the path of the config file
	char *input_and_output_path = (char *) malloc(len + 1 + 1); // +1 for extra char if needed, +1 for trailing zero

	if(argv[1][(strlen(argv[1])-1)] != '/') { // if path doesn't end in '/'
		strcpy(input_and_output_path, argv[1]);
		input_and_output_path[len] = '/'; // add '/' to the end of the path
		input_and_output_path[len + 1] = '\0'; // add the trailing zero
	}
	else strcpy(input_and_output_path, argv[1]);

	return input_and_output_path;
}

// read parameter values from config file:
void read_parameters(boost::property_tree::ptree pt)
{
	int i;
	float buff;
	float max_driver_mu;
	float predicted_clones_per_deme;

	// spatial parameters:
	log2_deme_carrying_capacity = pt.get <int> ("parameters.spatial_structure.log2_deme_carrying_capacity"); // log2 of deme carrying capacity
	
	// dispersal parameters:
	migration_type = pt.get <int> ("parameters.dispersal.migration_type"); // type of cell migration / deme fission
	init_migration_rate = pt.get <float> ("parameters.dispersal.init_migration_rate"); // initial migration / deme fission rate
	migration_edge_only = pt.get <int> ("parameters.dispersal.migration_edge_only"); // whether migration / deme fission occurs only at the tumour edge
	migration_rate_scales_with_K = pt.get <int> ("parameters.dispersal.migration_rate_scales_with_K"); // whether migration / deme fission rate should be divided by sqrt(K)

	// fitness effects:
	normal_birth_rate = pt.get <float> ("parameters.fitness_effects.normal_birth_rate"); // birth rate of normal cells, relative to initial cancer cell birth rate
	baseline_death_rate = pt.get <float> ("parameters.fitness_effects.baseline_death_rate"); // baseline death rate regardless of deme population size
	s_driver_birth = pt.get <float> ("parameters.fitness_effects.s_driver_birth"); // mean positive effect on birth rate per driver mutation
	s_driver_migration = pt.get <float> ("parameters.fitness_effects.s_driver_migration"); // mean positive effect on migration rate per driver mutation
	s_passenger = pt.get <float> ("parameters.fitness_effects.s_passenger"); // mean cost to birth rate per passenger mutation
	max_relative_birth_rate = pt.get <float> ("parameters.fitness_effects.max_relative_birth_rate"); // maximum birth rate, relative to initial birth rate
	max_relative_migration_rate = pt.get <float> ("parameters.fitness_effects.max_relative_migration_rate"); // maximum migration rate, relative to initial birth rate

	// mutation rates:
	mu_driver_birth = pt.get <float> ("parameters.mutation_rates.mu_driver_birth"); // rate of driver mutations affecting birth rate
	mu_driver_migration = pt.get <float> ("parameters.mutation_rates.mu_driver_migration"); // rate of driver mutations affecting migration rate
	mu_passenger = pt.get <float> ("parameters.mutation_rates.mu_passenger"); // passenger mutation rate
	passenger_pop_threshold = pt.get <int> ("parameters.mutation_rates.passenger_pop_threshold"); // population size at which passenger mutations stop occurring

	// seed:
	seed = pt.get <int> ("parameters.non_biological_parameters.seed"); // seed for random number generator

	// other parameters:
	record_matrix = pt.get <int> ("parameters.non_biological_parameters.record_matrix"); // whether to record distance matrix for all genotypes (not just driver genotypes) 
	biopsy_size_per_sample = pt.get <int> ("parameters.non_biological_parameters.biopsy_size_per_sample"); // max number of cells per biopsy sample
	init_pop = pt.get <int> ("parameters.non_biological_parameters.init_pop"); // initial number of cancer cells
	max_pop = pt.get <int> ("parameters.non_biological_parameters.max_pop"); // number of cancer cells at which the program halts
	max_time = pt.get <int> ("parameters.non_biological_parameters.max_time"); // max execution time (seconds) before the program halts
	max_generations = pt.get <int> ("parameters.non_biological_parameters.max_generations"); // max generations before the program halts
	write_grid = pt.get <int> ("parameters.non_biological_parameters.write_grid"); // whether to create images of the grid state
	write_demes_file = pt.get <int> ("parameters.non_biological_parameters.write_demes_file"); // whether to write deme states to file
	write_clones_file = pt.get <int> ("parameters.non_biological_parameters.write_clones_file"); // whether to write clone states to file
	write_phylo = pt.get <int> ("parameters.non_biological_parameters.write_phylo"); // whether to write phylogeny for all genotypes (not just driver genotypes) 
	calculate_total_diversity = pt.get <int> ("parameters.non_biological_parameters.calculate_total_diversity"); // whether to calculate diversity across all genotypes (not just driver genotypes) 
	matrix_max = pt.get <int> ("parameters.non_biological_parameters.matrix_max"); // maximum number of matrix columns (for assigning memory)

	K = pow(2, log2_deme_carrying_capacity); // deme carrying capacity

	if(mu_driver_birth < 0) mu_driver_birth = pow(10, mu_driver_birth);
	if(mu_driver_migration < 0) mu_driver_migration = pow(10, mu_driver_migration);
	if(mu_passenger < 0) mu_passenger = pow(10, mu_passenger);

	// preset migration rates:
	if(init_migration_rate < 0) {
		if(migration_type < 2) {
			if(fabs(normal_birth_rate - 0.9) < 0.01) {
				if(migration_edge_only == 0) init_migration_rate = set_init_migration_rate(K, init_migration_rate, -0.1593, -0.2868, 0.4646);
				else init_migration_rate = set_init_migration_rate(K, init_migration_rate, -0.2041, -0.14, 0.5761);
			}
			else init_migration_rate = -init_migration_rate * MIN(1, pow(K, -0.93));
		}
		else {
			if(migration_edge_only == 0) init_migration_rate = set_init_migration_rate(K, init_migration_rate, -0.15, -0.8928, -2.3445);
			else init_migration_rate = set_init_migration_rate(K, init_migration_rate, -0.0431, -1.7951, 0.0726);
		}
	}
	
	if(migration_rate_scales_with_K) init_migration_rate /= sqrt(K); // sqrt(K) corresponds to deme diameter

	if(init_migration_rate == 0 && init_pop <= K && (mu_driver_migration == 0 || s_driver_migration == 0)) { // cells cannot escape initial deme
		well_mixed = 1;
		dim_grid = 1;
		if(max_pop < 0) max_pop = 2 * K;
	}
	else if(max_pop < 0) { // grid is initially fully occupied
		filled_grid = 1;
		well_mixed = 0;
		dim_grid = ceil(MAX(3, sqrt((float)init_pop / K))); // width of the square grid (measured in demes)
		max_pop = ceil(dim_grid * dim_grid * 1.2 * K);
	}
	else { // expanding population
		filled_grid = 0;
		well_mixed = 0;
		buff = ((float)log2_deme_carrying_capacity / 20 + 1) * 4 * sqrt(max_pop / (PI*K));
		dim_grid = (int)MAX(3, buff); // width of the square grid (measured in demes)
	}
	if(max_pop <= init_pop) {
		printf("Error: max_pop (%d) <= init_pop (%d)\n", max_pop, init_pop);
		exit(1);
	}
	max_demes = dim_grid * dim_grid; // maximum number of demes (for assigning memory)

	if(passenger_pop_threshold < 0) passenger_pop_threshold = 2 * MAX(max_pop, init_pop);
	
	max_clones_per_deme = ceil(MIN(MAX(K + 5, 1.2 * K), max_pop)); // a bit larger than K

	max_driver_mu = MAX(mu_driver_migration, mu_driver_birth);
	
	max_genotypes = MIN(400 * MAX(mu_passenger * passenger_pop_threshold, max_driver_mu * max_pop), 1e8);
	max_genotypes = MAX(max_genotypes, 10);
	
	if(matrix_max <= 0) max_driver_genotypes = MIN(400 * max_driver_mu * max_pop, 1e8);
	else max_driver_genotypes = matrix_max;
	max_driver_genotypes = MAX(max_driver_genotypes, 10);

	predicted_clones_per_deme = MAX(MIN(ceil(K * MAX(max_driver_mu, mu_passenger) * 400), K), 1);

	max_clones = MIN(MAX(max_genotypes, max_demes * predicted_clones_per_deme), max_pop);

	for(i = 0; i < 21; i++) powers_list[i] = power(SET_SIZE, i); // lookup table for integer powers of SET_SIZE

	max_distinct_allele_freqs = max_pop; // max number of distinct allele frequencies

	dmax = 10; // any genotypes that differ by at least this many mutations are considered totally unrelated (i.e. the distance between them is 1)
	
	if(K > 16) use_clone_bintrees = 1; // use clone bintrees only if K is sufficiently large for it to be worthwhile
	else use_clone_bintrees = 0;

	max_bintree_deme_elements = ceil(1.1 * (float)max_demes / (SET_SIZE - 1));
	max_bintree_clone_elements_per_deme = ceil(1.1 * (float)max_clones_per_deme / (SET_SIZE - 1));
	// maximum number of bintree elements = arraysize * (1/SET_SIZE + 1/SET_SIZE^2 + ...) < arraysize/(SET_SIZE-1);
	// allow an extra 10% due to rounding up

	if(dim_grid < 2 || dim_grid > max_grid_for_output) write_grid = 0; // don't write any graphical output

	if(mu_passenger == 0) record_matrix = 0; // don't record matrix for all mutations if there are no passenger mutations

	length_of_max_layer_needed_array = MAX(max_demes, max_clones_per_deme) + 1;

	printf("Done reading parameters\n");
}

float set_init_migration_rate(int K, float init_migration_rate, float A, float B, float C)
{
	float mig_rate = -init_migration_rate * pow(10, A * pow(log10(K), 2) + B * log10(K) + C);
	
	if(migration_type == 3) mig_rate *= K;

	return mig_rate;
}

// initialise all variables:
void initialise(int *num_cells, int *num_clones, int *num_demes, int *num_matrix_cols, int *num_empty_cols, int init_driver_birth_mutations, int init_driver_mig_mutations, 
	int init_passengers, int *num_empty_driver_cols, int *num_driver_matrix_cols, int *next_driver_genotype_id, int *next_genotype_id, long *idum, int *num_extinct_genotypes, 
	int *num_empty_demes, int *num_extinct_driver_genotypes)
{
	int i, j, x, y;
	int index;
	int init_diameter;
	int remaining_pop;
	int max_layer_needed;

	for(i = 0; i <= MAX(max_demes, max_clones_per_deme); i++) max_layer_needed_array[i] = get_max_layer_needed(i);
	for(i = 0; i <= max_layer_needed_array[max_demes]; i++) min_deme_bintree_index_by_layer[i] = get_bintree_index(i, 0, max_demes);
	for(i = 0; i <= max_layer_needed_array[max_clones_per_deme]; i++) min_clone_bintree_index_by_layer[i] = get_bintree_index(i, 0, max_clones_per_deme);

	// counts:
	*num_demes = ceil((float)init_pop / K);
	*num_empty_demes = 0;
	*num_clones = *num_demes;
	*num_matrix_cols = 1;
	*num_driver_matrix_cols = 1;
	*num_empty_cols = 0;
	*num_empty_driver_cols = 0;
	*num_extinct_genotypes = 0;
	*num_extinct_driver_genotypes = 0;
	*num_cells = init_pop;

	init_diameter = (int)ceil(sqrt(*num_demes)); // initial diameter (measured in demes)
	printf("dim_grid %d, max_demes %d, num_demes = %d, init_diameter = %d\n", dim_grid, max_demes, *num_demes, init_diameter);
	printf("clone bintrees ");
	if(!use_clone_bintrees) printf("not ");
	printf("in use; max_pop %d; max_genotypes %d\n", max_pop, max_genotypes);

	// warn if dim_grid is too large for gnuplot to create images:
	if(write_grid && dim_grid > max_grid_for_output) {
		printf("Warning: grid size too large for graphical output; grids won't be plotted.\n");
		fprintf(error_log, "Warning: grid size too large for graphical output; grids won't be plotted.\n");
	}

	// genotypes:
	genotype_ints[POPULATION][0] = init_pop;
	genotype_ints[IDENTITY][0] = 0; // unique ID
	genotype_ints[DRIVER_IDENTITY][0] = 0; // unique ID
	genotype_ints[PARENT][0] = 0; // initial genotype has no parent
	genotype_ints[NUM_DRIVER_MUTATIONS][0] = init_driver_birth_mutations; // number of drivers
	genotype_ints[NUM_MIGRATION_MUTATIONS][0] = init_driver_mig_mutations; // number of drivers
	genotype_ints[IMMORTAL][0] = 1; // don't overwrite this genotype
	genotype_ints[NUM_PASSENGER_MUTATIONS][0] = init_passengers; // number of passengers
	*next_genotype_id = 1;
	// set genotype birth and miration rates:
	genotype_floats[BIRTH_RATE][0] = set_birth_rate(init_driver_birth_mutations, init_passengers, 1.0, idum);
	genotype_floats[MIGRATION_RATE][0] = set_migration_rate(init_driver_mig_mutations, init_migration_rate, idum);
	genotype_floats[ORIGIN_TIME][0] = 0; // generation at which genotype originated

	// driver genotypes:
	driver_genotype_ints[POPULATION][0] = init_pop;
	driver_genotype_ints[IDENTITY][0] = 0;
	driver_genotype_ints[DRIVER_IDENTITY][0] = 0;
	driver_genotype_ints[PARENT][0] = 0;
	driver_genotype_ints[NUM_DRIVER_MUTATIONS][0] = init_driver_birth_mutations; // number of drivers
	driver_genotype_ints[NUM_MIGRATION_MUTATIONS][0] = init_driver_mig_mutations; // number of drivers
	driver_genotype_ints[IMMORTAL][0] = 1; // don't overwrite this genotype
	driver_genotype_ints[NUM_PASSENGER_MUTATIONS][0] = init_passengers; // number of passengers
	*next_driver_genotype_id = 1;
	// set driver genotype birth and miration rates:
	driver_genotype_floats[BIRTH_RATE][0] = genotype_floats[BIRTH_RATE][0];
	driver_genotype_floats[MIGRATION_RATE][0] = genotype_floats[MIGRATION_RATE][0];
	driver_genotype_floats[ORIGIN_TIME][0] = 0; // generation at which genotype originated

	// demes:
	remaining_pop = init_pop;
	for(i=0; i < init_diameter; i++) for(j=0; j < init_diameter; j++) if(remaining_pop > 0) {
		index = i * init_diameter + j;
		x = dim_grid / 2 - init_diameter / 2 + i;
		y = dim_grid / 2 - init_diameter / 2 + j;
		deme_ints[POPULATION][index] = MIN(K, remaining_pop); // initially all cells belong to demes near centre of grid
		deme_ints[XCOORD][index] = x;
		deme_ints[YCOORD][index] = y;
		if(normal_birth_rate >= 0) deme_ints[NORMAL_CELLS][index] = K - deme_ints[POPULATION][index]; // if normal_birth_rate >= 0 then fill up the deme with normal cells
		else deme_ints[NORMAL_CELLS][index] = 0; // otherwise don't add any normal cells
		deme_ints[NUM_CLONES_IN_DEME][index] = 1;
		deme_floats[DEATH_RATE][index] = set_death_rate(index, *num_cells);
		deme_floats[MIGRATION_MODIFIER][index] = set_migration_modifier(deme_ints[NORMAL_CELLS][index], deme_ints[POPULATION][index]);
		deme_doubles[SUM_BIRTH_RATES][index] = (double)genotype_floats[BIRTH_RATE][0] * (double)deme_ints[POPULATION][index];	
		deme_doubles[SUM_MIGRATION_RATES][index] = (double)genotype_floats[MIGRATION_RATE][0] * (double)deme_ints[POPULATION][index];
		deme_doubles[SUM_RATES][index] = (double)deme_ints[POPULATION][index] * (double)deme_floats[DEATH_RATE][index] + deme_doubles[SUM_BIRTH_RATES][index];
		deme_doubles[SUM_RATES][index] += (double)deme_ints[NORMAL_CELLS][index] * ((double)normal_birth_rate + (double)deme_floats[DEATH_RATE][index]);
		if(migration_type == 1 || migration_type == 3) deme_doubles[SUM_RATES][index] += deme_doubles[SUM_MIGRATION_RATES][index];
		remaining_pop -= deme_ints[POPULATION][index]; // number of cancer cells still to be assigned to demes
	}

	// set deme bintree sums for layer 0:
	set_bintree_sums_layer0(bintree_deme_doubles, deme_doubles[SUM_RATES], *num_demes);

	// set deme bintree sums for subsequent layers:
	max_layer_needed = max_layer_needed_array[*num_demes];
	set_bintree_sums_subsequent_layers(bintree_deme_doubles, *num_demes, max_layer_needed, max_demes);

	// clones:
	for(index=0; index < *num_clones; index++) {
		clone_ints[POPULATION][index] = deme_ints[POPULATION][index];
		clone_ints[DEME][index] = index;
		clone_ints[GENOTYPE][index] = 0;
		clone_ints[DRIVER_GENOTYPE][index] = 0;
		clone_ints[INDEX_IN_DEME][index] = 0;
		clones_list_in_deme[index][0] = index;
		if(use_clone_bintrees) set_clone_in_deme(index, 0, 0, clone_ints[POPULATION][index]);
	}

	// grid:
	for(i=0; i<dim_grid; i++) for(j=0; j<dim_grid; j++) grid[i][j] = EMPTY; // means that grid squares have not yet been occupied
	for(i=0; i < init_diameter; i++) for(j=0; j < init_diameter; j++) {
		index = i * init_diameter + j;
		x = dim_grid / 2 - init_diameter / 2 + i;
		y = dim_grid / 2 - init_diameter / 2 + j;
		if(index < *num_demes) grid[x][y] = index;
	}

	// clone bintrees for each deme:
	if(use_clone_bintrees) reset_clone_bintree_sums(*num_demes, *num_clones);

	// matrix:
	if(record_matrix) for(i=0; i<matrix_max; i++) for(j=0; j<matrix_max; j++) matrix[i][j] = -99; // means that these entries of the matrix are not yet used
	matrix[0][0] = 0; // matrix comprises one entry, which is zero (the distance from the first genotype to itself)

	// driver matrix:
	if(matrix_max > 0) for(i=0; i<matrix_max; i++) for(j=0; j<matrix_max; j++) driver_matrix[i][j] = -99; // means that these entries of the driver matrix are not yet used
	driver_matrix[0][0] = 0; // driver matrix comprises one entry, which is zero (the distance from the first genotype to itself)

	printf("Initialised\n");
	printf("################################################################\n");
}

/////////////////////// simulation:

// run the simulation:
void run_sim(char *input_and_output_path, char *config_file_with_path)
{
	int i;
	int start_index, end_index;
	int trial_num;
	int num_demes, num_clones, num_cells, num_empty_demes;
	int num_matrix_cols, num_driver_matrix_cols; // number of columns in distance matrix
	int next_genotype_id, next_driver_genotype_id; // unique ID to assign to next genotype
	int num_empty_cols, num_empty_driver_cols; // number of columns in distance matrix that can be replaced (because the corresponding genotype has gone extinct)
	int num_extinct_genotypes, num_extinct_driver_genotypes; // number of genotypes that kept in the matrix and are extinct
	int init_driver_birth_mutations = 0, init_driver_mig_mutations = 0, init_passengers = 0; // initial number of drivers, migration 
	double gens_elapsed, gens_added; // counting time in cell generations
	double gens_timer_slow, gens_timer_grids, gens_timer_output; // timers for scheduling periodic tasks (e.g. writing to files)
	long idum; // seed of random number generator
	int chosen_clone, event_type;
	int num_samples_list[3] = {1, 2, 4};
	long t1 = (long)time(NULL); // for timing how long the program takes to run in seconds
	char *preamble_text = (char *) malloc(9999+1); // gnuplot code needed to set up each plot
	char *preamble_drivers_text = (char *) malloc(9999+1); // gnuplot code needed to set up each plot
	char *preamble_passengers_text = (char *) malloc(9999+1); // gnuplot code needed to set up each plot
	char *buffer_text_long = (char *) malloc(9999+1);
	char *buffer_text_short = (char *) malloc(20+1);
	int daughter_clone_nums[2];
	int parent_deme_num, chosen_deme, num_clones_in_deme;
	int max_layer_needed;
	int parent_geno_num, daughter_geno_num, parent_driver_geno_num, daughter_driver_geno_num;
	int new_birth_mutations[] = {0,0}, new_mig_mutations[] = {0,0}, new_passengers[] = {0,0};
	int new_mutations[2];
	int event_counter[7];
	int driver_counts[MAX_DRIVERS_TO_COUNT + 1];
	int cell_type, chosen;
	float temp_sum;
	int threshold_pop;
	int print_output;

	if(!filled_grid) threshold_pop = 5000; // output is written when population size first reaches this threshold
	else threshold_pop = 2 * max_pop; // unless the population size is not increasing

	// generate gnuplot code needed to set up plots:
	preamble_text = preamble(preamble_text, buffer_text_long);
	preamble_drivers_text = preamble_drivers(preamble_drivers_text, buffer_text_long);
	preamble_passengers_text = preamble(preamble_passengers_text, buffer_text_long);

	for(trial_num = 0; trial_num < MAX_TRIALS; trial_num++) {

		if(trial_num > 0) { // restart with a different random seed and empty output files
			seed += 1000;
			close_files();
		}

		idum = -(seed + 1);
		exit_code = 0;

		open_files(input_and_output_path);
		initiate_files(num_samples_list);

		/// initialise:
		initialise(&num_cells, &num_clones, &num_demes, &num_matrix_cols, &num_empty_cols, init_driver_birth_mutations, 
			init_driver_mig_mutations, init_passengers, &num_empty_driver_cols, &num_driver_matrix_cols, &next_driver_genotype_id, 
			&next_genotype_id, &idum, &num_extinct_genotypes, &num_empty_demes, &num_extinct_driver_genotypes);

		gens_timer_slow = 0;
		gens_timer_grids = 0;
		gens_timer_output = 0;
		gens_elapsed = 0;
		iterations = 0;
		for(i = 0; i < 7; i++) {
			event_counter[i] = 0;
			total_events[i] = 0;
		}

		// main loop:
		do{
			// optionally, passenger mutations stop occurring when population size reaches a threshold:
			if(num_cells >= passenger_pop_threshold) mu_passenger = 0;

			// periodically reset sums of rates, to get rid of accumulated rounding errors:
			if(!(iterations % 100000) && iterations > 0) {
				reset_deme_and_bintree_sums(num_demes, num_clones, 0);
				if(use_clone_bintrees) {
					check_clone_bintree_sums(FALSE);
					reset_clone_bintree_sums(num_demes, num_clones);
				}
			}

			// calculate general metrics, write to disc:
			if(gens_timer_output >= 1 || iterations == 0) {
				print_output = 0;
				main_calculations_and_output(&idum, num_demes, num_matrix_cols, num_driver_matrix_cols, gens_elapsed, driver_counts, num_cells, num_clones, 0, 
					num_samples_list, next_genotype_id, next_driver_genotype_id, event_counter, num_empty_cols, num_empty_driver_cols, num_extinct_genotypes, 
					num_extinct_driver_genotypes, num_empty_demes, t1, print_output, PHYLO_AND_POPS);
				gens_timer_slow += gens_timer_output;
				gens_timer_output = 0;
			}
			// calculate diversity metrics, write to disc and to screen:
			if(gens_timer_slow >= 10 || num_cells >= threshold_pop || iterations == 0) {
				if(gens_timer_slow >= 10) {
					print_output = 2;
					for(i=0; i<7; i++) total_events[i] += event_counter[i];
				}
				else print_output = 0;
				
				main_calculations_and_output(&idum, num_demes, num_matrix_cols, num_driver_matrix_cols, gens_elapsed, driver_counts, num_cells, num_clones, 0, 
					num_samples_list, next_genotype_id, next_driver_genotype_id, event_counter, num_empty_cols, num_empty_driver_cols, num_extinct_genotypes, 
					num_extinct_driver_genotypes, num_empty_demes, t1, print_output, DIVERSITIES);

				if(gens_timer_slow >= 10) {
					for(i=0; i<7; i++) event_counter[i] = 0; // reset event counters
					gens_timer_slow = 0;
				}

				if(num_cells >= threshold_pop) threshold_pop += 5000;
			}
			// periodically create images of the grid state:
			if(gens_timer_grids >= 10 || iterations == 0) {
				gens_timer_grids = 0;
				grids_output(preamble_text, preamble_drivers_text, preamble_passengers_text, gens_elapsed, num_clones, num_demes, input_and_output_path, buffer_text_short, buffer_text_long, FALSE);
			}

			// calculate sum of all rates (for calculating elapsed time):
			temp_sum = sum_of_all_rates(num_demes);

			// calculate elapsed time according to the Gillespie algorithm
			gens_added = expdev(&idum) / temp_sum; // number of generations elapsed in this loop
			
			gens_elapsed += gens_added;
			gens_timer_grids += gens_added; // timer for periodic events
			gens_timer_output += gens_added; // timer for periodic events

			// increase max_gens if needed (max_gens is needed for creating unique image file names):
			if(gens_elapsed > max_gens) max_gens *= 10;

			// choose bintree elements:
			max_layer_needed = max_layer_needed_array[num_demes];
			chosen = choose_bintree_doubles(bintree_deme_doubles, max_layer_needed, num_demes - 1, max_demes, &idum);

			// choose a deme, with weights determined by sums of rates:
			start_index = SET_SIZE * chosen;
			end_index = MIN(start_index + SET_SIZE - 1, num_demes - 1);
			chosen_deme = weighted_random_doubles(deme_doubles[SUM_RATES], &idum, start_index, end_index - start_index + 1);

			// choose whether event affects cancer cell or normal cell:
			buff_array[0] = (float)deme_doubles[SUM_BIRTH_RATES][chosen_deme] + deme_ints[POPULATION][chosen_deme] * deme_floats[DEATH_RATE][chosen_deme];
			if(migration_type == 1 || migration_type == 3) buff_array[0] += (float)deme_doubles[SUM_MIGRATION_RATES][chosen_deme];
			buff_array[1] = buff_array[0] + deme_ints[NORMAL_CELLS][chosen_deme] * (normal_birth_rate + deme_floats[DEATH_RATE][chosen_deme]);
			cell_type = weighted_random_known_sums_floats(buff_array, &idum, 2);
			
			// error check:
			check_rates_sum(buff_array[0], buff_array[1], chosen_deme, cell_type);
			if(exit_code != 0) break;

			// if the event concerns a cancer cell:
			if(cell_type == 0) {				
				
				// select clone and event_type:
				if(use_clone_bintrees) { // if using binary trees when selecting clone
					num_clones_in_deme = deme_ints[NUM_CLONES_IN_DEME][chosen_deme];
					max_layer_needed = max_layer_needed_array[num_clones_in_deme];
					// select birth or migation or death:
					chosen = choose_event_for_deme(chosen_deme, buff_array, &idum);

					if(chosen == DEATH_EVENT) { // if death:
						// select bintree element:
						chosen = choose_bintree_ints(bintree_clone_ints[chosen_deme], max_layer_needed, num_clones_in_deme - 1, max_clones_per_deme, &idum);
						// select clone:
						chosen_clone = choose_clone_with_bintree_by_population(chosen, chosen_deme, num_clones_in_deme, &idum);
						if(exit_code != 0) break;
						// set event_type:
						event_type = DEATH_EVENT;
					}
					else { // if other event:
						// select bintree element:
						chosen = choose_bintree_doubles(bintree_clone_doubles[chosen_deme], max_layer_needed, num_clones_in_deme - 1, max_clones_per_deme, &idum);
						// select clone:
						chosen_clone = choose_clone_with_bintree_by_rates(chosen, chosen_deme, num_clones_in_deme, &idum);
						if(exit_code != 0) break;
						// select event_type:
						event_type = choose_event_for_clone(FALSE, chosen_deme, chosen_clone, buff_array, &idum);
					}
				}
				else { // if not using binary trees when selecting clone
					// select a clone:
					chosen_clone = choose_clone_without_bintree(chosen_deme, buff_array, buff_array[0], &idum);
					if(exit_code != 0) break;
					// select an event type:
					event_type = choose_event_for_clone(TRUE, chosen_deme, chosen_clone, buff_array, &idum);
				}

				parent_geno_num = clone_ints[GENOTYPE][chosen_clone];
				parent_driver_geno_num = clone_ints[DRIVER_GENOTYPE][chosen_clone];
				parent_deme_num = clone_ints[DEME][chosen_clone];

				// cell division (and attempted migration or deme fission if appropriate migration_type):
				if(event_type == BIRTH_EVENT) {
					// cell division:
					cell_division(event_counter, &num_cells, parent_deme_num, new_passengers, new_mig_mutations, new_birth_mutations, new_mutations, &idum, chosen_clone, 
						parent_geno_num, daughter_clone_nums, &num_empty_cols, &num_matrix_cols, empty_cols, &num_clones, parent_driver_geno_num, &num_empty_driver_cols, 
						&num_driver_matrix_cols, empty_driver_cols, &next_driver_genotype_id, num_demes, &next_genotype_id, &num_extinct_genotypes, 
						&num_empty_demes, &num_extinct_driver_genotypes, gens_elapsed);
					// error check:
					check_geno_populations(chosen_clone, event_type);
					if(exit_code != 0) break;

					// attempted migration:
					if(migration_type == 0) {
						if(ran1(&idum) < 0.5) chosen_clone = daughter_clone_nums[0]; // migrating cell is chosen at random from the two daughters
						else chosen_clone = daughter_clone_nums[1];
						daughter_geno_num = clone_ints[GENOTYPE][chosen_clone];
						daughter_driver_geno_num = clone_ints[DRIVER_GENOTYPE][chosen_clone];

						// migration may be attempted, depdending on the genotype's migration rate:
						if(ran1(&idum) < genotype_floats[MIGRATION_RATE][daughter_geno_num]) cell_migration(event_counter, parent_deme_num, &idum, &num_demes, &num_clones, 
							chosen_clone, &num_cells, daughter_geno_num, daughter_driver_geno_num, &num_empty_demes, empty_cols, &num_empty_cols, &num_empty_driver_cols, empty_driver_cols, 
							&num_extinct_genotypes, &num_extinct_driver_genotypes);
						// error check:
						check_clone_populations(chosen_clone, event_type, parent_deme_num);
						if(exit_code != 0) break;
					}
					// attempted deme fission:
					else if(migration_type == 2 && deme_ints[POPULATION][parent_deme_num] + deme_ints[NORMAL_CELLS][parent_deme_num] >=K && 
						ran1(&idum) < deme_doubles[SUM_MIGRATION_RATES][parent_deme_num]) deme_fission(event_counter, parent_deme_num, &idum, &num_demes, &num_clones, 
						&num_cells, &num_empty_demes, &num_empty_cols, &num_empty_driver_cols, empty_cols, empty_driver_cols, &num_extinct_genotypes, &num_extinct_driver_genotypes, num_matrix_cols);
					// error check:
					check_geno_populations(chosen_clone, event_type);
				}
				
				// cell death:
				else if(event_type == DEATH_EVENT) cell_death(event_counter, &num_cells, parent_deme_num, parent_geno_num, empty_cols, &num_empty_cols, chosen_clone, &num_clones, 
					parent_driver_geno_num, &num_empty_driver_cols, empty_driver_cols, num_demes, &num_extinct_genotypes, &num_empty_demes, &num_extinct_driver_genotypes);
				
				// cell migration or deme fission:
				else if(deme_ints[POPULATION][parent_deme_num] > 1) {
					if(migration_type == 1) cell_migration(event_counter, parent_deme_num, &idum, &num_demes, &num_clones, chosen_clone, &num_cells, parent_geno_num, parent_driver_geno_num, &num_empty_demes, 
						empty_cols, &num_empty_cols, &num_empty_driver_cols, empty_driver_cols, &num_extinct_genotypes, &num_extinct_driver_genotypes);
					else if(deme_ints[POPULATION][parent_deme_num] + deme_ints[NORMAL_CELLS][parent_deme_num] >=K) deme_fission(event_counter, parent_deme_num, &idum, &num_demes, &num_clones, &num_cells, &num_empty_demes, 
						&num_empty_cols, &num_empty_driver_cols, empty_cols, empty_driver_cols, &num_extinct_genotypes, &num_extinct_driver_genotypes, num_matrix_cols);
				}

				// error check:
				check_geno_populations(chosen_clone, event_type);
				if(exit_code != 0) break; // program has achieved a halting condition
			}

			else { // if the event concerns a normal cell:
				// assign weights for choosing which event to perform:
				buff_array[0] = normal_birth_rate;
				buff_array[1] = buff_array[0] + deme_floats[DEATH_RATE][chosen_deme];

				// determine event type (birth or death):
				event_type = weighted_random_known_sums_floats(buff_array, &idum, 2);

				// cell division:
				if(event_type == BIRTH_EVENT) add_or_remove_normal_cell(1, chosen_deme, num_cells, event_counter, num_demes);
				// cell death:
				else add_or_remove_normal_cell(-1, chosen_deme, num_cells, event_counter, num_demes);
				// error check:
				check_normal_pops(chosen_deme);
				if(exit_code != 0) break; // program has achieved a halting condition
			}

			// error check:
			check_genotype_counts(num_matrix_cols, num_empty_cols, num_extinct_genotypes, num_driver_matrix_cols, num_empty_driver_cols, num_extinct_driver_genotypes, cell_type, event_type, chosen_deme);
			if(exit_code != 0) break; // program has achieved a halting condition

			iterations++;
		
		}while(num_cells < max_pop && num_cells > 0 && (long)time(NULL)-t1 < max_time && gens_elapsed < max_generations);

		end_of_loop_output(num_cells, gens_elapsed, t1);

		if(num_cells >= MIN(20, max_pop) || filled_grid) break;
	}

	if(num_cells > 0) {
		main_calculations_and_output(&idum, num_demes, num_matrix_cols, num_driver_matrix_cols, gens_elapsed, driver_counts, num_cells, num_clones, 1, 
			num_samples_list, next_genotype_id, next_driver_genotype_id, event_counter, num_empty_cols, num_empty_driver_cols, num_extinct_genotypes, 
			num_extinct_driver_genotypes, num_empty_demes, t1, 1, PHYLO_AND_POPS + DIVERSITIES + GENOPROPS);
		grids_output(preamble_text, preamble_drivers_text, preamble_passengers_text, gens_elapsed, num_clones, num_demes, input_and_output_path, buffer_text_short, buffer_text_long, TRUE);
	}

	final_output(trial_num, num_cells, t1);
	close_files();

	free(preamble_text);
	free(preamble_drivers_text);
	free(preamble_passengers_text);
	free(buffer_text_long);
	free(buffer_text_short);
}

/////////////////////// choose update type:

// randomly choose an element of a bintree:
int choose_bintree_doubles(double *bintree, int max_layer_needed, int max_array_index, int array_size, long *idum)
{
	int layer;
	int chosen = 0;
	int start_index, end_index;
	int min_bintree_index_this_layer, min_bintree_index_previous_layer = 0, final_bintree_index;
	
	for(layer = max_layer_needed; layer >= 0; layer--) { // loop over all bintree layers, starting at the top
		if(array_size == max_demes) min_bintree_index_this_layer = min_deme_bintree_index_by_layer[layer];
		else if(array_size == max_clones_per_deme) min_bintree_index_this_layer = min_clone_bintree_index_by_layer[layer];
		else min_bintree_index_this_layer = get_bintree_index(layer, 0, array_size);
		
		start_index = min_bintree_index_this_layer + SET_SIZE * (chosen - min_bintree_index_previous_layer);

		if(array_size == max_demes) {
			final_bintree_index = get_deme_bintree_index(layer, max_array_index);
			end_index = MIN(start_index + SET_SIZE - 1, final_bintree_index);
		}
		else if(array_size == max_clones_per_deme) {
			final_bintree_index = get_clone_bintree_index(layer, max_array_index);
			end_index = MIN(start_index + SET_SIZE - 1, final_bintree_index);
		}
		else {
			final_bintree_index = get_bintree_index(layer, max_array_index, array_size);
			end_index = MIN(start_index + SET_SIZE - 1, final_bintree_index);
		}
		// choose a bintree element from the current layer, with weights determined by sums of rates:
		chosen = weighted_random_doubles(bintree, idum, start_index, end_index - start_index + 1);
		min_bintree_index_previous_layer = min_bintree_index_this_layer;
	}

	return chosen;
}

// randomly choose an element of a bintree:
int choose_bintree_ints(int *bintree, int max_layer_needed, int max_array_index, int array_size, long *idum)
{
	int layer;
	int chosen = 0;
	int start_index, end_index;
	int min_bintree_index_this_layer, min_bintree_index_previous_layer = 0, final_bintree_index;
	
	for(layer = max_layer_needed; layer >= 0; layer--) { // loop over all bintree layers, starting at the top
		if(array_size == max_demes) min_bintree_index_this_layer = min_deme_bintree_index_by_layer[layer];
		else if(array_size == max_clones_per_deme) min_bintree_index_this_layer = min_clone_bintree_index_by_layer[layer];
		else min_bintree_index_this_layer = get_bintree_index(layer, 0, array_size);
		
		start_index = min_bintree_index_this_layer + SET_SIZE * (chosen - min_bintree_index_previous_layer);
		
		if(array_size == max_demes) {
			final_bintree_index = get_deme_bintree_index(layer, max_array_index);
			end_index = MIN(start_index + SET_SIZE - 1, final_bintree_index);
		}
		else if(array_size == max_clones_per_deme) {
			final_bintree_index = get_clone_bintree_index(layer, max_array_index);
			end_index = MIN(start_index + SET_SIZE - 1, final_bintree_index);
		}
		else {
			final_bintree_index = get_bintree_index(layer, max_array_index, array_size);
			end_index = MIN(start_index + SET_SIZE - 1, final_bintree_index);
		}
		// choose a bintree element from the current layer, with weights determined by sums of rates:
		chosen = weighted_random_ints(bintree, idum, start_index, end_index - start_index + 1);
		min_bintree_index_previous_layer = min_bintree_index_this_layer;
	}

	return chosen;
}

// randomly choose a clone within a deme:
int choose_clone_without_bintree(int chosen_deme, float *buff_array, float deme_sum_cancer, long *idum)
{
	int i, clone_index, genotype_index, pop, chosen;

	for(i = 0; i < deme_ints[NUM_CLONES_IN_DEME][chosen_deme]; i++) { // loop over all clones in the chosen deme
		clone_index = clones_list_in_deme[chosen_deme][i];
		genotype_index = clone_ints[GENOTYPE][clone_index];
		pop = clone_ints[POPULATION][clone_index];
		// buff_array[i] contains cumulative sum of all cancer cell rates in clones 0..i
		buff_array[i] = pop * (genotype_floats[BIRTH_RATE][genotype_index] + deme_floats[DEATH_RATE][chosen_deme]);
		if(migration_type == 1 || migration_type == 3) buff_array[i] += pop * genotype_floats[MIGRATION_RATE][genotype_index];
		if(i > 0) buff_array[i] += buff_array[i - 1];
	}

	// error check:
	check_deme_sums(chosen_deme, deme_sum_cancer);

	chosen = weighted_random_known_sums_floats(buff_array, idum, deme_ints[NUM_CLONES_IN_DEME][chosen_deme]);
	return clones_list_in_deme[chosen_deme][chosen];
}

int choose_clone_with_bintree_by_population(int bintree_index, int chosen_deme, int num_clones_in_deme, long *idum)
{
	int chosen_clone;
	int start_index = SET_SIZE * bintree_index;
	int end_index = MIN(start_index + SET_SIZE - 1, num_clones_in_deme - 1);
	
	chosen_clone = weighted_random_ints(clones_pops_in_deme[chosen_deme], idum, start_index, end_index - start_index + 1);

	// error check:
	check_chosen_clone(chosen_clone, num_clones_in_deme);

	return clones_list_in_deme[chosen_deme][chosen_clone];
}

int choose_clone_with_bintree_by_rates(int bintree_index, int chosen_deme, int num_clones_in_deme, long *idum)
{
	int chosen_clone;
	int start_index = SET_SIZE * bintree_index;
	int end_index = MIN(start_index + SET_SIZE - 1, num_clones_in_deme - 1);
	
	chosen_clone = weighted_random_doubles(clones_rates_in_deme[chosen_deme], idum, start_index, end_index - start_index + 1);

	// error check:
	check_chosen_clone(chosen_clone, num_clones_in_deme);

	return clones_list_in_deme[chosen_deme][chosen_clone];
}

// randomly choose an event type:
// if include_death = 0 then function returns 0 for birth, or 2 for migration;
// if include_death = 1 then function returns 0 for birth, 1 for death, or 2 for migration;
int choose_event_for_clone(bool include_death, int chosen_deme, int chosen_clone, float *buff_array, long *idum)
{
	int buff_array_length = 0;
	int genotype_index = clone_ints[GENOTYPE][chosen_clone];
	int res;

	// assign weights for choosing which event to perform (birth, death or migration):
	buff_array[0] = genotype_floats[BIRTH_RATE][genotype_index];
	buff_array_length++;

	if(include_death) {
		buff_array[buff_array_length] = buff_array[buff_array_length - 1] + deme_floats[DEATH_RATE][chosen_deme];
		buff_array_length++;
	}

	if(migration_type == 1 || migration_type == 3) {
		buff_array[buff_array_length] = buff_array[buff_array_length - 1] + genotype_floats[MIGRATION_RATE][clone_ints[GENOTYPE][chosen_clone]];
		buff_array_length++;
	}

	// determine event type:
	res = weighted_random_known_sums_floats(buff_array, idum, buff_array_length);

	if(include_death) return res;
	else return res * 2; // so that migration corresponds to 2
}

// randomly choose an event type:
// returns 0 for birth or migration, 1 for death;
int choose_event_for_deme(int chosen_deme, float *buff_array, long *idum)
{
	buff_array[0] = deme_doubles[SUM_BIRTH_RATES][chosen_deme];
	if(migration_type == 1 || migration_type == 3) buff_array[0] += deme_doubles[SUM_MIGRATION_RATES][chosen_deme];
	buff_array[1] = buff_array[0] + deme_floats[DEATH_RATE][chosen_deme] * deme_ints[POPULATION][chosen_deme];

	// determine event type:
	return weighted_random_known_sums_floats(buff_array, idum, 2);
}

/////////////////////// cell events:

// division of a cancer cell:
void cell_division(int *event_counter, int *num_cells, int parent_deme_num, int *new_passengers, int *new_mig_mutations, int *new_birth_mutations, 
	int *new_mutations, long *idum, int parent_clone, int parent_geno_num, int *daughter_clone_nums, int *num_empty_cols, int *num_matrix_cols, int *empty_cols, 
	int *num_clones, int parent_driver_geno_num, int *num_empty_driver_cols, int *num_driver_matrix_cols, int *empty_driver_cols, int *next_driver_genotype_id, 
	int num_demes, int *next_genotype_id, int *num_extinct_genotypes, int *num_empty_demes, int *num_extinct_driver_genotypes, float gens_elapsed)
{
	int index, start_index, end_index;
	int daughter_geno_num, daughter_driver_geno_num;
	float new_birth_rate, new_migration_rate;
	int daughter_driver_id;

	// record the daughter's clone number (for use in subsequent functions);
	// unless changed by mutations, both daughter cells will belong to the parent's clone:
	daughter_clone_nums[0] = parent_clone;
	daughter_clone_nums[1] = parent_clone;

	// error check:
	check_geno_populations(parent_clone, 0);
	if(exit_code != 0) return;

	event_counter[BIRTH_EVENT]++;

	// increment total and deme populations:
	*num_cells = *num_cells + 1;
	increment_or_decrement_deme(1, parent_deme_num, *num_cells, num_demes, num_empty_demes);

	// choose numbers of mutations of each type:
	choose_number_mutations(new_passengers, new_mig_mutations, new_birth_mutations, new_mutations, idum);

	// if no new mutations then simply increment parent clone, genotype, and driver genotype populations:
	if(!(new_mutations[0]) && !(new_mutations[1])) {
		increment_or_decrement_clone(parent_clone, parent_deme_num, num_clones, +1, num_demes);
		increment_or_decrement_genotype(genotype_ints, parent_geno_num, empty_cols, num_empty_cols, +1, num_extinct_genotypes);
		increment_or_decrement_genotype(driver_genotype_ints, parent_driver_geno_num, empty_driver_cols, num_empty_driver_cols, +1, num_extinct_driver_genotypes);
		daughter_geno_num = parent_geno_num;
		daughter_driver_geno_num = parent_driver_geno_num;
	}

	else { // it at least one daughter has at least one mutation:
		
		// if only the first daughter has one or more new mutations then index = 0;
		// if only the second daughter has one or more new mutations then index = 1;
		// if both daughters have one or more new mutations then index = 0 and 1;
		start_index = 1;
		if(new_mutations[0]) start_index = 0;
		end_index = 0;
		if(new_mutations[1]) end_index = 1;

		for(index = start_index; index <= end_index; index++) { // loop over the mutated daughter cell(s)
			event_counter[MUTATION_EVENT]++;

			new_birth_rate = set_birth_rate(new_birth_mutations[index], new_passengers[index], genotype_floats[BIRTH_RATE][parent_geno_num], idum);
			new_migration_rate = set_migration_rate(new_mig_mutations[index], genotype_floats[MIGRATION_RATE][parent_geno_num], idum);
			
			// if the daughter has at least one driver mutation then create a new driver genotype and add a column to the driver distance matrix:
			if(new_birth_mutations[index] || new_mig_mutations[index]) {
				daughter_driver_geno_num = select_genotype_index(num_empty_driver_cols, num_driver_matrix_cols, empty_driver_cols);
				daughter_driver_id = *next_driver_genotype_id;
				create_genotype(driver_genotype_ints, driver_genotype_floats, num_driver_matrix_cols, daughter_driver_geno_num, parent_driver_geno_num, 
					next_driver_genotype_id, daughter_driver_id, new_birth_rate, new_migration_rate, new_passengers[index], new_birth_mutations[index], new_mig_mutations[index], gens_elapsed);
				if(matrix_max > 0) create_column(driver_matrix, *num_driver_matrix_cols, parent_driver_geno_num, daughter_driver_geno_num, new_birth_mutations[index] + new_mig_mutations[index]);
			}
			else {
				daughter_driver_geno_num = parent_driver_geno_num; // if the daughter has no driver mutations then it belongs to its parent's driver genotype
				daughter_driver_id = driver_genotype_ints[DRIVER_IDENTITY][parent_driver_geno_num];
			}
			
			// determine daughter genotype index (either at end of the array, or replacing an extinct entry):
			daughter_geno_num = select_genotype_index(num_empty_cols, num_matrix_cols, empty_cols);
			// add new column to genotype distance matrix (if required):
			if(record_matrix) create_column(matrix, *num_matrix_cols, parent_geno_num, daughter_geno_num, new_mutations[index]);
			// create a new genotype containing only the daughter cell:
			create_genotype(genotype_ints, genotype_floats, num_matrix_cols, daughter_geno_num, parent_geno_num, next_genotype_id, daughter_driver_id, new_birth_rate, new_migration_rate, 
				new_passengers[index], new_birth_mutations[index], new_mig_mutations[index], gens_elapsed);

			// create a clone containing only the daughter cell:
			create_clone(daughter_geno_num, num_clones, parent_deme_num, daughter_driver_geno_num, num_demes, 1);
			if(exit_code != 0) return;
			// record the daughter's new clone number (for use in subsequent functions):
			daughter_clone_nums[index] = *num_clones - 1;
		}

		// if both daughters have only passenger mutations then daughter is added to parent's driver genotype:
		if(!(new_birth_mutations[0] + new_mig_mutations[0]) && !(new_birth_mutations[1] + new_mig_mutations[1])) driver_genotype_ints[POPULATION][parent_driver_geno_num]++;
		// if both daughters have driver mutations then decrement parent's driver genotype:
		else if(new_birth_mutations[0] + new_mig_mutations[0] > 0 && new_birth_mutations[1] + new_mig_mutations[1] > 0) increment_or_decrement_genotype(driver_genotype_ints, 
			parent_driver_geno_num, empty_driver_cols, num_empty_driver_cols, -1, num_extinct_driver_genotypes);
			
		// if both daughters have mutations then remove parent from genotype and clone:
		if(new_mutations[0] && new_mutations[1]) {
			increment_or_decrement_genotype(genotype_ints, parent_geno_num, empty_cols, num_empty_cols, -1, num_extinct_genotypes);
			if(clone_ints[POPULATION][parent_clone] == 1) daughter_clone_nums[1] = parent_clone; // parent clone is about to go extinct and its clone number will be reassigned to its daughter
			increment_or_decrement_clone(parent_clone, parent_deme_num, num_clones, -1, num_demes);
		}
	}
	// error check:
	check_geno_populations(parent_clone, 0);
	check_matrix_cols1(*num_matrix_cols, *num_driver_matrix_cols, *num_cells);
	if(exit_code != 0) return;
}

// death of a cancer cell:
void cell_death(int *event_counter, int *num_cells, int parent_deme_num, int parent_geno_num, int *empty_cols, int *num_empty_cols, int chosen_clone, int *num_clones, int parent_driver_geno_num, 
	int *num_empty_driver_cols, int *empty_driver_cols, int num_demes, int *num_extinct_genotypes, int *num_empty_demes, int *num_extinct_driver_genotypes)
{
	event_counter[DEATH_EVENT]++;

	// remove from cell count, deme, and genotype:
	*num_cells = *num_cells - 1;
	increment_or_decrement_deme(-1, parent_deme_num, *num_cells, num_demes, num_empty_demes);
	increment_or_decrement_genotype(genotype_ints, parent_geno_num, empty_cols, num_empty_cols, -1, num_extinct_genotypes);
	increment_or_decrement_genotype(driver_genotype_ints, parent_driver_geno_num, empty_driver_cols, num_empty_driver_cols, -1, num_extinct_driver_genotypes);

	// remove from old clone, and update deme and bintree sums of rates:
	increment_or_decrement_clone(chosen_clone, parent_deme_num, num_clones, -1, num_demes);
}

// migration of a cancer cell:
void cell_migration(int *event_counter, int origin_deme_num, long *idum, int *num_demes, int *num_clones, int clone_num, int *num_cells,
	int geno_num, int driver_geno_num, int *num_empty_demes, int *empty_cols, int *num_empty_cols, int *num_empty_driver_cols, int *empty_driver_cols, 
	int *num_extinct_genotypes, int *num_extinct_driver_genotypes)
{
	int new_deme_num;
	int existing_clone_index;
	int new_x, new_y;

	// choose the cooridinates of the destination deme:
	choose_grid_square(origin_deme_num, idum, &new_x, &new_y);

	// check if the new cell has gone beyond the edge of the grid:
	if(new_x < 0 || new_x >= dim_grid || new_y < 0 || new_y >= dim_grid) {
		if(!well_mixed && !filled_grid) {
			printf("Outgrown the grid!\n");
			fprintf(error_log, "Outgrown the grid!\n");
			exit_code = 1; // program will terminate
			return;
		}
		else if(filled_grid) {
			cell_death(event_counter, num_cells, origin_deme_num, geno_num, empty_cols, num_empty_cols, clone_num, num_clones, driver_geno_num, 
				num_empty_driver_cols, empty_driver_cols, *num_demes, num_extinct_genotypes, num_empty_demes, num_extinct_driver_genotypes);
			return;
		}
	}

	// if migrating cell stays within same deme then the migration attempts fails:
	if(new_x == deme_ints[XCOORD][origin_deme_num] && new_y == deme_ints[YCOORD][origin_deme_num]) return;
	
	new_deme_num = grid[new_x][new_y]; // index of the destination deme

	// if migration should occur only at the edge, and the destination deme is already too full then the migration attempts fails:
	if(new_deme_num != EMPTY) if(migration_edge_only && ran1(idum) < deme_floats[MIGRATION_MODIFIER][new_deme_num]) return;

	// if new deme has not yet been occupied then initiate the deme:
	if(new_deme_num == EMPTY) {
		new_deme_num = *num_demes;
		create_deme(new_x, new_y, num_demes, *num_cells);
	}
	// if the destination deme has been occupied but currently contains no cancer cells then decrement the count of empty demes:
	else if(deme_ints[POPULATION][new_deme_num] == 0) *num_empty_demes = *num_empty_demes - 1;

	event_counter[MIGRATION_EVENT]++;

	// add to new deme:
	increment_or_decrement_deme(1, new_deme_num, *num_cells, *num_demes, num_empty_demes);
	
	// determine whether the destination deme already contains a clone that matches the migrating cell's genotype:
	existing_clone_index = select_clone_in_deme(new_deme_num, geno_num);
	// if no matching clone exists then create a new clone:
	if(existing_clone_index == -1) create_clone(geno_num, num_clones, new_deme_num, driver_geno_num, *num_demes, 1);
	// if a matching clone already exists then add the migrating cell to that clone, and update deme and bintree sums of rates:
	else increment_or_decrement_clone(existing_clone_index, new_deme_num, num_clones, +1, *num_demes);
	if(exit_code != 0) return;

	// remove from old deme and clone, and update deme and bintree sums of rates:
	increment_or_decrement_deme(-1, origin_deme_num, *num_cells, *num_demes, num_empty_demes);
	increment_or_decrement_clone(clone_num, origin_deme_num, num_clones, -1, *num_demes);

	// error check:
	check_clones_in_deme1(origin_deme_num, clone_num);
	if(exit_code != 0) return;
}

// fission of a deme:
void deme_fission(int *event_counter, int origin_deme_num, long *idum, int *num_demes, int *num_clones, int *num_cells, int *num_empty_demes, int *num_empty_cols, int *num_empty_driver_cols, 
	int *empty_cols, int *empty_driver_cols, int *num_extinct_genotypes, int *num_extinct_driver_genotypes, int num_matrix_cols)
{
	int new_deme_index = *num_demes;
	int old_x = deme_ints[XCOORD][origin_deme_num], old_y = deme_ints[YCOORD][origin_deme_num];
	int x_to_fill, y_to_fill;
	int dividing_beyond_the_edge = 0;

	// find the coordinates of the site to be filled:
	get_deme_coordinates(&x_to_fill, &y_to_fill, old_x, old_y, idum);

	// when budging is allowed, exit if the site to be filled is already occupied, 
	// unless it is allowed to budge demes beyond the edge of the grid (i.e. when the grid is always entirely filled):
	if(!migration_edge_only) {
		if(grid[x_to_fill][y_to_fill] != EMPTY) {
			if(!filled_grid) {
				printf("Outgrown the grid! num_cells %d\n", *num_cells);
				fprintf(error_log, "Outgrown the grid! num_cells %d\n", *num_cells);
				exit_code = 1; // program will terminate
				return;
			}
			// test whether deme is on the edge of the grid and its daughter deme will be placed beyond the edge:
			else if(old_x == x_to_fill && old_y == y_to_fill) dividing_beyond_the_edge = 1;
			// if a deme is pushed off the edge of the grid (not dividing beyond the edge) then remove it and all of its contents:
			else {
				if(origin_deme_num == *num_demes - 1) origin_deme_num = grid[x_to_fill][y_to_fill]; // anticipating change in deme indexing, as a result of remove_deme
				remove_deme(grid[x_to_fill][y_to_fill], num_cells, num_clones, num_demes, num_empty_demes, empty_cols, num_empty_cols, num_empty_driver_cols, empty_driver_cols, 
					num_extinct_genotypes, num_extinct_driver_genotypes);
			}
		}
	}
	// when budging is not allowed, exit if a deme wants to divide beyond the edge of the grid, 
	// unless it is allowed to do so (i.e. when the grid is always entirely filled):
	else {
		if(x_to_fill < 0 || x_to_fill >= dim_grid || y_to_fill < 0 || y_to_fill >= dim_grid) {
			if(!filled_grid) {
				printf("Outgrown the grid! num_cells %d\n", *num_cells);
				fprintf(error_log, "Outgrown the grid! num_cells %d\n", *num_cells);
				exit_code = 1; // program will terminate
				return;
			}
			else dividing_beyond_the_edge = 1;
		}
		// if budging is not allowed and the neighbour deme is already occupied then fission attempt fails:
		else if(grid[x_to_fill][y_to_fill] != EMPTY) return;
	}

	// budge demes, starting at the far end of the line:
	if(!migration_edge_only && !dividing_beyond_the_edge) budge_demes(old_x, old_y, &x_to_fill, &y_to_fill);

	// create new deme in empty neighbouring grid square, unless deme is dividing beyond the edge:
	if(!dividing_beyond_the_edge) create_deme(x_to_fill, y_to_fill, num_demes, *num_cells);

	// move cells from the dividing deme into the newly created deme (or remove cells if deme is dividing beyond the edge):
	move_cells(idum, origin_deme_num, dividing_beyond_the_edge, new_deme_index, num_cells, num_demes, num_empty_demes, num_clones, 
		num_empty_cols, num_empty_driver_cols, num_extinct_genotypes, num_extinct_driver_genotypes, event_counter);
	
	event_counter[FISSION_EVENT]++;
}

///// genotype or driver genotype:

// choose numbers of mutations of each type (passenger, birth rate, migration rate):
void choose_number_mutations(int *new_passengers, int *new_mig_mutations, int *new_birth_mutations, int *new_mutations, long *idum)
{
	int i;

	for(i=0; i<2; i++) { // loop over both daughter cells
		new_birth_mutations[i] = poisson(idum, mu_driver_birth); // driver mutations
		new_passengers[i] = poisson(idum, mu_passenger); // passenger mutations
		new_mig_mutations[i] = poisson(idum, mu_driver_migration); // migration rate mutations
		new_mutations[i] = new_birth_mutations[i] + new_passengers[i] + new_mig_mutations[i];
	}
}

// determine new genotype or driver genotype index (either at end of the array, or replacing an extinct entry):
int select_genotype_index(int *num_empty_cols, int *num_matrix_cols, int *empty_cols)
{
	int daughter_geno_num;

	if(*num_empty_cols) { // new column replaces empty column
		daughter_geno_num = empty_cols[*num_empty_cols - 1];
		*num_empty_cols = *num_empty_cols - 1;
	}
	else { // new column is added to end of matrix
		daughter_geno_num = *num_matrix_cols;
		*num_matrix_cols = *num_matrix_cols + 1;
	}

	return daughter_geno_num;
}

// create a new genotype or driver genotype:
void create_genotype(int **geno_or_driver_ints, float **geno_or_driver_floats, int *num_matrix_cols, int daughter_geno_num, int parent_geno_num, int *next_genotype_id, int daughter_driver_id, 
	float new_birth_rate, float new_migration_rate, int new_passengers, int new_birth_mutations, int new_mig_mutations, float gens_elapsed)
{
	// create new genotype:
	geno_or_driver_ints[POPULATION][daughter_geno_num] = 1; // initialise population

	// record new mutations:
	geno_or_driver_ints[NUM_PASSENGER_MUTATIONS][daughter_geno_num] = geno_or_driver_ints[NUM_PASSENGER_MUTATIONS][parent_geno_num] + new_passengers;
	geno_or_driver_ints[NUM_DRIVER_MUTATIONS][daughter_geno_num] = geno_or_driver_ints[NUM_DRIVER_MUTATIONS][parent_geno_num] + new_birth_mutations;
	geno_or_driver_ints[NUM_MIGRATION_MUTATIONS][daughter_geno_num] = geno_or_driver_ints[NUM_MIGRATION_MUTATIONS][parent_geno_num] + new_mig_mutations;
	
	// set birth and migration rates:
	geno_or_driver_floats[BIRTH_RATE][daughter_geno_num] = new_birth_rate;
	geno_or_driver_floats[MIGRATION_RATE][daughter_geno_num] = new_migration_rate;

	// set unique ID:
	geno_or_driver_ints[IDENTITY][daughter_geno_num] = *next_genotype_id;
	*next_genotype_id = *next_genotype_id + 1;

	// set driver ID:
	geno_or_driver_ints[DRIVER_IDENTITY][daughter_geno_num] = daughter_driver_id;

	// record parent's ID:
	geno_or_driver_ints[PARENT][daughter_geno_num] = geno_or_driver_ints[IDENTITY][parent_geno_num];

	// parent genotype record will never be overwritten but daughter genotype record may be overwritten:
	geno_or_driver_ints[IMMORTAL][parent_geno_num] = 1; // don't overwrite this genotype
	geno_or_driver_ints[IMMORTAL][daughter_geno_num] = 0;

	geno_or_driver_floats[ORIGIN_TIME][daughter_geno_num] = gens_elapsed;
}

// increment or decrement a genotype or driver genotype:
void increment_or_decrement_genotype(int **geno_or_driver_ints, int parent_geno_num, int *empty_cols, int *num_empty_cols, int change, int *num_extinct_genotypes)
{
	geno_or_driver_ints[POPULATION][parent_geno_num] += change;
	
	// if the genotype is becoming extinct and can be overwritten:
	if(geno_or_driver_ints[POPULATION][parent_geno_num] == 0 && geno_or_driver_ints[IMMORTAL][parent_geno_num] == 0) {
		empty_cols[*num_empty_cols] = parent_geno_num; // record that this genotype can be overwritten
		*num_empty_cols = *num_empty_cols + 1; // add to the count of extinct genotypes that can be overwritten
	}

	// if the genotype is becoming extinct but should not be overwritten then add to the count of extinct genotypes that should not be overwritten:
	if(geno_or_driver_ints[POPULATION][parent_geno_num] == 0 && geno_or_driver_ints[IMMORTAL][parent_geno_num] == 1) *num_extinct_genotypes += 1;

	// if genotype population has reached a reasonably large size then it will never be overwritten:
	if(geno_or_driver_ints[POPULATION][parent_geno_num] >= 10) geno_or_driver_ints[IMMORTAL][parent_geno_num] = 1;
}

// create new column of distance matrix or driver distance matrix:
void create_column(int **either_matrix, int num_matrix_cols, int parent_geno_num, int daughter_geno_num, int num_mutations)
{ // this function updates matrix[][] or driver_matrix[][] only
	int i;
	int buff[num_matrix_cols + 1];

	// record values of new column, based on parent column:
	for(i=0; i<= parent_geno_num; i++) buff[i] = either_matrix[i][parent_geno_num] + num_mutations;
	for(i=parent_geno_num+1; i< num_matrix_cols; i++) buff[i] = either_matrix[parent_geno_num][i] + num_mutations;

	// either add a new column or replace a column corresponding to an extinct genotype:
	for(i=0; i<= daughter_geno_num; i++) either_matrix[i][daughter_geno_num] = buff[i];
	for(i=daughter_geno_num+1; i< num_matrix_cols; i++) either_matrix[daughter_geno_num][i] = buff[i];
	if(parent_geno_num < daughter_geno_num) either_matrix[parent_geno_num][daughter_geno_num] = num_mutations;
	else either_matrix[daughter_geno_num][parent_geno_num] = num_mutations;
	either_matrix[daughter_geno_num][daughter_geno_num] = 0;
}

///// clone:

// determine whether a deme already contains a clone with a particular genotype:
int select_clone_in_deme(int deme_index, int genotype_index)
{
	int i, clone_index;

	// check all clones to see if the "new" clone created by migration already exists, and if it does exist then return its index:
	for(i = 0; i < deme_ints[NUM_CLONES_IN_DEME][deme_index]; i++) {
		clone_index = clones_list_in_deme[deme_index][i];
		if(clone_ints[GENOTYPE][clone_index] == genotype_index) return clone_index;
	}
	
	return -1; // there is no matching clone
}

// create a clone, and update deme and bintree sums of rates:
void create_clone(int daughter_geno_num, int *num_clones, int daughter_deme_num, int daughter_driver_geno_num, int num_demes, int initial_clone_pop)
{
	int daughter_clone_num = *num_clones;
	int max_layer_needed;
	int num_clones_in_deme = deme_ints[NUM_CLONES_IN_DEME][daughter_deme_num];

	// create new clone:
	clone_ints[POPULATION][daughter_clone_num] = initial_clone_pop; // initialise population
	clone_ints[DEME][daughter_clone_num] = daughter_deme_num; // point to deme index
	clone_ints[GENOTYPE][daughter_clone_num] = daughter_geno_num; // point to genotype index
	clone_ints[DRIVER_GENOTYPE][daughter_clone_num] = daughter_driver_geno_num; // point to driver genotype index

	update_deme_bintree_phenotype_rates(daughter_clone_num, daughter_deme_num, initial_clone_pop, num_demes);

	// add to list of clones in deme:
	clone_ints[INDEX_IN_DEME][daughter_clone_num] = num_clones_in_deme;
	clones_list_in_deme[daughter_deme_num][num_clones_in_deme] = daughter_clone_num;
	if(use_clone_bintrees) set_clone_in_deme(daughter_deme_num, daughter_geno_num, num_clones_in_deme, initial_clone_pop);

	deme_ints[NUM_CLONES_IN_DEME][daughter_deme_num]++;

	// error check:
	check_clones_in_deme2(num_clones_in_deme);
	if(exit_code != 0) return;
	
	*num_clones = *num_clones + 1; // increment count of clones

	num_clones_in_deme++;

	if(use_clone_bintrees) {
		max_layer_needed = max_layer_needed_array[num_clones_in_deme];

		extend_bintree_int(bintree_clone_ints[daughter_deme_num], num_clones_in_deme - 1, max_layer_needed, max_clones_per_deme);
		update_bintree_layers_int(clone_ints[POPULATION][daughter_clone_num], bintree_clone_ints[daughter_deme_num], num_clones_in_deme - 1, max_layer_needed, max_clones_per_deme);

		extend_bintree(bintree_clone_doubles[daughter_deme_num], num_clones_in_deme - 1, max_layer_needed, max_clones_per_deme);
		update_bintree_layers(clones_rates_in_deme[daughter_deme_num][num_clones_in_deme - 1], bintree_clone_doubles[daughter_deme_num], num_clones_in_deme - 1, max_layer_needed, max_clones_per_deme);
	}
}

// increment or decrement a clone, and update deme and bintree sums of rates for changes in birth and migration rates:
void increment_or_decrement_clone(int chosen_clone, int deme_index, int *num_clones, int change, int num_demes)
{
	int num_clones_in_deme = deme_ints[NUM_CLONES_IN_DEME][deme_index];
	int max_layer_needed = max_layer_needed_array[num_clones_in_deme];
	int index_in_deme = clone_ints[INDEX_IN_DEME][chosen_clone];
	double rates_change;
	int old_pop = clone_ints[POPULATION][chosen_clone];

	// decrement clone population:
	clone_ints[POPULATION][chosen_clone] += change;

	update_deme_bintree_phenotype_rates(chosen_clone, deme_index, change, num_demes);

	if(use_clone_bintrees) {
		rates_change = change * clones_rates_in_deme[deme_index][index_in_deme] / old_pop;

		clones_pops_in_deme[deme_index][index_in_deme] += change;
		clones_rates_in_deme[deme_index][index_in_deme] += rates_change;

		update_bintree_layers_int(change, bintree_clone_ints[deme_index], index_in_deme, max_layer_needed, max_clones_per_deme);
		
		update_bintree_layers(rates_change, bintree_clone_doubles[deme_index], index_in_deme, max_layer_needed, max_clones_per_deme);
	}

	// if clone is becoming extinct:
	if(clone_ints[POPULATION][chosen_clone] == 0) remove_clone(chosen_clone, deme_index, num_clones, num_clones_in_deme);
}

// remove a clone from the list of clones, and update the lists of clones in demes:
void remove_clone(int chosen_clone, int deme_index, int *num_clones, int num_clones_in_deme)
{
	int i;
	int index_in_deme = clone_ints[INDEX_IN_DEME][chosen_clone];
	int max_layer_needed = max_layer_needed_array[num_clones_in_deme];

	int last_clone_in_deme = clones_list_in_deme[deme_index][num_clones_in_deme - 1];
	int index_in_deme_of_last_clone_in_deme = clone_ints[INDEX_IN_DEME][last_clone_in_deme];

	int deme_of_last_clone = clone_ints[DEME][*num_clones - 1];
	int index_in_deme_of_last_clone;

	if(index_in_deme < num_clones_in_deme - 1) {
		// replace the clone in the list of clones in the deme, unless it's at the end of the list:
		clones_list_in_deme[deme_index][index_in_deme] = clones_list_in_deme[deme_index][num_clones_in_deme - 1];
		
		if(use_clone_bintrees) {
			clones_pops_in_deme[deme_index][index_in_deme] = clones_pops_in_deme[deme_index][num_clones_in_deme - 1];
			clones_rates_in_deme[deme_index][index_in_deme] = clones_rates_in_deme[deme_index][num_clones_in_deme - 1];
		}
		
		// update the list of clones in the deme in light of the replacement:
		clone_ints[INDEX_IN_DEME][last_clone_in_deme] = index_in_deme;

		// update bintree in light of the replacement:
		if(use_clone_bintrees && index_in_deme / SET_SIZE != index_in_deme_of_last_clone_in_deme / SET_SIZE) {
			update_bintree_layers_int(clone_ints[POPULATION][last_clone_in_deme], bintree_clone_ints[deme_index], index_in_deme, max_layer_needed, max_clones_per_deme);
			update_bintree_layers_int(-clone_ints[POPULATION][last_clone_in_deme], bintree_clone_ints[deme_index], index_in_deme_of_last_clone_in_deme, max_layer_needed, max_clones_per_deme);
			
			update_bintree_layers(clones_rates_in_deme[deme_index][num_clones_in_deme - 1], bintree_clone_doubles[deme_index], index_in_deme, max_layer_needed, max_clones_per_deme);
			update_bintree_layers(-clones_rates_in_deme[deme_index][num_clones_in_deme - 1], bintree_clone_doubles[deme_index], index_in_deme_of_last_clone_in_deme, max_layer_needed, max_clones_per_deme);
		}
	}
	deme_ints[NUM_CLONES_IN_DEME][deme_index]--;

	if(chosen_clone < *num_clones - 1) {
		// replace the clone in clones array, unless it's at the end of the array:
		for(i = 0; i < NUM_CLONE_INT_PROPS; i++) clone_ints[i][chosen_clone] = clone_ints[i][*num_clones - 1];

		// update the list of clones in the deme that contained the clone at the end of the array:
		index_in_deme_of_last_clone = clone_ints[INDEX_IN_DEME][*num_clones - 1];
		clones_list_in_deme[deme_of_last_clone][index_in_deme_of_last_clone] = chosen_clone;
	}
	*num_clones = *num_clones - 1; // decrement the count of clones
}

///// deme:

// move cells from the dividing deme into the newly created deme (or remove cells if deme is dividing beyond the edge):
void move_cells(long *idum, int origin_deme_num, int dividing_beyond_the_edge, int new_deme_index, int *num_cells, int *num_demes, int *num_empty_demes, int *num_clones, 
	int *num_empty_cols, int *num_empty_driver_cols, int *num_extinct_genotypes, int *num_extinct_driver_genotypes, int *event_counter)
{
	int deme_pop_all_cells = deme_ints[POPULATION][origin_deme_num] + deme_ints[NORMAL_CELLS][origin_deme_num];
	int num_already_sampled_from_deme = 0;
	int previous_clone_pops = 0;
	int new_deme_pop_all_cells = stochastic_round((float)deme_pop_all_cells/2, idum); // number of cells to be moved to the new deme
	int i;
	int clone_num, clone_pop, geno_num, driver_geno_num;
	int num_draws, clone_sample_size;

	for(i = -1; i < deme_ints[NUM_CLONES_IN_DEME][origin_deme_num]; i++) { // loop over all clones and the normal cells in the origin deme
		if(i >= 0) {
			clone_num = clones_list_in_deme[origin_deme_num][i];
			clone_pop = clone_ints[POPULATION][clone_num];
			geno_num = clone_ints[GENOTYPE][clone_num];
			driver_geno_num = clone_ints[DRIVER_GENOTYPE][clone_num];
		}
		else clone_pop = deme_ints[NORMAL_CELLS][origin_deme_num];
		// clone_sample_size is drawn from a hypergeometric distribution (i.e. sampling without replacement):
		num_draws = MAX(new_deme_pop_all_cells - num_already_sampled_from_deme, 0);
		clone_sample_size = hypergeometric(idum, clone_pop, deme_pop_all_cells - previous_clone_pops - clone_pop, num_draws);

		num_already_sampled_from_deme += clone_sample_size;
		previous_clone_pops += clone_pop;

		if(clone_sample_size > 0) {
			if(i >= 0) { // for a clone:
				increment_or_decrement_deme(-clone_sample_size, origin_deme_num, *num_cells, *num_demes, num_empty_demes); // remove cells from origin deme
				if(!dividing_beyond_the_edge) increment_or_decrement_deme(+clone_sample_size, new_deme_index, *num_cells, *num_demes, num_empty_demes); // add cells to new deme

				if(!dividing_beyond_the_edge) create_clone(geno_num, num_clones, new_deme_index, driver_geno_num, *num_demes, clone_sample_size); // create clone in new deme
				increment_or_decrement_clone(clone_num, origin_deme_num, num_clones, -clone_sample_size, *num_demes); // remove cells from clone in origin deme

				if(dividing_beyond_the_edge) {
					*num_cells = *num_cells - clone_sample_size;
					increment_or_decrement_genotype(genotype_ints, geno_num, empty_cols, num_empty_cols, -clone_sample_size, num_extinct_genotypes);
					increment_or_decrement_genotype(driver_genotype_ints, driver_geno_num, empty_driver_cols, num_empty_driver_cols, -clone_sample_size, num_extinct_driver_genotypes);
				}
			}
			else { // for normal cells:
				add_or_remove_normal_cell(-clone_sample_size, origin_deme_num, *num_cells, event_counter, *num_demes); // remove normal cells from origin deme
				if(!dividing_beyond_the_edge) add_or_remove_normal_cell(+clone_sample_size, new_deme_index, *num_cells, event_counter, *num_demes); // add normal cells to new deme
			}
		}
	}
}

// budge demes, starting at the far end of the line:
void budge_demes(int old_x, int old_y, int *x_to_fill, int *y_to_fill)
{
	int x_to_move, y_to_move;
	float xdiff, ydiff;
	float hypoten;

	do {
		xdiff = old_x - *x_to_fill; // difference between coordinate of dividing deme and coordinate of the site to be filled
		ydiff = old_y - *y_to_fill;
		hypoten = sqrt(xdiff * xdiff + ydiff * ydiff); // Euclidean distance between sites

		// coordinates of the site to be budged into the site to be filled:
		if(fabs(xdiff/hypoten) >= ROOT2-1) x_to_move = *x_to_fill + SIGN(xdiff);
		else x_to_move = *x_to_fill;
		if(fabs(ydiff/hypoten) >= ROOT2-1) y_to_move = *y_to_fill + SIGN(ydiff);
		else y_to_move = *y_to_fill;

		if(x_to_move == old_x && y_to_move == old_y) break; // no more budging to be done
		
		if(grid[x_to_move][y_to_move] != EMPTY) { // update the deme's coordinates:
			deme_ints[XCOORD][grid[x_to_move][y_to_move]] = *x_to_fill;
			deme_ints[YCOORD][grid[x_to_move][y_to_move]] = *y_to_fill;
		}
		grid[*x_to_fill][*y_to_fill] = grid[x_to_move][y_to_move]; // make the grid square point to a different deme index

		*x_to_fill = x_to_move;
		*y_to_fill = y_to_move;
	} while(1);
}

// get coordinates of first deme to be filled during deme fission:
void get_deme_coordinates(int *x_to_fill, int *y_to_fill, int old_x, int old_y, long *idum)
{
	float theta, tan_theta;
	int direction, quadrant;
	int x_add = 0, y_add = 0;

	// if demes can be budged then the first site to be filled is on the edge of the grid
	if(!migration_edge_only) {
		theta = ran1(idum) * 2 * PI; // random angle along which demes will be budged
		tan_theta = tan(theta);
		quadrant = which_quadrant(old_x, old_y, theta, tan_theta, dim_grid); // quadrant corresponding to angle theta
		
		if(quadrant == 1) { // left quadrant
			*x_to_fill = 0; // coordinates of site that will receive a cell through budging or proliferation
			*y_to_fill = old_y - (float)old_x / tan_theta;
		}
		else if(quadrant == 2) { // top quadrant
			*x_to_fill = old_x + tan_theta * (dim_grid - old_y); // coordinates of site that will receive a cell through budging or proliferation
			*y_to_fill = dim_grid - 1;
		}
		else if(quadrant == 3) { // right quadrant
			*x_to_fill = dim_grid - 1; // coordinates of site that will receive a cell through budging or proliferation
			*y_to_fill = old_y + (float)(dim_grid - old_x) / tan_theta;
		}
		else { // bottom quadrant
			*x_to_fill = old_x - tan_theta * old_y; // coordinates of site that will receive a cell through budging or proliferation
			*y_to_fill = 0;
		}
	}
	// if budging is not allowed then the site to be filled is one of the four nearest neighbours
	else {
		direction = floor(4 * ran1(idum)); // randomly choose a direction;
		// directions are right, left, up, down:
		if(direction == 0) x_add = 1;
		else if(direction == 1) x_add = -1;
		else if(direction == 2) y_add = 1;
		else y_add = -1;

		*x_to_fill = old_x + x_add; // coordinates of site that will receive a cell through budging or proliferation
		*y_to_fill = old_y + y_add;
	}
}

// choose the cooridinates of a destination deme:
void choose_grid_square(int old_deme_index, long *idum, int *new_x, int *new_y)
{
	int old_x = deme_ints[XCOORD][old_deme_index];
	int old_y = deme_ints[YCOORD][old_deme_index];
	int direction = floor(4 * ran1(idum)); // randomly choose a direction

	// directions are right, left, up, down:
	if(direction == 0) {*new_x = old_x + 1; *new_y = old_y;}
	else if(direction == 1) {*new_x = old_x - 1; *new_y = old_y;}
	else if(direction == 2) {*new_x = old_x; *new_y = old_y + 1;}
	else {*new_x = old_x; *new_y = old_y - 1;}
}

// initialise a deme before it becomes occupied for the first time:
void create_deme(int new_x, int new_y, int *num_demes, int num_cells)
{
	int max_layer_needed;
	int new_deme_index = *num_demes; // index of the new deme
	
	grid[new_x][new_y] = new_deme_index; // point from grid coordinates to the new deme index
	
	*num_demes = *num_demes + 1; // increment count of demes

	// initialise deme properties (before any cancer cell arrives):
	deme_ints[POPULATION][new_deme_index] = 0;
	deme_ints[XCOORD][new_deme_index] = new_x;
	deme_ints[YCOORD][new_deme_index] = new_y;
	if(normal_birth_rate < 0 || migration_type >= 2) deme_ints[NORMAL_CELLS][new_deme_index] = 0; // if normal_birth_rate < 0 or if deme is created via fission then initial normal cell count = 0
	else deme_ints[NORMAL_CELLS][new_deme_index] = K; // otherwise initial normal cell count = K
	deme_ints[NUM_CLONES_IN_DEME][new_deme_index] = 0;
	deme_floats[DEATH_RATE][new_deme_index] = set_death_rate(new_deme_index, num_cells);
	deme_floats[MIGRATION_MODIFIER][new_deme_index] = set_migration_modifier(deme_ints[NORMAL_CELLS][new_deme_index], deme_ints[POPULATION][new_deme_index]);
	deme_doubles[SUM_BIRTH_RATES][new_deme_index] = 0;
	deme_doubles[SUM_MIGRATION_RATES][new_deme_index] = 0;
	deme_doubles[SUM_RATES][new_deme_index] = (double)deme_ints[NORMAL_CELLS][new_deme_index] * ((double)normal_birth_rate + (double)deme_floats[DEATH_RATE][new_deme_index]);

	max_layer_needed = max_layer_needed_array[*num_demes];

	extend_bintree(bintree_deme_doubles, new_deme_index, max_layer_needed, max_demes);

	update_bintree_layers(deme_doubles[SUM_RATES][new_deme_index], bintree_deme_doubles, new_deme_index, max_layer_needed, max_demes);
}

// remove a deme and all of its contents:
void remove_deme(int deme_index, int *num_cells, int *num_clones, int *num_demes, int *num_empty_demes, int *empty_cols, int *num_empty_cols, int *num_empty_driver_cols, int *empty_driver_cols, 
	int *num_extinct_genotypes, int *num_extinct_driver_genotypes)
{
	int i, clone_i;
	int clone_num, geno_num, driver_geno_num;
	int deme_of_last_clone;
	int max_layer_needed = max_layer_needed_array[*num_demes];

	*num_cells = *num_cells - deme_ints[POPULATION][deme_index];

	// remove clones and decrement genotypes and driver genotypes:
	for(clone_i = 0; clone_i < deme_ints[NUM_CLONES_IN_DEME][deme_index]; clone_i++) {
		clone_num = clones_list_in_deme[deme_index][clone_i];
		geno_num = clone_ints[GENOTYPE][clone_num];
		driver_geno_num = clone_ints[DRIVER_GENOTYPE][clone_num];

		// subtract clone population from genotype and driver genotype populations:
		increment_or_decrement_genotype(genotype_ints, geno_num, empty_cols, num_empty_cols, -clone_ints[POPULATION][clone_num], num_extinct_genotypes);
		increment_or_decrement_genotype(driver_genotype_ints, driver_geno_num, empty_driver_cols, num_empty_driver_cols, -clone_ints[POPULATION][clone_num], num_extinct_driver_genotypes);

		// replace the clone in clones array, and in list of clones in relevant deme, unless it's at the end of the array:
		if(clone_num < *num_clones - 1) {
			// copy all properties from final clone to clone being replaced:
			for(i=0; i<NUM_CLONE_INT_PROPS; i++) clone_ints[i][clone_num] = clone_ints[i][*num_clones - 1];

			// update the list of clones in the deme that contained the clone that has changed index:
			deme_of_last_clone = clone_ints[DEME][*num_clones - 1];
			for(i = 0; i < deme_ints[NUM_CLONES_IN_DEME][deme_of_last_clone]; i++) {
				if(clones_list_in_deme[deme_of_last_clone][i] == *num_clones - 1) {
					clones_list_in_deme[deme_of_last_clone][i] = clone_num;
					break;
				}
			}
		}
		*num_clones = *num_clones - 1; // decrement the count of clones
	}

	// if the deme cancer cell count was zero then remove it from the count of empty demes:
	if(deme_ints[POPULATION][deme_index] == 0) *num_empty_demes = *num_empty_demes - 1;

	// record the site as being empty:
	grid[deme_ints[XCOORD][deme_index]][deme_ints[YCOORD][deme_index]] = EMPTY;

	// update deme and bintree rates:
	update_bintree_layers(-deme_doubles[SUM_RATES][deme_index], bintree_deme_doubles, deme_index, max_layer_needed, max_demes);

	// replace the deme in each demes array and update the deme number of the clones, unless the deme is at the end of the array:
	if(deme_index < *num_demes - 1) {
		for(i=0; i<deme_ints[NUM_CLONES_IN_DEME][*num_demes - 1]; i++) {
			clones_list_in_deme[deme_index][i] = clones_list_in_deme[*num_demes - 1][i];
			clone_ints[DEME][clones_list_in_deme[*num_demes - 1][i]] = deme_index;
		}
		
		grid[deme_ints[XCOORD][*num_demes - 1]][deme_ints[YCOORD][*num_demes - 1]] = deme_index;

		for(i=0; i<NUM_DEME_INT_PROPS; i++) deme_ints[i][deme_index] = deme_ints[i][*num_demes - 1];
		for(i=0; i<NUM_DEME_FLOAT_PROPS; i++) deme_floats[i][deme_index] = deme_floats[i][*num_demes - 1];
		for(i=0; i<NUM_DEME_DOUBLE_PROPS; i++) deme_doubles[i][deme_index] = deme_doubles[i][*num_demes - 1];
	}
	*num_demes = *num_demes - 1;
}

// increment or decrement a deme:
void increment_or_decrement_deme(int change, int deme_index, int num_cells, int num_demes, int *num_empty_demes)
{
	// add or subtract from deme cancer cell count:
	deme_ints[POPULATION][deme_index] += change;

	// update deme and bintree rates:
	update_deme_bintree_environmental_rates(change, 0, num_demes, deme_index, num_cells);

	// if the deme cancer cell count has reached zero then add to the count of empty demes:
	if(deme_ints[POPULATION][deme_index] == 0) *num_empty_demes = *num_empty_demes + 1;
}

// increment or decrement a deme's normal cell count:
void add_or_remove_normal_cell(int change, int deme_index, int num_cells, int *event_counter, int num_demes)
{
	if(change > 0) event_counter[NORMALBIRTH_EVENT]++;
	else event_counter[NORMALDEATH_EVENT]++;

	// add or subtract from deme normal cell count:
	deme_ints[NORMAL_CELLS][deme_index] += change;

	// update deme and bintree rates:
	update_deme_bintree_environmental_rates(0, change, num_demes, deme_index, num_cells);

	// error check:
	check_normal_pops(deme_index);
}

/////////////////////// memory handling:

// assign memory:
void assign_memory()
{
	printf("max_clones %d; max_genotypes %d; max_driver_genotypes %d; max_demes %d\n", max_clones, max_genotypes, max_driver_genotypes, max_demes);
	printf("dim_grid %d; matrix_max %d; max_distinct_allele_freqs %d\n", dim_grid, matrix_max, max_distinct_allele_freqs);
	printf("max_clones_per_deme %d; max_bintree_clone_elements_per_deme %d\n", max_clones_per_deme, max_bintree_clone_elements_per_deme);
	
	// two-dim arrays:
	mallocArray_int(&clone_ints, NUM_CLONE_INT_PROPS, max_clones);
	mallocArray_int(&genotype_ints, NUM_GENOTYPE_INT_PROPS, max_genotypes);
	mallocArray_int(&genotype_relatives, 4, max_genotypes);
	mallocArray_float(&genotype_floats, NUM_GENOTYPE_FLOAT_PROPS, max_genotypes);
	mallocArray_int(&driver_genotype_ints, NUM_GENOTYPE_INT_PROPS, max_driver_genotypes);
	mallocArray_float(&driver_genotype_floats, NUM_GENOTYPE_FLOAT_PROPS, max_driver_genotypes);
	mallocArray_int(&deme_ints, NUM_DEME_INT_PROPS, max_demes);
	mallocArray_float(&deme_floats, NUM_DEME_FLOAT_PROPS, max_demes);
	mallocArray_double(&deme_doubles, NUM_DEME_DOUBLE_PROPS, max_demes);
	mallocArray_int(&grid, dim_grid, dim_grid);
	if(matrix_max > 0) mallocArray_int(&driver_matrix, matrix_max, matrix_max);
	else mallocArray_int(&driver_matrix, 1, 1); // if driver matrix isn't needed then allocate it minimal memory
	mallocArray_int(&freq_table, 2, max_distinct_allele_freqs);
	mallocArray_int(&clones_list_in_deme, max_demes, max_clones_per_deme);
	mallocArray_float(&depth_diversity, 3, 12); // first dimension: 1, 2 or 3 samples; second dimension: depths 0-10 and random sampling
	mallocArray_float(&depth_diversity_bigsample, 3, 12); // first dimension: 1, 2 or 3 samples; second dimension: depths 0-10 and random sampling
	if(record_matrix) mallocArray_int(&matrix, matrix_max, matrix_max);
	else mallocArray_int(&matrix, 1, 1); // if matrix isn't needed then allocate it minimal memory
	if(use_clone_bintrees) {
		mallocArray_double(&bintree_clone_doubles, max_demes, max_bintree_clone_elements_per_deme);
		mallocArray_int(&bintree_clone_ints, max_demes, max_bintree_clone_elements_per_deme);
		mallocArray_int(&clones_pops_in_deme, max_demes, max_clones_per_deme);
		mallocArray_double(&clones_rates_in_deme, max_demes, max_clones_per_deme);
	}
	else { // if clone bintrees aren't needed then allocate them minimal memory
		mallocArray_double(&bintree_clone_doubles, 1, 1);
		mallocArray_int(&bintree_clone_ints, 1, 1);
		mallocArray_int(&clones_pops_in_deme, 1, 1);
		mallocArray_double(&clones_rates_in_deme, 1, 1);
	}

	// one-dim arrays:
	empty_cols = (int *) malloc(max_genotypes * sizeof *empty_cols);
	if(empty_cols == NULL) {printf("Memory problem: empty_cols\n"); exit(1);}
	empty_driver_cols = (int *) malloc(max_driver_genotypes * sizeof *empty_driver_cols);
	if(empty_driver_cols == NULL) {printf("Memory problem: empty_driver_cols\n"); exit(1);}
	within_deme_diversity = (float *) malloc(max_demes * sizeof *within_deme_diversity);
	if(within_deme_diversity == NULL) {printf("Memory problem: within_deme_diversity\n"); exit(1);}
	within_deme_driver_diversity = (float *) malloc(max_demes * sizeof *within_deme_driver_diversity);
	if(within_deme_driver_diversity == NULL) {printf("Memory problem: within_deme_driver_diversity\n"); exit(1);}
	allele_count = (int *) malloc(max_genotypes * sizeof *allele_count);
	if(allele_count == NULL) {printf("Memory problem: allele_count\n"); exit(1);}
	bintree_deme_doubles = (double *) malloc(max_bintree_deme_elements * sizeof *bintree_deme_doubles);
	if(bintree_deme_doubles == NULL) {printf("Memory problem: bintree_deme_doubles\n"); exit(1);}
	position_in_edge_list = (int *) malloc(max_demes * sizeof *position_in_edge_list);
	if(position_in_edge_list == NULL) {printf("Memory problem: position_in_edge_list\n"); exit(1);}
	demes_at_edge = (int *) malloc(max_demes * sizeof *demes_at_edge);
	if(demes_at_edge == NULL) {printf("Memory problem: demes_at_edge\n"); exit(1);}
	sides_at_edge = (int *) malloc(max_demes * sizeof *sides_at_edge);
	if(sides_at_edge == NULL) {printf("Memory problem: sides_at_edge\n"); exit(1);}
	genotype_edge_pop = (int *) malloc(max_genotypes * sizeof *genotype_edge_pop);
	if(genotype_edge_pop == NULL) {printf("Memory problem: genotype_edge_pop\n"); exit(1);}
	at_edge = (int *) malloc(max_demes * sizeof *at_edge);
	if(at_edge == NULL) {printf("Memory problem: at_edge\n"); exit(1);}
	buff_array = (float *) malloc(max_clones_per_deme * sizeof *buff_array);
	if(buff_array == NULL) {printf("Memory problem: buff_array\n"); exit(1);}
	max_layer_needed_array = (int *) malloc(length_of_max_layer_needed_array * sizeof *max_layer_needed_array);
	if(max_layer_needed_array == NULL) {printf("Memory problem: max_layer_needed_array\n"); exit(1);}
	min_deme_bintree_index_by_layer = (int *) malloc((get_max_layer_needed(max_demes) + 1) * sizeof *min_deme_bintree_index_by_layer);
	if(min_deme_bintree_index_by_layer == NULL) {printf("Memory problem: min_deme_bintree_index_by_layer\n"); exit(1);}
	min_clone_bintree_index_by_layer = (int *) malloc((get_max_layer_needed(max_clones_per_deme) + 1) * sizeof *min_clone_bintree_index_by_layer);
	if(min_clone_bintree_index_by_layer == NULL) {printf("Memory problem: min_clone_bintree_index_by_layer\n"); exit(1);}
	buff_small_double = (double *) malloc(MAX(SET_SIZE, 3) * sizeof *buff_small_double);
	if(buff_small_double == NULL) {printf("Memory problem: buff_small_double\n"); exit(1);}
	buff_small_float = (float *) malloc(MAX(SET_SIZE, 3) * sizeof *buff_small_float);
	if(buff_small_float == NULL) {printf("Memory problem: buff_small_float\n"); exit(1);}
	buff_array_int = (int *) malloc(MAX(MAX(SET_SIZE, 3), max_genotypes) * sizeof *buff_array_int);
	if(buff_array_int == NULL) {printf("Memory problem: buff_array_int\n"); exit(1);}

	printf("Assigned memory\n");
}

// free memory:
void free_memory()
{
	freeArray_int(clone_ints, NUM_CLONE_INT_PROPS);
	freeArray_float(genotype_floats, NUM_GENOTYPE_FLOAT_PROPS);
	freeArray_float(driver_genotype_floats, NUM_GENOTYPE_FLOAT_PROPS);
	freeArray_int(deme_ints, NUM_DEME_INT_PROPS);
	freeArray_float(deme_floats, NUM_DEME_FLOAT_PROPS);
	freeArray_double(deme_doubles, NUM_DEME_DOUBLE_PROPS);
	freeArray_int(grid, dim_grid);
	if(record_matrix) freeArray_int(matrix, matrix_max);
	else freeArray_int(matrix, 1);
	if(matrix_max > 0) freeArray_int(driver_matrix, matrix_max);
	else freeArray_int(driver_matrix, 1);
	free(empty_cols);
	free(empty_driver_cols);
	free(within_deme_diversity);
	free(within_deme_driver_diversity);
	free(bintree_deme_doubles);
	freeArray_int(clones_list_in_deme, max_demes);
	freeArray_float(depth_diversity, 3);
	freeArray_float(depth_diversity_bigsample, 3);
	freeArray_int(genotype_ints, NUM_GENOTYPE_INT_PROPS);
	freeArray_int(driver_genotype_ints, NUM_GENOTYPE_INT_PROPS);
	free(allele_count);
	freeArray_int(freq_table, 2);
	freeArray_int(genotype_relatives, 4);
	free(position_in_edge_list);
	free(demes_at_edge);
	free(sides_at_edge);
	free(genotype_edge_pop);
	free(at_edge);
	free(buff_array);
	free(max_layer_needed_array);
	free(buff_small_float);
	free(buff_small_double);
	free(buff_array_int);
	free(min_deme_bintree_index_by_layer);
	free(min_clone_bintree_index_by_layer);
	if(use_clone_bintrees) {
		freeArray_double(bintree_clone_doubles, max_demes);
		freeArray_int(bintree_clone_ints, max_demes);
		freeArray_int(clones_pops_in_deme, max_demes);
		freeArray_double(clones_rates_in_deme, max_demes);
	}
	else {
		freeArray_double(bintree_clone_doubles, 1);
		freeArray_int(bintree_clone_ints, 1);
		freeArray_int(clones_pops_in_deme, 1);
		freeArray_double(clones_rates_in_deme, 1);
	}
}

/////////////////////// file handling:

// open output files:
void open_files(char *input_and_output_path)
{
	char filebuff[200];

	// create or open output files:
	sprintf(filebuff, "%sparameters.dat", input_and_output_path);
	output_parameters=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput.dat", input_and_output_path);
	output_pops=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_diversities.dat", input_and_output_path);
	output_diversities=fopen(filebuff, "w+");
	sprintf(filebuff, "%sdemes.dat", input_and_output_path);
	output_demes=fopen(filebuff, "w+");
	sprintf(filebuff, "%sclones.dat", input_and_output_path);
	output_clones=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_genotype_counts.dat", input_and_output_path);
	output_genotype_counts=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_driver_genotype_counts.dat", input_and_output_path);
	output_driver_genotype_counts=fopen(filebuff, "w+");
	sprintf(filebuff, "%sphylo.dat", input_and_output_path);
	output_phylo=fopen(filebuff, "w+");
	sprintf(filebuff, "%sdriver_phylo.dat", input_and_output_path);
	output_driver_phylo=fopen(filebuff, "w+");
	sprintf(filebuff, "%smatrix.dat", input_and_output_path);
	output_matrix=fopen(filebuff, "w+");
	sprintf(filebuff, "%sdriver_matrix.dat", input_and_output_path);
	output_driver_matrix=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_popgrid.dat", input_and_output_path);
	output_popgrid=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_passengersgrid.dat", input_and_output_path);
	output_passengersgrid=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_normalcellsgrid.dat", input_and_output_path);
	output_normalcellsgrid=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_deathrates_grid.dat", input_and_output_path);
	output_deathrates_grid=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_birthratesgrid.dat", input_and_output_path);
	output_birthratesgrid=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_migrationratesgrid.dat", input_and_output_path);
	output_migrationratesgrid=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_driversgrid.dat", input_and_output_path);
	output_driversgrid=fopen(filebuff, "w+");
	sprintf(filebuff, "%serror_log.dat", input_and_output_path);
	error_log=fopen(filebuff, "w+");
	sprintf(filebuff, "%ssample_size_log.dat", input_and_output_path);
	sample_size_log=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_allele_counts.dat", input_and_output_path);
	output_allele_counts=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_driver_allele_counts.dat", input_and_output_path);
	output_driver_allele_counts=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_genotype_properties.dat", input_and_output_path);
	output_genotype_properties=fopen(filebuff, "w+");
	sprintf(filebuff, "%soutput_driver_genotype_properties.dat", input_and_output_path);
	output_driver_genotype_properties=fopen(filebuff, "w+");

	// check that files were successfully opened:
	if (output_parameters == NULL || output_pops == NULL || output_diversities == NULL || output_demes == NULL || output_clones == NULL || output_genotype_counts == NULL || output_driver_genotype_counts == NULL || 
		output_phylo == NULL || output_driver_phylo == NULL || output_matrix == NULL || output_driver_matrix == NULL || output_popgrid == NULL || output_passengersgrid == NULL || output_normalcellsgrid == NULL || 
		output_deathrates_grid == NULL || output_birthratesgrid == NULL || output_migrationratesgrid == NULL || output_driversgrid == NULL || error_log == NULL || sample_size_log == NULL || 
		output_allele_counts == NULL || output_driver_allele_counts == NULL || output_genotype_properties == NULL || output_driver_genotype_properties == NULL) {
		perror("Failed to open output file");
		if(error_log != NULL) fprintf(error_log, "Failed to open output file\n");
		exit(1);
	}

	// open gnuplot pipe:
	gp = popen("gnuplot -persist", "w"); // 'gp' is the pipe descriptor
	if (gp == NULL) {
		printf("Error opening pipe to gnuplot.\n");
		fprintf(error_log, "Error opening pipe to gnuplot.\n");
		write_grid = 0;
	}
}

// write column headings and parameter values to files:
void initiate_files(int *num_samples_list)
{
	int i, j;

	// write parameter values to file:
	fprintf(output_parameters, "K\tmigration_type\tinit_migration_rate\tmigration_edge_only\tmigration_rate_scales_with_K\tmu_driver_birth\tmu_passenger\tmu_driver_migration\t");
	fprintf(output_parameters, "normal_birth_rate\tbaseline_death_rate\ts_driver_birth\ts_passenger\ts_driver_migration\tmax_relative_birth_rate\tmax_relative_migration_rate\tinit_pop\tseed\tbiopsy_size_per_sample\n");
	fprintf(output_parameters, "%d\t%d\t%e\t%d\t%d\t%e\t%e\t%e\t", K, migration_type, init_migration_rate, migration_edge_only, migration_rate_scales_with_K, 
		mu_driver_birth, mu_passenger, mu_driver_migration);
	fprintf(output_parameters, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\n", normal_birth_rate, baseline_death_rate, s_driver_birth, s_passenger, s_driver_migration, max_relative_birth_rate, max_relative_migration_rate, 
		init_pop, seed, biopsy_size_per_sample);
	fflush(output_parameters);

	// write column names to output_pops file:
	fprintf(output_pops, "Generation\tNumCells\tTotalBirths\tTotalDeaths\tTotalMigrations\tPassengers\tDrivers\tVarPassengers\tVarDrivers\t");
	fprintf(output_pops, "MeanBirthRate\tMeanMigrationRate\tVarBirthRate\tVarMigRate\tNumClones\tNumGenotypes\tNumDriverGenotypes\tNumDemes\tNumNon-EmptyDemes\tMatrixColumns\tDriverMatrixColumns\tEverGenotypes");
	for(i=0; i<=MAX_DRIVERS_TO_COUNT; i++) fprintf(output_pops, "\tCellsWith%dDrivers", i);
	fprintf(output_pops, "\n");

	fprintf(output_diversities, "Generation\tNumCells");
	if(calculate_total_diversity) {
		fprintf(output_diversities, "\tDiversity\tEdgeDiversity\tAlphaDiversity");
	}
	fprintf(output_diversities, "\tDriverDiversity\tDriverEdgeDiversity\tDriverAlphaDiversity");
	if(!filled_grid && biopsy_size_per_sample > 0) for(i=0; i<3; i++) {
		if(!well_mixed) for(j=0; j<=10; j++) fprintf(output_diversities, "\tDriverDiversityFrom%dSamplesAtDepth%d", num_samples_list[i], j);
		fprintf(output_diversities, "\tDriverDiversityFrom%dRandomSamples", num_samples_list[i]);
		if(!well_mixed) for(j=0; j<=10; j++) fprintf(output_diversities, "\tDriverDiversityFrom%dBigSamplesAtDepth%d", num_samples_list[i], j);
		fprintf(output_diversities, "\tDriverDiversityFrom%dBigRandomSamples", num_samples_list[i]);
	}
	fprintf(output_diversities, "\n");
	
	// write column names to other output files:
	fprintf(output_demes, "Generation\tX\tY\tPopulation\tNormalCells\tDeathRate\tDiversity\tDriverDiversity\n");
	fprintf(output_clones, "Generation\tClone\tDeme\tGenotype\tDriverGenotype\tParent\tDriverParent\tX\tY\tNormalCells\tPopulation\tBirthRate\t");
	fprintf(output_clones, "MigrationRate\tDeathRate\tMigrationModifier\tDriverMutations\tMigrationMutations\tPassengerMutations\n");
	fprintf(output_phylo, "NumCells\tNumSamples\tCellsPerSample\tSampleDepth\tGeneration\tIdentity\tParent\tPopulation\tBirthRate\tMigrationRate\tOriginTime\tDriverMutations\tMigrationMutations\tPassengerMutations\n");
	fprintf(output_driver_phylo, "NumCells\tNumSamples\tCellsPerSample\tSampleDepth\tGeneration\tIdentity\tDriverIdentity\tParent\tPopulation\tBirthRate\tMigrationRate\tOriginTime\tDriverMutations\tMigrationMutations\tPassengerMutations\n");
	fprintf(sample_size_log, "Generation\tNumSamples\tDepth\tTargetSampleSize\tActualSampleSize\n");
	fprintf(output_allele_counts, "Generation\tSize\tFrequency\tCount\n");
	fprintf(output_driver_allele_counts, "Generation\tSize\tFrequency\tCount\n");
	fprintf(output_genotype_counts, "Generation\tSize\tFrequency\tCount\n");
	fprintf(output_driver_genotype_counts, "Generation\tSize\tFrequency\tCount\n");
	fprintf(output_genotype_properties, "Population\tParent\tIdentity\tDriverIdentity\tDriverMutations\tMigrationMutations\tImmortal\tPassengerMutations\tBirthRate\tMigrationRate\tOriginTime\tDescendants\n");
	fprintf(output_driver_genotype_properties, "Population\tParent\tIdentity\tDriverIdentity\tDriverMutations\tMigrationMutations\tImmortal\tPassengerMutations\tBirthRate\tMigrationRate\tOriginTime\tDescendants\n");
}

// calculate metrics and write to files:
void main_calculations_and_output(long *idum, int num_demes, int num_matrix_cols, int num_driver_matrix_cols, float gens_elapsed, int *driver_counts, int num_cells, int num_clones, 
	int record_phylogenies, int *num_samples_list, int next_genotype_id, int next_driver_genotype_id, int *event_counter, int num_empty_cols, int num_empty_driver_cols, int num_extinct_genotypes, 
	int num_extinct_driver_genotypes, int num_empty_demes, long t1, int print_output, int which_parts)
{
	int i;
	int rao;
	float centre_X, centre_Y;
	int num_freqs;
	double sum_birth_rates, sum_death_rates, sum_migration_rates; // sums of rates of all cancer cells
	double sum_normal_birth_rates, sum_normal_death_rates; // sums of rates of all normal cells
	float mean_num_passengers, mean_num_drivers, var_num_passengers, var_num_drivers;
	float variance_mig_rate, variance_birth_rate, diversity, alpha_diversity, edge_diversity, driver_diversity, alpha_driver_diversity, edge_driver_diversity;

	if(matrix_max > 0) rao = 1;
	else rao = 0;

	// calculate sums of rates for all cancer cells and for all normal cells:
	calculate_sums_of_rates(&sum_death_rates, &sum_birth_rates, &sum_migration_rates, &sum_normal_birth_rates, &sum_normal_death_rates,
		num_demes, num_matrix_cols, gens_elapsed);

	// calculate means and variances of numbers of mutations per cell:
	calculate_mutation_metrics(&mean_num_passengers, &mean_num_drivers, &var_num_passengers, &var_num_drivers, driver_counts, num_matrix_cols, num_cells);

	// calculate variance of cell birth and migration rates:
	variance_birth_rate = calculate_variance_of_rate(num_matrix_cols, num_cells, sum_birth_rates, genotype_floats[BIRTH_RATE]);
	variance_mig_rate = calculate_variance_of_rate(num_matrix_cols, num_cells, sum_migration_rates, genotype_floats[MIGRATION_RATE]);

	if(component(which_parts, PHYLO_AND_POPS)) {
		// write phylogenetic data of driver genotypes to file:
		write_output_phylo(output_driver_phylo, num_driver_matrix_cols, gens_elapsed, driver_genotype_ints[POPULATION], driver_genotype_ints, driver_genotype_floats, 1, -1, -1, num_cells);

		// write population metrics to file:
		write_output_pops(output_pops, num_cells, event_counter, mean_num_passengers, mean_num_drivers, sum_birth_rates, sum_migration_rates, var_num_passengers, var_num_drivers, variance_birth_rate, 
			variance_mig_rate, num_matrix_cols, num_empty_cols, num_driver_matrix_cols, num_empty_driver_cols, num_extinct_genotypes, num_extinct_driver_genotypes, num_demes, num_empty_demes, gens_elapsed, 
			num_clones, driver_counts, next_genotype_id);
	}

	if(component(which_parts, DIVERSITIES)) {
		// calculate diversity metrics of driver genotypes:
		get_diversity_metrics(rao, &driver_diversity, &edge_driver_diversity, &alpha_driver_diversity, within_deme_driver_diversity, num_demes, num_cells, num_driver_matrix_cols, num_clones, 
			clone_ints[DRIVER_GENOTYPE], driver_matrix, dmax, driver_genotype_ints[POPULATION], idum, gens_elapsed, position_in_edge_list, demes_at_edge, sides_at_edge, genotype_edge_pop, 1);
		// calculate diversity metrics of genotypes:
		if(calculate_total_diversity) get_diversity_metrics(0, &diversity, &edge_diversity, &alpha_diversity, within_deme_diversity, num_demes, num_cells, num_matrix_cols, num_clones, 
			clone_ints[GENOTYPE], matrix, dmax, genotype_ints[POPULATION], idum, gens_elapsed, position_in_edge_list, demes_at_edge, sides_at_edge, genotype_edge_pop, 0);

		if(!filled_grid && biopsy_size_per_sample > 0) {
			// find the tumour's centre of gravity:
			centre_X = centre_of_gravity(1, num_cells);
			centre_Y = centre_of_gravity(2, num_cells);

			// calculate driver diversity and record driver phylogenies of sampled cells:
			for(i = 0; i < 3; i++) { // loop over number of samples (1, 2 or 3)
				get_biopsy_data(rao, depth_diversity[i], num_samples_list[i], biopsy_size_per_sample, num_demes, num_cells, num_driver_matrix_cols, num_clones, clone_ints[DRIVER_GENOTYPE], driver_matrix, 
					dmax, idum, sample_size_log, gens_elapsed, output_driver_phylo, 1, record_phylogenies, centre_X, centre_Y, driver_genotype_ints, driver_genotype_floats, at_edge);
				get_biopsy_data(rao, depth_diversity_bigsample[i], num_samples_list[i], 10*biopsy_size_per_sample, num_demes, num_cells, num_driver_matrix_cols, num_clones, clone_ints[DRIVER_GENOTYPE], driver_matrix, 
					dmax, idum, sample_size_log, gens_elapsed, output_driver_phylo, 1, record_phylogenies, centre_X, centre_Y, driver_genotype_ints, driver_genotype_floats, at_edge);
			}
		}

		// write diversity metrics to file:
		write_diversities(output_diversities, diversity, alpha_diversity, edge_diversity, driver_diversity, alpha_driver_diversity, edge_driver_diversity, depth_diversity, 
			depth_diversity_bigsample, gens_elapsed, num_cells);

		// write to other output files:
		write_other_files(output_demes, output_clones, output_genotype_counts, output_driver_genotype_counts, output_phylo, gens_elapsed, num_demes, num_matrix_cols, num_driver_matrix_cols, num_clones, 
			within_deme_diversity, within_deme_driver_diversity, num_cells);

		if(mu_passenger > 0) {
			get_relatives(num_matrix_cols, next_genotype_id, genotype_ints);
			get_allele_frequencies(num_matrix_cols, allele_count, genotype_ints);
			get_frequency_table(num_matrix_cols, allele_count, freq_table, &num_freqs);
			write_frequency_table(output_allele_counts, freq_table, num_cells, gens_elapsed, num_freqs);

			get_frequency_table(num_matrix_cols, genotype_ints[POPULATION], freq_table, &num_freqs);
			write_frequency_table(output_genotype_counts, freq_table, num_cells, gens_elapsed, num_freqs);
		}
		
		get_relatives(num_driver_matrix_cols, next_driver_genotype_id, driver_genotype_ints);
		get_allele_frequencies(num_driver_matrix_cols, allele_count, driver_genotype_ints);
		get_frequency_table(num_driver_matrix_cols, allele_count, freq_table, &num_freqs);
		write_frequency_table(output_driver_allele_counts, freq_table, num_cells, gens_elapsed, num_freqs);

		get_frequency_table(num_driver_matrix_cols, driver_genotype_ints[POPULATION], freq_table, &num_freqs);
		write_frequency_table(output_driver_genotype_counts, freq_table, num_cells, gens_elapsed, num_freqs);
	}

	if(component(which_parts, GENOPROPS)) {
		// write properties of genotypes and driver genotypes:
		get_relatives(num_matrix_cols, next_genotype_id, genotype_ints);
		get_allele_frequencies(num_matrix_cols, allele_count, genotype_ints);
		write_genotypes(output_genotype_properties, num_matrix_cols, allele_count, genotype_ints, genotype_floats);
		get_relatives(num_driver_matrix_cols, next_driver_genotype_id, driver_genotype_ints);
		get_allele_frequencies(num_driver_matrix_cols, allele_count, driver_genotype_ints);
		write_genotypes(output_driver_genotype_properties, num_driver_matrix_cols, allele_count, driver_genotype_ints, driver_genotype_floats);
		// write relatedness matrices:
		if(record_matrix && mu_passenger > 0) write_matrix_to_file(output_matrix, matrix, num_matrix_cols);
		if(matrix_max > 0) write_matrix_to_file(output_driver_matrix, driver_matrix, num_driver_matrix_cols);
		
		if(mu_passenger > 0) {
			// record genotype phylogenies of sampled cells (but don't calculate diversity):
			if(write_phylo > 0.5 && !filled_grid) for(i = 0; i < 3; i++) { // loop over number of samples (1, 2 or 3)
				get_biopsy_data(1, depth_diversity[i], num_samples_list[i], biopsy_size_per_sample, num_demes, num_cells, num_matrix_cols, num_clones, clone_ints[GENOTYPE], matrix, 
					dmax, idum, sample_size_log, gens_elapsed, output_phylo, 0, record_phylogenies, centre_X, centre_Y, genotype_ints, genotype_floats, at_edge);
				get_biopsy_data(1, depth_diversity_bigsample[i], num_samples_list[i], 10*biopsy_size_per_sample, num_demes, num_cells, num_matrix_cols, num_clones, clone_ints[GENOTYPE], matrix, 
					dmax, idum, sample_size_log, gens_elapsed, output_phylo, 0, record_phylogenies, centre_X, centre_Y, genotype_ints, genotype_floats, at_edge);
			}
		}
	}

	if(print_output > 0) print_to_screen(gens_elapsed, num_cells, t1, num_demes, num_matrix_cols, num_clones, mean_num_passengers, mean_num_drivers, diversity, driver_diversity, sum_birth_rates,
		sum_death_rates, sum_migration_rates, num_empty_cols, alpha_diversity, alpha_driver_diversity, event_counter, num_driver_matrix_cols, sum_normal_death_rates, 
		sum_normal_birth_rates, edge_driver_diversity, num_empty_driver_cols, num_empty_demes, num_extinct_genotypes, num_extinct_driver_genotypes, next_genotype_id);
	
	if(print_output == 2) {
		if(iterations > 0) {
			printf("----------------------------------------------------------------\n");
			printf("%ld seconds elapsed\n", (long)time(NULL)-t1);
		}
		printf("################################################################\n");
	}
}

// write and/or plot grids:
void grids_output(char *preamble_text, char *preamble_drivers_text, char *preamble_passengers_text, float gens_elapsed, int num_clones, int num_demes, char *input_and_output_path, 
	char *buffer_text_short, char *buffer_text_long, bool to_file)
{
	if(write_grid) {
		plot_birth_rate_grid(gp, preamble_text, gens_elapsed, num_clones, input_and_output_path, buffer_text_short, buffer_text_long);
		plot_drivers_grid(gp, preamble_drivers_text, gens_elapsed, input_and_output_path, buffer_text_short, buffer_text_long);
		plot_pops_grid(gp, preamble_text, gens_elapsed, input_and_output_path, buffer_text_short, buffer_text_long);
		plot_passengers_grid(gp, preamble_passengers_text, gens_elapsed, num_clones, num_demes, input_and_output_path, buffer_text_short, buffer_text_long);
		plot_migration_grid(gp, preamble_text, gens_elapsed, num_clones, input_and_output_path, buffer_text_short, buffer_text_long);
	}

	if(to_file) {
		write_pops_grid_to_file(output_popgrid);
		write_passengers_grid_to_file(output_passengersgrid, num_clones, num_demes);
		write_normalcells_grid_to_file(output_normalcellsgrid);
		write_deathrates_grid_to_file(output_deathrates_grid);
		write_rates_grid_to_file(output_birthratesgrid, deme_doubles[SUM_BIRTH_RATES], deme_ints[POPULATION]);
		write_rates_grid_to_file(output_migrationratesgrid, deme_doubles[SUM_MIGRATION_RATES], deme_ints[POPULATION]);
		write_drivers_grid_to_file(output_driversgrid);
	}
}

// write and print at the end of a single simulation (not necessarily end of execution):
void end_of_loop_output(int num_cells, float gens_elapsed, long t1)
{
	printf("\n");
	if((long)time(NULL)-t1 >= max_time) {
		printf("Reached time limit\n");
		fprintf(error_log, "Reached time limit\n");
		exit_code += 2;
	}
	if(num_cells >= max_pop) {
		printf("Reached max pop (%d cells)\n", num_cells);
		fprintf(error_log, "Reached max pop (%d cells)\n", num_cells);
	}
	if(num_cells <= 0) {
		printf("No more cells\n");
		fprintf(error_log, "No more cells\n");
		exit_code += 4;
	}
	if(gens_elapsed >= max_generations) {
		printf("Reached max generations (%f >= %d)\n", gens_elapsed, max_generations);
		fprintf(error_log, "Reached max generations\n");
	}
	printf("\n");
}

// write and print at end of execution:
void final_output(int trial_num, int num_cells, long t1)
{
	fprintf(error_log, "Ran %d trials (max permitted %d trials)\n", trial_num + 1, MAX_TRIALS);
	fprintf(error_log, "Reached %d iterations and %d cells\n", iterations, num_cells);
	printf("Completed in %ld seconds.\n", (long)time(NULL)-t1);
	fprintf(error_log, "Completed in %ld seconds.\n", (long)time(NULL)-t1);
	printf("Exit code %d\n", exit_code);
	fprintf(error_log, "Exit code %d\n", exit_code);
}

// close output files:
void close_files()
{
	fclose(output_parameters);
	fclose(output_pops);
	fclose(output_diversities);
	fclose(output_demes);
	fclose(output_clones);
	fclose(output_genotype_counts);
	fclose(output_driver_genotype_counts);
	fclose(output_phylo);
	fclose(output_driver_phylo);
	fclose(output_matrix);
	fclose(output_driver_matrix);
	fclose(output_popgrid);
	fclose(output_passengersgrid);
	fclose(output_normalcellsgrid);
	fclose(output_deathrates_grid);
	fclose(output_birthratesgrid);
	fclose(output_migrationratesgrid);
	fclose(output_driversgrid);
	fclose(error_log);
	fclose(sample_size_log);
	fclose(output_allele_counts);
	fclose(output_driver_allele_counts);
	fclose(output_genotype_properties);
	fclose(output_driver_genotype_properties);
}

/////////////////////// write system state to screen or file:

// print updates to screen:
void print_to_screen(float gens_elapsed, int num_cells, long t1, int num_demes, int num_matrix_cols, int num_clones, float mean_num_passengers, 
	float mean_num_drivers, float diversity, float driver_diversity, double sum_birth_rates, double sum_death_rates, double sum_migration_rates, int num_empty_cols, 
	float alpha_diversity, float alpha_driver_diversity, int *event_counter, int num_driver_matrix_cols, double sum_normal_death_rates, double sum_normal_birth_rates, float edge_diversity, 
	int num_empty_driver_cols, int num_empty_demes, int num_extinct_genotypes, int num_extinct_driver_genotypes, int next_genotype_id)
{
	int i;
	int normal_pop = 0;

	for(i=0; i<num_demes; i++) normal_pop += deme_ints[NORMAL_CELLS][i];

	printf("Seed = %d; K = %d\n", seed, K);
	printf("Migration type = %d, initial rate = %e, edge only = %d\n", migration_type, init_migration_rate, migration_edge_only);
	printf("Mut rates: %f (birth), %f (mig), %f (pass)\n", mu_driver_birth, mu_driver_migration, mu_passenger);
	printf("Mut effects: %f (birth), %f (mig), %f (pass)\n", s_driver_birth, s_driver_migration, s_passenger);
	printf("----------------------------------------------------------------\n");
	printf("%d generations, %d iterations\n", (int)(gens_elapsed), iterations);
	printf("%d cells, %d clones, %d genotypes, %d driver genotypes\n", num_cells, num_clones, num_matrix_cols - num_empty_cols - num_extinct_genotypes, 
		num_driver_matrix_cols - num_empty_driver_cols - num_extinct_driver_genotypes);
	printf("%d matrix columns; %d genotypes ever created\n", num_matrix_cols, next_genotype_id);
	printf("%d demes (%d not empty), %d bintree layers", num_demes, num_demes - num_empty_demes, max_layer_needed_array[num_demes] + 1);
	printf("\n");
	printf("----------------------------------------------------------------\n");
	printf("Mean number of passengers, drivers = %f, %f\n", mean_num_passengers, mean_num_drivers);
	if(calculate_total_diversity) printf("Diversity = %f (alpha = %f, beta = %f)\n", diversity, alpha_diversity, diversity / alpha_diversity);
	printf("Driver diversity = %f (alpha = %f, beta = %f)\n", driver_diversity, alpha_driver_diversity, driver_diversity / alpha_driver_diversity);
	if(migration_type < 2) printf("Mean birth, death, mig rates = %f, %f, %f\n", sum_birth_rates / num_cells, sum_death_rates / num_cells, sum_migration_rates / num_cells);
	else printf("Mean birth, death, fission rates = %f, %f, %.3e\n", sum_birth_rates / num_cells, sum_death_rates / num_cells, sum_migration_rates / num_cells);
	if(normal_pop > 0) printf("Mean normal birth, death rates = %f, %f\n", sum_normal_birth_rates / normal_pop, sum_normal_death_rates / normal_pop);
	else printf("No normal cells in tumour demes\n");
	if(iterations > 0) {
		printf("----------------------------------------------------------------\n");
		printf("%d births, %d deaths, ", event_counter[BIRTH_EVENT], event_counter[DEATH_EVENT]);
		if(migration_type < 2) printf("%d migrations, ", event_counter[MIGRATION_EVENT]);
		else printf("%d fissions, ", event_counter[FISSION_EVENT]);
		printf("%d mutations\n", event_counter[MUTATION_EVENT]);
		printf("%d normal cell births, %d normal cell deaths\n", event_counter[NORMALBIRTH_EVENT], event_counter[NORMALDEATH_EVENT]);
	}

	// error check:
	check_matrix_cols2(num_matrix_cols, num_driver_matrix_cols, num_empty_cols, num_empty_driver_cols, num_extinct_genotypes, num_extinct_driver_genotypes);
}

// write overall population metrics to file:
void write_output_pops(FILE *output_pops, int num_cells, int *event_counter, float mean_num_passengers, float mean_num_drivers, double sum_birth_rates, double sum_migration_rates, 
	float var_num_passengers, float var_num_drivers, float variance_birth_rate, float variance_mig_rate, int num_matrix_cols, int num_empty_cols, int num_driver_matrix_cols, int num_empty_driver_cols, 
	int num_extinct_genotypes, int num_extinct_driver_genotypes, int num_demes, int num_empty_demes, float gens_elapsed, int num_clones, int *driver_counts, int next_genotype_id)
{
	int i;
	int mig_or_fission_count;

	if(migration_type < 2) mig_or_fission_count = total_events[MIGRATION_EVENT] + event_counter[MIGRATION_EVENT];
	else mig_or_fission_count = total_events[FISSION_EVENT] + event_counter[FISSION_EVENT];

	fprintf(output_pops, "%f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t", gens_elapsed, num_cells, total_events[BIRTH_EVENT] + event_counter[BIRTH_EVENT], total_events[DEATH_EVENT] + event_counter[DEATH_EVENT], 
		mig_or_fission_count, mean_num_passengers, mean_num_drivers, var_num_passengers, var_num_drivers);
	fprintf(output_pops, "%e\t", sum_birth_rates / num_cells);
	fprintf(output_pops, "%e\t", sum_migration_rates / num_cells);
	fprintf(output_pops, "%e\t%e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", variance_birth_rate, variance_mig_rate, num_clones, num_matrix_cols - num_empty_cols - num_extinct_genotypes, 
		num_driver_matrix_cols - num_empty_driver_cols - num_extinct_driver_genotypes, num_demes, num_demes - num_empty_demes, num_matrix_cols, num_driver_matrix_cols, next_genotype_id);
	for(i=0; i<=MAX_DRIVERS_TO_COUNT; i++) fprintf(output_pops, "\t%d", driver_counts[i]);
	fprintf(output_pops, "\n");
	fflush(output_pops);
}

// write diversity metrics to file:
void write_diversities(FILE *output_diversities, float diversity, float alpha_diversity, float edge_diversity, float driver_diversity, float alpha_driver_diversity, float edge_driver_diversity, 
	float **depth_diversity, float **depth_diversity_bigsample, float gens_elapsed, int num_cells)
{
	int i, j;
	int min_depth;

	if(well_mixed) min_depth = 11;
	else if(filled_grid) min_depth = 12;
	else min_depth = 0;

	fprintf(output_diversities, "%f\t%d\t", gens_elapsed, num_cells);
	if(calculate_total_diversity) fprintf(output_diversities, "%f\t%f\t%f\t", diversity, edge_diversity, alpha_diversity);
	fprintf(output_diversities, "%f\t%f\t%f", driver_diversity, edge_driver_diversity, alpha_driver_diversity);
	for(i=0; i<3; i++) for(j=min_depth; j<=11; j++) fprintf(output_diversities, "\t%f", depth_diversity[i][j]);
	for(i=0; i<3; i++) for(j=min_depth; j<=11; j++) fprintf(output_diversities, "\t%f", depth_diversity_bigsample[i][j]);
	fprintf(output_diversities, "\n");
}

// write deme, genotype, and phylogeny properties to file:
void write_other_files(FILE *output_demes, FILE *output_clones, FILE *output_genotype_counts, FILE *output_driver_genotype_counts, FILE* output_phylo, float gens_elapsed, 
	int num_demes, int num_matrix_cols, int num_driver_matrix_cols, int num_clones, float *within_deme_diversity, float *within_deme_driver_diversity, int num_cells)
{
	int i;
	int geno_num, deme_num, driver_geno_num;

	// write deme sizes, locations, and diversities to file:
	if(write_demes_file) for(i=0; i<num_demes; i++) fprintf(output_demes, "%f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n", gens_elapsed, deme_ints[XCOORD][i], deme_ints[YCOORD][i],
		deme_ints[POPULATION][i], deme_ints[NORMAL_CELLS][i], deme_floats[DEATH_RATE][i], within_deme_diversity[i], within_deme_driver_diversity[i]);

	// write clone properties to file:
	if(write_clones_file) for(i=0; i<num_clones; i++) {
		geno_num = clone_ints[GENOTYPE][i];
		driver_geno_num = clone_ints[DRIVER_GENOTYPE][i];
		deme_num = clone_ints[DEME][i];
		fprintf(output_clones, "%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\n", gens_elapsed, i, deme_num, genotype_ints[IDENTITY][geno_num], 
			driver_genotype_ints[IDENTITY][driver_geno_num], genotype_ints[PARENT][geno_num], driver_genotype_ints[PARENT][driver_geno_num], 
			deme_ints[XCOORD][deme_num], deme_ints[YCOORD][deme_num], deme_ints[NORMAL_CELLS][deme_num], clone_ints[POPULATION][i], genotype_floats[BIRTH_RATE][geno_num], 
			genotype_floats[MIGRATION_RATE][geno_num], deme_floats[DEATH_RATE][deme_num], deme_floats[MIGRATION_MODIFIER][deme_num], genotype_ints[NUM_DRIVER_MUTATIONS][geno_num], 
			genotype_ints[NUM_MIGRATION_MUTATIONS][geno_num], genotype_ints[NUM_PASSENGER_MUTATIONS][geno_num]);
	}

	// write phylogenetic data of genotypes to file:
	if(write_phylo > 0) write_output_phylo(output_phylo, num_matrix_cols, gens_elapsed, genotype_ints[POPULATION], genotype_ints, genotype_floats, 1, -1, -1, num_cells);
}

// write phylogenetic data of genotypes to file:
void write_output_phylo(FILE* output, int num_cols, float gens_elapsed, int *populations, int **genotype_or_driver_ints, float **genotype_or_driver_floats, 
	int samples, int biopsy_size_per_sample, int depth, int num_cells)
{
	int i;

	for(i=0; i<num_cols; i++) if(genotype_or_driver_ints[IMMORTAL][i] == 1) fprintf(output, "%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%d\n", num_cells, samples, biopsy_size_per_sample, depth, gens_elapsed, 
		genotype_or_driver_ints[IDENTITY][i], genotype_or_driver_ints[DRIVER_IDENTITY][i], genotype_or_driver_ints[PARENT][i], populations[i], genotype_or_driver_floats[BIRTH_RATE][i], genotype_or_driver_floats[MIGRATION_RATE][i], 
		genotype_or_driver_floats[ORIGIN_TIME][i], genotype_or_driver_ints[NUM_DRIVER_MUTATIONS][i], genotype_or_driver_ints[NUM_MIGRATION_MUTATIONS][i], genotype_or_driver_ints[NUM_PASSENGER_MUTATIONS][i]);
}

void write_frequency_table(FILE *filename, int **freq_table, int denominator, float gens_elapsed, int array_length)
{
	int i;

	for(i = 0; i < array_length; i++) fprintf(filename, "%f\t%d\t%f\t%d\n", gens_elapsed, freq_table[0][i], (float)freq_table[0][i] / denominator, freq_table[1][i]);
}

void write_genotypes(FILE *output_genotype_properties, int num_matrix_cols, int *allele_count, int **genotype_or_driver_ints, float **genotype_or_driver_floats)
{
	int i;

	for(i = 0; i < num_matrix_cols; i++) if(allele_count[i] > 0) {
		fprintf(output_genotype_properties, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", genotype_or_driver_ints[POPULATION][i], genotype_or_driver_ints[PARENT][i], genotype_or_driver_ints[IDENTITY][i], genotype_or_driver_ints[DRIVER_IDENTITY][i],
			genotype_or_driver_ints[NUM_DRIVER_MUTATIONS][i], genotype_or_driver_ints[NUM_MIGRATION_MUTATIONS][i], genotype_or_driver_ints[IMMORTAL][i], genotype_or_driver_ints[NUM_PASSENGER_MUTATIONS][i]);
		fprintf(output_genotype_properties, "%f\t%f\t%f\t", genotype_or_driver_floats[BIRTH_RATE][i], genotype_or_driver_floats[MIGRATION_RATE][i], genotype_or_driver_floats[ORIGIN_TIME][i]);
		fprintf(output_genotype_properties, "%d\n", allele_count[i]);
	}
}

// write distance matrix to file:
void write_matrix_to_file(FILE *output_matrix, int **either_matrix, int num_matrix_cols)
{
	int i, j;

	for(i = 0; i < num_matrix_cols; i++) {
		fprintf(output_matrix, "%d", either_matrix[i][0]);
		for(j = 1; j < num_matrix_cols; j++) fprintf(output_matrix, "\t%d", either_matrix[i][j]);
		fprintf(output_matrix, "\n");
		}
}

// write grid of deme populations to file:
void write_pops_grid_to_file(FILE *output_popgrid)
{
	int i, j;
	
	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			if(grid[i][j] == EMPTY) fprintf(output_popgrid, "NA");
			else fprintf(output_popgrid, "%d", deme_ints[POPULATION][grid[i][j]]);
			if(j < dim_grid - 1) fprintf(output_popgrid, "\t");
		}
		fprintf(output_popgrid, "\n");
		}
}

// write grid of mean passenger mutation counts to file:
void write_passengers_grid_to_file(FILE *output_passengersgrid, int num_clones, int num_demes)
{
	int i, j;
	int clone_num, clone_index, deme_num, geno_num, passenger_count, deme_pop;

	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			deme_num = grid[i][j];
			
			if(deme_num == EMPTY) fprintf(output_passengersgrid, "NA");
			else if(deme_ints[POPULATION][deme_num] == 0) fprintf(output_passengersgrid, "NA");
			else {
				deme_pop = deme_ints[POPULATION][deme_num];
				passenger_count = 0;
				for(clone_num=0; clone_num<deme_ints[NUM_CLONES_IN_DEME][deme_num]; clone_num++) {
					clone_index = clones_list_in_deme[deme_num][clone_num];
					geno_num = clone_ints[GENOTYPE][clone_index];
					passenger_count += clone_ints[POPULATION][clone_index] * genotype_ints[NUM_PASSENGER_MUTATIONS][geno_num];
				}
				fprintf(output_passengersgrid, "%f", (float)passenger_count / deme_pop);
			}
			if(j < dim_grid - 1) fprintf(output_passengersgrid, "\t");
		}
		fprintf(output_passengersgrid, "\n");
	}
}

// write grid of normal cell populations to file:
void write_normalcells_grid_to_file(FILE *output_normalcellsgrid)
{
	int i, j;
	
	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			if(grid[i][j] == EMPTY) {
				if(normal_birth_rate > 0) fprintf(output_normalcellsgrid, "%d", K);
				else fprintf(output_normalcellsgrid, "0");
			}
			else fprintf(output_normalcellsgrid, "%d", deme_ints[NORMAL_CELLS][grid[i][j]]);
			if(j < dim_grid - 1) fprintf(output_normalcellsgrid, "\t");
		}
		fprintf(output_normalcellsgrid, "\n");
	}
}

// write grid of death rates to file:
void write_deathrates_grid_to_file(FILE *output_deathrates_grid)
{
	int i, j;
	int deme_num;
	
	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			deme_num = grid[i][j];
			if(deme_num == EMPTY) fprintf(output_deathrates_grid, "NA");
			else fprintf(output_deathrates_grid, "%f", deme_floats[DEATH_RATE][deme_num]);
			if(j < dim_grid - 1) fprintf(output_deathrates_grid, "\t");
		}
		fprintf(output_deathrates_grid, "\n");
		}
}

// write grid of birth or migration rates to file:
void write_rates_grid_to_file(FILE *output, double *deme_rates, int *deme_pops)
{
	int i, j;
	int deme_num;

	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			deme_num = grid[i][j];
			if(deme_num == EMPTY) fprintf(output, "NA");
			else if(deme_ints[POPULATION][deme_num] == 0) fprintf(output, "NA");
			else fprintf(output, "%f", deme_rates[grid[i][j]] / deme_pops[grid[i][j]]);
			if(j < dim_grid - 1) fprintf(output, "\t");
		}
		fprintf(output, "\n");
		}
}

// write grid of birth or migration rates to file:
void write_drivers_grid_to_file(FILE *output)
{
	int i, j;
	int dominant;
	int max_clone_pop;
	int deme_num, clone_num, clone_index;

	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			deme_num = grid[i][j];
			max_clone_pop = 0;
			if(deme_num != EMPTY) for(clone_num=0; clone_num<deme_ints[NUM_CLONES_IN_DEME][deme_num]; clone_num++) {
				clone_index = clones_list_in_deme[deme_num][clone_num];
				if(clone_ints[POPULATION][clone_index] > max_clone_pop) {
					dominant = clone_ints[DRIVER_GENOTYPE][clone_index];
					max_clone_pop = clone_ints[POPULATION][clone_index];
				}
			}
			if(max_clone_pop == 0 || deme_num == EMPTY) fprintf(output, "NA"); // never occupied
			else fprintf(output, "%d", driver_genotype_ints[IDENTITY][dominant]);
			if(j < dim_grid - 1) fprintf(output, "\t");
		}
		fprintf(output, "\n");
	}
}

// create grid image, coloured by birth rate:
void plot_birth_rate_grid(FILE *gp, char *preamble_text, float gens_elapsed, int num_clones,
	char *input_and_output_path, char *buffer_text_short, char *buffer_text_long)
{
	int i,j;
	float val_to_write;
	float min_rate = 0.9, max_rate = 2;
	int deme_index;
	float buff;

	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set output '");
	fprintf(gp, "%sBirthRateGrids/birth_rate_grid%s%d.png'\n", input_and_output_path, zeroes((int)(gens_elapsed), max_gens, buffer_text_short, buffer_text_long), 
		(int)(gens_elapsed));
	fprintf(gp, preamble_text, 0);

	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			deme_index = grid[i][j];
			if(deme_index == EMPTY) val_to_write = 4.5; // never occupied
			else if(deme_ints[POPULATION][deme_index] == 0) val_to_write = 4.5; // empty
			else {
				buff = 3 * (deme_doubles[SUM_BIRTH_RATES][deme_index] / deme_ints[POPULATION][deme_index] - min_rate) / (max_rate - min_rate);
				buff = MAX(0, buff);
				val_to_write = MIN(3, buff);
			}
			sprintf(buffer_text_short, "%f ", val_to_write); // write its state
			fprintf(gp, buffer_text_short, 0);
		}
		fprintf(gp, "\n");
	}
	fprintf(gp, "e\n");
	fprintf(gp, "e\n");
}

// create grid image, coloured by driver genotype:
void plot_drivers_grid(FILE *gp, char *preamble_drivers_text, float gens_elapsed, char *input_and_output_path, char *buffer_text_short, char *buffer_text_long)
{
	int i,j,clone_index;
	float val_to_write;
	int dominant;
	int max_clone_pop;
	int deme_num,clone_num;

	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set output '");
	fprintf(gp, "%sDriverGrids/driver_grid%s%d.png'\n", input_and_output_path, zeroes((int)(gens_elapsed), max_gens, buffer_text_short, buffer_text_long), 
		(int)(gens_elapsed));
	fprintf(gp, preamble_drivers_text, 0);

	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			deme_num = grid[i][j];
			max_clone_pop = 0;
			if(deme_num != EMPTY) for(clone_num=0; clone_num<deme_ints[NUM_CLONES_IN_DEME][deme_num]; clone_num++) {
				clone_index = clones_list_in_deme[deme_num][clone_num];
				if(clone_ints[POPULATION][clone_index] > max_clone_pop) {
					dominant = clone_ints[DRIVER_GENOTYPE][clone_index];
					max_clone_pop = clone_ints[POPULATION][clone_index];
				}
			}
			if(max_clone_pop == 0 || deme_num == EMPTY) val_to_write = 26.5; // never occupied
			else val_to_write = 0.5 + driver_genotype_ints[IDENTITY][dominant] % 26;
			sprintf(buffer_text_short, "%f ", val_to_write); // write its state
			fprintf(gp, buffer_text_short, 0);
		}
		fprintf(gp, "\n");
	}
	fprintf(gp, "e\n");
	fprintf(gp, "e\n");
}

// create grid image, coloured by cancer cell population density:
void plot_pops_grid(FILE *gp, char *preamble_text, float gens_elapsed, char *input_and_output_path, char *buffer_text_short, char * buffer_text_long)
{
	int i,j;
	float val_to_write, buff;

	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set output '");
	fprintf(gp, "%sPopsGrids/pops_grid%s%d.png'\n", input_and_output_path, zeroes((int)(gens_elapsed), max_gens, buffer_text_short, buffer_text_long), 
		(int)(gens_elapsed));
	fprintf(gp, preamble_text, 0);

	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			if(grid[i][j] == EMPTY) val_to_write = 4.5; // never occupied
			else if(deme_ints[POPULATION][grid[i][j]] == 0) val_to_write = 4.5;
			else {
				buff = 1.5 * (float)deme_ints[POPULATION][grid[i][j]] / K;
				val_to_write = MIN(3, buff);
			}
			sprintf(buffer_text_short, "%f ", val_to_write); // write its state
			fprintf(gp, buffer_text_short, 0);
		}
		fprintf(gp, "\n");
	}
	fprintf(gp, "e\n");
	fprintf(gp, "e\n");
}

// create grid image, coloured by number of passenger mutations:
void plot_passengers_grid(FILE *gp, char *preamble_text, float gens_elapsed, int num_clones, int num_demes, char *input_and_output_path, 
	char *buffer_text_short, char *buffer_text_long)
{
	int i,j;
	float val_to_write, buff;
	int clone_num, clone_index, deme_num, geno_num, passenger_count, deme_pop;

	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set output '");
	fprintf(gp, "%sPassengersGrids/passengers_grid%s%d.png'\n", input_and_output_path, zeroes((int)(gens_elapsed), max_gens, buffer_text_short, buffer_text_long), 
		(int)(gens_elapsed));
	fprintf(gp, preamble_text, 0);

	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			deme_num = grid[i][j];
			if(deme_num == EMPTY) val_to_write = 4.5; // never occupied
			else if(deme_ints[POPULATION][grid[i][j]] == 0) val_to_write = 4.5; // empty
			else {
				deme_pop = deme_ints[POPULATION][deme_num];
				passenger_count = 0;
				for(clone_num=0; clone_num<deme_ints[NUM_CLONES_IN_DEME][deme_num]; clone_num++) {
					clone_index = clones_list_in_deme[deme_num][clone_num];
					geno_num = clone_ints[GENOTYPE][clone_index];
					passenger_count += clone_ints[POPULATION][clone_index] * genotype_ints[NUM_PASSENGER_MUTATIONS][geno_num];
				}
				buff = 0.5 * (float)passenger_count / deme_pop;
				val_to_write = MIN(3.99, buff);
			}
			sprintf(buffer_text_short, "%f ", val_to_write); // write its state
			fprintf(gp, buffer_text_short, 0);
		}
		fprintf(gp, "\n");
	}
	fprintf(gp, "e\n");
	fprintf(gp, "e\n");
}

// create grid image, coloured by migration rate:
void plot_migration_grid(FILE *gp, char *preamble_text, float gens_elapsed, int num_clones, char *input_and_output_path, char *buffer_text_short, char *buffer_text_long)
{
	int i,j;
	float val_to_write, buff;
	int deme_index;

	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set output '");
	fprintf(gp, "%sMigrationRateGrids/migration_rate_grid%s%d.png'\n", input_and_output_path, zeroes((int)(gens_elapsed), max_gens, buffer_text_short, buffer_text_long), 
		(int)(gens_elapsed));
	fprintf(gp, preamble_text, 0);

	for(i=0; i<dim_grid; i++) {
		for(j=0; j<dim_grid; j++) {
			deme_index = grid[i][j];
			if(deme_index == EMPTY) val_to_write = 4.5; // never occupied
			else if(deme_ints[POPULATION][deme_index] == 0) val_to_write = 4.5; // empty
			else if(init_migration_rate == 0) val_to_write = 0; // no migration
			else {
				buff = 2 * ((float)deme_doubles[SUM_MIGRATION_RATES][deme_index] / (deme_ints[POPULATION][deme_index] * init_migration_rate) - 1);
				val_to_write = MIN(3, buff); // write mean cell migration rate
			}
			sprintf(buffer_text_short, "%f ", val_to_write); // write the deme's state
			fprintf(gp, buffer_text_short, 0);
		}
		fprintf(gp, "\n");
	}
	fprintf(gp, "e\n");
	fprintf(gp, "e\n");
}

/////////////////////// set rates:

// set genotype birth rate:
float set_birth_rate(int new_birth_mutations, int new_passengers, float parent_birth_rate, long *idum)
{
	float birth_rate = parent_birth_rate;
	int i;

	for(i = 0; i < new_passengers; i++) birth_rate = birth_rate / (1 + s_passenger);
	if(max_relative_birth_rate >= 0) for(i = 0; i < new_birth_mutations; i++) birth_rate = birth_rate * (1 + s_driver_birth * (1 - birth_rate / max_relative_birth_rate) * expdev(idum));
	else for(i = 0; i < new_birth_mutations; i++) birth_rate = birth_rate * (1 + s_driver_birth * expdev(idum));
	
	if(birth_rate >= baseline_death_rate + density_dept_death_rate) {
		printf("Birth rate exceeds death rate at iterations %d\n", iterations);
		fprintf(error_log, "Birth rate exceeds death rate at iterations %d\n", iterations);
		exit_code = 1;
	}

	return birth_rate;
}

// set genotype migration rate:
float set_migration_rate(int new_mig_mutations, float parent_mig_rate, long *idum)
{
	float mig_rate = parent_mig_rate;
	int i;

	if(max_relative_migration_rate >= 0) for(i = 0; i < new_mig_mutations; i++) mig_rate = mig_rate * (1 + s_driver_migration * (1 - mig_rate / (max_relative_migration_rate * init_migration_rate)) * expdev(idum));
	else for(i = 0; i < new_mig_mutations; i++) mig_rate = mig_rate * (1 + s_driver_migration * expdev(idum));

	return mig_rate;
}

// set deme death rate:
float set_death_rate(int deme_index, int num_cells)
{
	int total_pop_in_deme = deme_ints[NORMAL_CELLS][deme_index] + deme_ints[POPULATION][deme_index];

	if(total_pop_in_deme <= K) return baseline_death_rate;
	else return baseline_death_rate + density_dept_death_rate; // arbitrary value >> birth rate
}

// set deme migration modifier (probability that deme can be invaded):
float set_migration_modifier(int normal_cells, int population)
{
	if(migration_edge_only) return((float)population / K);
	else return 0.0;
}

/////////////////////// calculate or reset sums of rates:

// reset deme and deme bintree sums by summing lower level rates:
void reset_deme_and_bintree_sums(int num_demes, int num_clones, int print_always)
{
	int i, j;
	double old_sum_demes, sum_demes, sum_bintree;
	int count_cells, count_normal;
	int min_bintree_index, max_bintree_index;
	int deme_index, clone_index, genotype_index;
	int max_layer_needed = max_layer_needed_array[num_demes];

	// reset sums of rates for demes, by summing rates of constituent clones:
	old_sum_demes = 0;
	sum_demes = 0;
	count_cells = 0;
	count_normal = 0;
	for(deme_index = 0; deme_index < num_demes; deme_index++) { // loop over all demes
		deme_doubles[SUM_BIRTH_RATES][deme_index] = 0;
		deme_doubles[SUM_MIGRATION_RATES][deme_index] = 0;
		for(j = 0; j < deme_ints[NUM_CLONES_IN_DEME][deme_index]; j++) { // loop over clones in the deme
			clone_index = clones_list_in_deme[deme_index][j];
			genotype_index = clone_ints[GENOTYPE][clone_index];
			deme_doubles[SUM_BIRTH_RATES][deme_index] += (double)genotype_floats[BIRTH_RATE][genotype_index] * (double)clone_ints[POPULATION][clone_index];
			deme_doubles[SUM_MIGRATION_RATES][deme_index] += (double)genotype_floats[MIGRATION_RATE][genotype_index] * (double)clone_ints[POPULATION][clone_index];
		}
		old_sum_demes += deme_doubles[SUM_RATES][deme_index];
		deme_doubles[SUM_RATES][deme_index] = deme_doubles[SUM_BIRTH_RATES][deme_index] + (double)deme_ints[POPULATION][deme_index] * (double)deme_floats[DEATH_RATE][deme_index];
		if(migration_type == 1 || migration_type == 3) deme_doubles[SUM_RATES][deme_index] += deme_doubles[SUM_MIGRATION_RATES][deme_index];
		deme_doubles[SUM_RATES][deme_index] += (double)deme_ints[NORMAL_CELLS][deme_index] * ((double)normal_birth_rate + (double)deme_floats[DEATH_RATE][deme_index]);
		sum_demes += deme_doubles[SUM_RATES][deme_index];
		count_cells += deme_ints[POPULATION][deme_index];
		count_normal += deme_ints[NORMAL_CELLS][deme_index];
	}

	// find the old sum of rates across the highest bintree layer:
	sum_bintree = 0;
	min_bintree_index = min_deme_bintree_index_by_layer[max_layer_needed]; // smallest bintree index in the layer
	max_bintree_index = get_deme_bintree_index(max_layer_needed, num_demes - 1); // largest bintree index in the layer
	for(i = min_bintree_index; i <= max_bintree_index; i++) sum_bintree += bintree_deme_doubles[i];
	
	// warn if the rounding error exceeds a threshold:
	if(fabs(sum_bintree - sum_demes) > 0.05 || fabs(old_sum_demes - sum_demes) > 0.05 || print_always == 1) {
		printf("Rounding error %f or %f > 0.05\n", fabs(sum_bintree - sum_demes), fabs(old_sum_demes - sum_demes));
		printf("sums %.2f (bintree), %.2f (demes old), %.2f (demes new) at iterations %d\n", sum_bintree, old_sum_demes, sum_demes, iterations);
		printf("cells %d, normal cells %d\n", count_cells, count_normal);
		fprintf(error_log, "Rounding error %f or %f > 0.05 at iterations %d\n", fabs(sum_bintree - sum_demes), fabs(old_sum_demes - sum_demes), iterations);
	}

	// reset bintree:
	set_bintree_sums_layer0(bintree_deme_doubles, deme_doubles[SUM_RATES], num_demes);
	set_bintree_sums_subsequent_layers(bintree_deme_doubles, num_demes, max_layer_needed, max_demes);
}

// set values for one clone in an array of clones in a deme:
void set_clone_in_deme(int deme_num, int geno_num, int index_in_deme, int clone_pop)
{
	clones_pops_in_deme[deme_num][index_in_deme] = clone_pop;
	clones_rates_in_deme[deme_num][index_in_deme] = (double)clone_pop * genotype_floats[BIRTH_RATE][geno_num];
	if(migration_type == 1 || migration_type == 3) clones_rates_in_deme[deme_num][index_in_deme] += (double)clone_pop * genotype_floats[MIGRATION_RATE][geno_num];
}

// reset clone and clone bintree sums by summing lower level rates:
void reset_clone_bintree_sums(int num_demes, int num_clones)
{
	int i, max_layer_needed, num_clones_in_deme, chosen_clone, deme_num, geno_num, index_in_deme, clone_pop, deme_pop;
	float old_rate;
	int old_pop;
	float threshold;

	for(deme_num = 0; deme_num < num_demes; deme_num++) {
		num_clones_in_deme = deme_ints[NUM_CLONES_IN_DEME][deme_num];
		for(index_in_deme = 0; index_in_deme < num_clones_in_deme; index_in_deme++) {
			chosen_clone = clones_list_in_deme[deme_num][index_in_deme];
			geno_num = clone_ints[GENOTYPE][chosen_clone];
			clone_pop = clone_ints[POPULATION][chosen_clone];
			old_rate = clones_rates_in_deme[deme_num][index_in_deme];
			old_pop = clones_pops_in_deme[deme_num][index_in_deme];
			set_clone_in_deme(deme_num, geno_num, index_in_deme, clone_pop);
			deme_pop = deme_ints[POPULATION][deme_num];
			threshold = MAX(deme_pop * 5e-7, 0.05);
			// check that the new values match the old ones:
			if(fabs(clones_rates_in_deme[deme_num][index_in_deme] - old_rate) > threshold) {
				printf("rate %d,%d was %f; now %f (POP %d, BIRTH_RATE %f\n", deme_num, index_in_deme, old_rate, clones_rates_in_deme[deme_num][index_in_deme], clone_pop, genotype_floats[BIRTH_RATE][geno_num]);
				fprintf(error_log, "rate %d,%d was %f; now %f\n", deme_num, index_in_deme, old_rate, clones_rates_in_deme[deme_num][index_in_deme]);
				exit_code = 1;
			}
			if(clones_pops_in_deme[deme_num][index_in_deme] != old_pop) {
				printf("pop %d,%d was %d; now %d\n", deme_num, index_in_deme, old_pop, clones_pops_in_deme[deme_num][index_in_deme]);
				fprintf(error_log, "pop %d,%d was %d; now %d\n", deme_num, index_in_deme, old_pop, clones_pops_in_deme[deme_num][index_in_deme]);
				exit_code = 1;
			}
		}
	}

	for(i = 0; i < num_demes; i++) {
		num_clones_in_deme = deme_ints[NUM_CLONES_IN_DEME][i];

		// set clone bintree sums for layer 0:
		set_bintree_sums_layer0(bintree_clone_doubles[i], clones_rates_in_deme[i], num_clones_in_deme);
		set_bintree_sums_layer0_int(bintree_clone_ints[i], clones_pops_in_deme[i], num_clones_in_deme);

		// set clone bintree sums for subsequent layers:
		max_layer_needed = max_layer_needed_array[num_clones_in_deme];
		set_bintree_sums_subsequent_layers(bintree_clone_doubles[i], num_clones_in_deme, max_layer_needed, max_clones_per_deme);
		set_bintree_sums_subsequent_layers_int(bintree_clone_ints[i], num_clones_in_deme, max_layer_needed, max_clones_per_deme);
	}
}

// calculate summs of death, birth and migration rates:
void calculate_sums_of_rates(double *sum_death_rates, double *sum_birth_rates, double *sum_migration_rates, double *sum_normal_birth_rates, double *sum_normal_death_rates,
	int num_demes, int num_matrix_cols, float gens_elapsed)
{
	int i;

	*sum_death_rates = 0;
	*sum_migration_rates = 0;
	*sum_birth_rates = 0;
	*sum_normal_birth_rates = 0;
	*sum_normal_death_rates = 0;
	for(i=0; i<num_demes; i++) {
		*sum_death_rates += deme_ints[POPULATION][i] * deme_floats[DEATH_RATE][i];
		*sum_migration_rates += deme_doubles[SUM_MIGRATION_RATES][i];
		*sum_birth_rates += deme_doubles[SUM_BIRTH_RATES][i];
		*sum_normal_birth_rates += deme_ints[NORMAL_CELLS][i] * normal_birth_rate;
		*sum_normal_death_rates += deme_ints[NORMAL_CELLS][i] * deme_floats[DEATH_RATE][i];
	}
}

float sum_of_all_rates(int num_demes)
{
	float temp_sum = 0;
	int i;
	int max_layer_needed = max_layer_needed_array[num_demes];
	int min_bintree_index = min_deme_bintree_index_by_layer[max_layer_needed]; // smallest bintree index in the highest layer
	int max_bintree_index = get_deme_bintree_index(max_layer_needed, num_demes - 1); // largest bintree index in the highest layer
	
	for(i = min_bintree_index; i <= max_bintree_index; i++) temp_sum += (float)bintree_deme_doubles[i];

	return temp_sum;
}

// update deme and bintree rates for changes in population size:
void update_deme_bintree_environmental_rates(int change_pop, int change_normal_cells, int num_demes, int deme_index, int num_cells)
{
	int max_layer_needed = max_layer_needed_array[num_demes];

	// subtract the old SUM_RATES from bintree elements at each layer:
	update_bintree_layers(-deme_doubles[SUM_RATES][deme_index], bintree_deme_doubles, deme_index, max_layer_needed, max_demes);

	// subtract the old death rates from the deme SUM_RATES:
	deme_doubles[SUM_RATES][deme_index] -= ((double)deme_ints[POPULATION][deme_index] - (double)change_pop) * (double)deme_floats[DEATH_RATE][deme_index];
	deme_doubles[SUM_RATES][deme_index] -= ((double)deme_ints[NORMAL_CELLS][deme_index] - (double)change_normal_cells) * (double)deme_floats[DEATH_RATE][deme_index];

	// update death rate and migration modifier, based on the new cell counts:
	deme_floats[DEATH_RATE][deme_index] = set_death_rate(deme_index, num_cells);
	deme_floats[MIGRATION_MODIFIER][deme_index] = set_migration_modifier(deme_ints[NORMAL_CELLS][deme_index], deme_ints[POPULATION][deme_index]);

	// add the new death rates to the deme SUM_RATES:
	deme_doubles[SUM_RATES][deme_index] += (double)deme_ints[POPULATION][deme_index] * (double)deme_floats[DEATH_RATE][deme_index];
	deme_doubles[SUM_RATES][deme_index] += (double)deme_ints[NORMAL_CELLS][deme_index] * (double)deme_floats[DEATH_RATE][deme_index];

	// adjust SUM_RATES for any change in normal cells:
	deme_doubles[SUM_RATES][deme_index] += (double)change_normal_cells * (double)normal_birth_rate;

	// add the new SUM_RATES to bintree elements at each layer:
	update_bintree_layers(deme_doubles[SUM_RATES][deme_index], bintree_deme_doubles, deme_index, max_layer_needed, max_demes);
}

// update deme and bintree rates for changes in birth or migration rates:
void update_deme_bintree_phenotype_rates(int chosen_clone, int deme_index, int change, int num_demes)
{
	int genotype_index = clone_ints[GENOTYPE][chosen_clone];
	int max_layer_needed = max_layer_needed_array[num_demes]; // find highest bintree layer needed
	double sum_addend;

	deme_doubles[SUM_BIRTH_RATES][deme_index] += (double)change * (double)genotype_floats[BIRTH_RATE][genotype_index];
	deme_doubles[SUM_MIGRATION_RATES][deme_index] += (double)change * (double)genotype_floats[MIGRATION_RATE][genotype_index];

	sum_addend = (double)change * (double)genotype_floats[BIRTH_RATE][genotype_index];
	if(migration_type == 1 || migration_type == 3) sum_addend += (double)change * (double)genotype_floats[MIGRATION_RATE][genotype_index];
	deme_doubles[SUM_RATES][deme_index] += sum_addend;

	// add the new SUM_RATES to bintree elements at each layer:
	update_bintree_layers(sum_addend, bintree_deme_doubles, deme_index, max_layer_needed, max_demes);
}

/////////////////////// checks:

// check clone bintrees for each deme:
void check_clone_bintree_sums(bool full)
{
	int i, j, k, layer, max_layer_needed, start_index, end_index, chosen_deme, chosen_clone, parent_geno_num, num_clones_in_deme;
	double sum_array, sum_bintree, summand;

	for(i = 0; i < dim_grid; i++) {
		for(j = 0; j < dim_grid; j++) {
			chosen_deme = grid[i][j];
			if(chosen_deme != EMPTY) {

				num_clones_in_deme = deme_ints[NUM_CLONES_IN_DEME][chosen_deme];

				if(full) for(k = 0; k < num_clones_in_deme; k++) printf("%d ", clones_pops_in_deme[chosen_deme][k]);
				if(full) if(num_clones_in_deme > 0) printf("= ");
				// sum populations of all clones in the deme:
				sum_array = 0;
				for(k = 0; k < num_clones_in_deme; k++) {
					chosen_clone = clones_list_in_deme[chosen_deme][k];
					if(full) printf("%d ", clone_ints[POPULATION][chosen_clone]);
					sum_array += clone_ints[POPULATION][chosen_clone];
				}
				if(num_clones_in_deme > 0) {
					max_layer_needed = max_layer_needed_array[num_clones_in_deme];
					if(full) printf("; ");
					for(layer = 0; layer <= max_layer_needed; layer++) {
						if(full) printf("deme %d, layer %d: ", chosen_deme, layer);
						start_index = min_clone_bintree_index_by_layer[layer];
						end_index = get_clone_bintree_index(layer, num_clones_in_deme - 1);
						// sum across bintrees:
						sum_bintree = 0;
						for(k = start_index; k <= end_index; k++) {
							if(full) printf("%d ", bintree_clone_ints[chosen_deme][k]);
							sum_bintree += bintree_clone_ints[chosen_deme][k];
						}
						// check that the sums are equal:
						if(sum_array != sum_bintree) {
							printf("\nPop sums deme %d: %f != %f (iterations %d)\n", chosen_deme, sum_array, sum_bintree, iterations);
							fprintf(error_log, "Pop sums deme %d: %f != %f (iterations %d)\n", chosen_deme, sum_array, sum_bintree, iterations);
							exit_code = 1;
						}
					}
					if(full) printf("\n");
				}

				// sum of birth and migration rates:
				if(full) for(k = 0; k < num_clones_in_deme; k++) printf("%f ", clones_rates_in_deme[chosen_deme][k]);
				if(full) if(num_clones_in_deme > 0) printf("= ");
				// sum across clones in the deme:
				sum_array = 0;
				for(k = 0; k < num_clones_in_deme; k++) {
					chosen_clone = clones_list_in_deme[chosen_deme][k];
					parent_geno_num = clone_ints[GENOTYPE][chosen_clone];
					summand = clone_ints[POPULATION][chosen_clone] * genotype_floats[BIRTH_RATE][parent_geno_num];
					if(migration_type == 1 || migration_type == 3) summand += clone_ints[POPULATION][chosen_clone] * genotype_floats[MIGRATION_RATE][parent_geno_num];
					if(full) printf("%f ", summand);
					sum_array += summand;
				}
				// sum at each bintree layer:
				if(num_clones_in_deme > 0) {
					max_layer_needed = max_layer_needed_array[num_clones_in_deme];
					if(full) printf("\n");

					for(layer = 0; layer <= max_layer_needed; layer++) {
						if(full) printf("\ndeme %d, layer %d: ", chosen_deme, layer);
						start_index = min_clone_bintree_index_by_layer[layer];
						end_index = get_clone_bintree_index(layer, num_clones_in_deme - 1);
						sum_bintree = 0;
						for(k = start_index; k <= end_index; k++) {
							if(full) printf("%f ", bintree_clone_doubles[chosen_deme][k]);
							sum_bintree += bintree_clone_doubles[chosen_deme][k];
						}
						if(fabs((sum_array - sum_bintree) / sum_bintree) > 1e-3 || fabs(sum_array - sum_bintree) > 0.1) {
							printf("\nRate sums: %f != %f (iterations %d; chosen_deme %d, layer %d)\n", sum_array, sum_bintree, iterations, chosen_deme, layer);
							fprintf(error_log, "\nRate sums: %f != %f (iterations %d)\n", sum_array, sum_bintree, iterations);
							exit_code = 1;
						}
					}
					if(full) printf("\n");
				}
			}
		}
	}
}

void check_deme_sums(int chosen_deme, float deme_sum_cancer)
{
	if(fabs(buff_array[deme_ints[NUM_CLONES_IN_DEME][chosen_deme] - 1] -  deme_sum_cancer) / deme_sum_cancer > 1E-5) {
		printf("Deme sums don't match (%f and %f)\n", buff_array[deme_ints[NUM_CLONES_IN_DEME][chosen_deme] - 1], deme_sum_cancer);
		fprintf(error_log, "Deme sums don't match (%f and %f)\n", buff_array[deme_ints[NUM_CLONES_IN_DEME][chosen_deme] - 1], deme_sum_cancer);
		exit_code = 1;
	}
}

void check_rates_sum(float b0, float b1, int chosen_deme, int cell_type)
{
	if(b0 + b1 <= 0) {
		printf("Chosen deme rates <= 0; iterations %d; cell type = %d; buff_array = %f (%f+%d*%f), %f; SUM_RATES = %f\n", iterations, cell_type, b0, 
			deme_doubles[SUM_BIRTH_RATES][chosen_deme], deme_ints[POPULATION][chosen_deme], deme_floats[DEATH_RATE][chosen_deme], buff_array[1], deme_doubles[SUM_RATES][chosen_deme]);
		fprintf(error_log, "Chosen deme rates <= 0; iterations %d; cell type = %d; buff_array = %f (%f+%d*%f), %f; SUM_RATES = %f\n", iterations, cell_type, buff_array[0], 
			deme_doubles[SUM_BIRTH_RATES][chosen_deme], deme_ints[POPULATION][chosen_deme], deme_floats[DEATH_RATE][chosen_deme], buff_array[1], deme_doubles[SUM_RATES][chosen_deme]);
		exit_code = 1;
	}
}

void check_geno_populations(int chosen_clone, int event_type)
{
	if(mu_passenger == 0 && passenger_pop_threshold >= max_pop && genotype_ints[POPULATION][clone_ints[GENOTYPE][chosen_clone]] != driver_genotype_ints[POPULATION][clone_ints[DRIVER_GENOTYPE][chosen_clone]]) {
		printf("Unequal populations (%d and %d) event_type = %d", genotype_ints[POPULATION][clone_ints[GENOTYPE][chosen_clone]], 
			driver_genotype_ints[POPULATION][clone_ints[DRIVER_GENOTYPE][chosen_clone]], event_type);
		fprintf(error_log, "Unequal populations (%d and %d) event_type = %d", genotype_ints[POPULATION][clone_ints[GENOTYPE][chosen_clone]], 
			driver_genotype_ints[POPULATION][clone_ints[DRIVER_GENOTYPE][chosen_clone]], event_type);
		exit_code = 1;
	}
}

void check_matrix_cols1(int num_matrix_cols, int num_driver_matrix_cols, int num_cells)
{
	if(record_matrix && num_matrix_cols >= matrix_max - 1) {
		printf("Number of genotypes >= max matrix size; Iterations %d, num_cells %d\n", iterations, num_cells);
		fprintf(error_log, "Number of genotypes >= max matrix size; Iterations %d, num_cells %d\n", iterations, num_cells);
		exit_code = 1;
	}
	if(matrix_max > 0 && num_driver_matrix_cols >= matrix_max - 1) {
		printf("Number of driver genotypes >= max matrix size; Iterations %d, num_cells %d\n", iterations, num_cells);
		fprintf(error_log, "Number of genotypes >= max matrix size; Iterations %d, num_cells %d\n", iterations, num_cells);
		exit_code = 1;
	}
	if(num_matrix_cols >= max_genotypes - 1) {
		printf("num_matrix_cols is %d; allocated array size is %d; Iterations %d, num_cells %d\n", num_matrix_cols, max_genotypes, iterations, num_cells);
		fprintf(error_log, "num_matrix_cols is %d; allocated array size is %d; Iterations %d, num_cells %d\n", num_matrix_cols, max_genotypes, iterations, num_cells);
		exit_code = 1;
	}
	if(num_driver_matrix_cols >= max_driver_genotypes - 1) {
		printf("num_driver_matrix_cols is %d; allocated array size is %d; Iterations %d, num_cells %d\n", num_driver_matrix_cols, max_driver_genotypes, iterations, num_cells);
		fprintf(error_log, "num_driver_matrix_cols is %d; allocated array size is %d; Iterations %d, num_cells %d\n", num_driver_matrix_cols, max_driver_genotypes, iterations, num_cells);
		exit_code = 1;
	}
}

void check_matrix_cols2(int num_matrix_cols, int num_driver_matrix_cols, int num_empty_cols, int num_empty_driver_cols, int num_extinct_genotypes, int num_extinct_driver_genotypes)
{
	if(num_matrix_cols - num_empty_cols - num_extinct_genotypes < 1) {
		printf("Number of genotypes < 1 (%d - %d - %d)\n", num_matrix_cols, num_empty_cols, num_extinct_genotypes);
		fprintf(error_log, "Number of genotypes < 1 (%d - %d - %d)\n", num_matrix_cols, num_empty_cols, num_extinct_genotypes);
		exit_code = 1;
	}
	if(num_driver_matrix_cols - num_empty_driver_cols - num_extinct_driver_genotypes < 1) {
		printf("Number of driver genotypes < 1 (%d - %d - %d)\n", num_driver_matrix_cols, num_empty_driver_cols, num_extinct_driver_genotypes);
		fprintf(error_log, "Number of driver genotypes < 1 (%d - %d - %d)\n", num_driver_matrix_cols, num_empty_driver_cols, num_extinct_driver_genotypes);
		exit_code = 1;
	}
}

void check_clones_in_deme1(int origin_deme_num, int clone_num)
{
	if(deme_ints[NUM_CLONES_IN_DEME][origin_deme_num] > deme_ints[POPULATION][origin_deme_num]) {
		printf("NUM_CLONES_IN_DEME = %d > POPULATION = %d\n", deme_ints[NUM_CLONES_IN_DEME][origin_deme_num], deme_ints[POPULATION][origin_deme_num]);
		printf("iterations %d, SUM_RATES = %f, SUM_BIRTH_RATES = %f\n", iterations, deme_doubles[SUM_RATES][origin_deme_num], deme_doubles[SUM_BIRTH_RATES][origin_deme_num]);
		printf("origin_deme_num %d, clone_num %d\n", origin_deme_num, clone_num);
		fprintf(error_log, "NUM_CLONES_IN_DEME = %d > POPULATION = %d\n", deme_ints[NUM_CLONES_IN_DEME][origin_deme_num], deme_ints[POPULATION][origin_deme_num]);
		exit_code = 1;
	}
}

void check_clones_in_deme2(int num_clones_in_deme)
{
	if(num_clones_in_deme > max_clones_per_deme) {
		printf("Too many clones in deme (%d > %d)\n", num_clones_in_deme, max_clones_per_deme);
		fprintf(error_log, "Too many clones in deme\n");
		exit_code = 1;
	}
}

void check_chosen_clone(int chosen_clone, int num_clones_in_deme)
{
	if(chosen_clone > num_clones_in_deme) {
		printf("chosen_clone > num_clones_in_deme at iterations %d\n", iterations);
		fprintf(error_log, "chosen_clone > num_clones_in_deme at iterations %d\n", iterations);
		exit_code = 1;
	}
}

void check_clone_populations(int chosen_clone, int event_type, int parent_deme_num)
{
	if(clone_ints[POPULATION][chosen_clone] > deme_ints[POPULATION][parent_deme_num]) {
		printf("iterations %d, event_type %d, chosen_clone %d, parent_deme_num %d\n", iterations, event_type, chosen_clone, parent_deme_num);
		printf("clone POPULATION %d > deme POPULATION %d\n", clone_ints[POPULATION][chosen_clone], deme_ints[POPULATION][parent_deme_num]);
		fprintf(error_log, "iterations %d, event_type %d, chosen_clone %d, parent_deme_num %d\n", iterations, event_type, chosen_clone, parent_deme_num);
		fprintf(error_log, "clone POPULATION %d > deme POPULATION %d\n", clone_ints[POPULATION][chosen_clone], deme_ints[POPULATION][parent_deme_num]);
		exit_code = 1;
	}
}

void check_normal_pops(int deme_index)
{
	if(deme_ints[NORMAL_CELLS][deme_index] > DBL_MAX) {
		printf("Normal cell count is infinite\n");
		fprintf(error_log, "Normal cell count is infinite\n");
		exit_code = 1;
	}
	if(deme_ints[NORMAL_CELLS][deme_index] < 0) {
		printf("Normal pop < 0\n");
		fprintf(error_log, "Normal pop < 0\n");
		exit_code = 1;
	}
}

void check_genotype_counts(int num_matrix_cols, int num_empty_cols, int num_extinct_genotypes, int num_driver_matrix_cols, int num_empty_driver_cols, int num_extinct_driver_genotypes, 
	int cell_type, int event_type, int chosen_deme)
{
	if(mu_passenger == 0 && passenger_pop_threshold >= max_pop && num_matrix_cols - num_empty_cols - num_extinct_genotypes != num_driver_matrix_cols - num_empty_driver_cols - num_extinct_driver_genotypes) {
		printf("Inequality! Iterations = %d\n; cell type = %d; event type = %d\n", iterations, cell_type, event_type);
		printf("chosen_deme = %d; POPULATION = %d\n", chosen_deme, deme_ints[POPULATION][chosen_deme]);
		printf("num_matrix_cols - num_empty_cols - num_extinct_genotypes = %d - %d - %d, num_driver_matrix_cols - num_empty_driver_cols - num_extinct_driver_genotypes = %d - %d - %d\n", 
			num_matrix_cols, num_empty_cols, num_extinct_genotypes, num_driver_matrix_cols, num_empty_driver_cols, num_extinct_driver_genotypes);
		fprintf(error_log, "Inequality! Iterations = %d\n; cell type = %d; event type = %d\n", iterations, cell_type, event_type);
		exit_code = 1;
	}
}

/////////////////////// generic bintree functions:

// reset bintree sums for all elements of layer 0 by summing over array elements:
void set_bintree_sums_layer0(double *bintree, double *array, int num_array_elements)
{
	int i, j, start_index, end_index;

	for(i = 0; i < ceil((float)num_array_elements/SET_SIZE); i++) { // loop through the bintree elements
		bintree[i] = 0;
		start_index = SET_SIZE * i; // index of the first deme in the bintree element
		end_index = MIN(start_index + SET_SIZE - 1, num_array_elements - 1); // index of the last deme in the bintree element
		// sum over constituent demes:
		for(j = start_index; j <= end_index; j++) bintree[i] += array[j];
	}
}

// reset bintree sums for all elements of layer 0 by summing over array elements:
void set_bintree_sums_layer0_int(int *bintree, int *array, int num_array_elements)
{
	int i, j, start_index, end_index;

	for(i = 0; i < ceil((float)num_array_elements/SET_SIZE); i++) { // loop through the bintree elements
		bintree[i] = 0;
		start_index = SET_SIZE * i; // index of the first deme in the bintree element
		end_index = MIN(start_index + SET_SIZE - 1, num_array_elements - 1); // index of the last deme in the bintree element
		// sum over constituent demes:
		for(j = start_index; j <= end_index; j++) bintree[i] += array[j];
	}
}

// reset bintree sums for all elements of layers > 0 by summing over array elements:
void set_bintree_sums_subsequent_layers(double *bintree, int num_array_elements, int max_layer_needed, int array_length)
{
	int i, j, layer, start_index, end_index, bintree_index;
	int num_elements_previous_layer, min_bintree_index_this_layer, min_bintree_index_previous_layer = 0;

	num_elements_previous_layer = num_array_elements; // number of bintree elements in the layer below layer 0
	for(layer = 1; layer <= max_layer_needed; layer++) {
		num_elements_previous_layer = ceil((float)num_elements_previous_layer / SET_SIZE); // number of bintree elements in previous layer
		if(array_length == max_demes) min_bintree_index_this_layer = min_deme_bintree_index_by_layer[layer];
		else if(array_length == max_clones_per_deme) min_bintree_index_this_layer = min_clone_bintree_index_by_layer[layer];
		else min_bintree_index_this_layer = get_bintree_index(layer, 0, array_length); // index of the very first bintree element from the current layer
		for(i = 0; i < ceil((float)num_elements_previous_layer/SET_SIZE); i++) { // loop through the bintree elements in the current layer
			bintree_index = min_bintree_index_this_layer + i;
			bintree[bintree_index] = 0;
			start_index = min_bintree_index_previous_layer + SET_SIZE * i; // index of the first bintree element from the previous layer that is in the current bintree element
			end_index = MIN(start_index + SET_SIZE - 1, min_bintree_index_previous_layer + num_elements_previous_layer - 1);
			// sum over constituent bintree elements from previous layer:
			for(j = start_index; j <= end_index; j++) bintree[bintree_index] += bintree[j];
		}
		min_bintree_index_previous_layer = min_bintree_index_this_layer; // index of the very first bintree element from the previous layer
	}
}

// reset bintree sums for all elements of layers > 0 by summing over array elements:
void set_bintree_sums_subsequent_layers_int(int *bintree, int num_array_elements, int max_layer_needed, int array_length)
{
	int i, j, layer, start_index, end_index, bintree_index;
	int num_elements_previous_layer, min_bintree_index_this_layer, min_bintree_index_previous_layer = 0;

	num_elements_previous_layer = num_array_elements; // number of bintree elements in the layer below layer 0
	for(layer = 1; layer <= max_layer_needed; layer++) {
		num_elements_previous_layer = ceil((float)num_elements_previous_layer / SET_SIZE); // number of bintree elements in previous layer
		if(array_length == max_demes) min_bintree_index_this_layer = min_deme_bintree_index_by_layer[layer];
		else if(array_length == max_clones_per_deme) min_bintree_index_this_layer = min_clone_bintree_index_by_layer[layer];
		else min_bintree_index_this_layer = get_bintree_index(layer, 0, array_length); // index of the very first bintree element from the current layer
		for(i = 0; i < ceil((float)num_elements_previous_layer/SET_SIZE); i++) { // loop through the bintree elements in the current layer
			bintree_index = min_bintree_index_this_layer + i;
			bintree[bintree_index] = 0;
			start_index = min_bintree_index_previous_layer + SET_SIZE * i; // index of the first bintree element from the previous layer that is in the current bintree element
			end_index = MIN(start_index + SET_SIZE - 1, min_bintree_index_previous_layer + num_elements_previous_layer - 1);
			// sum over constituent bintree elements from previous layer:
			for(j = start_index; j <= end_index; j++) bintree[bintree_index] += bintree[j];
		}
		min_bintree_index_previous_layer = min_bintree_index_this_layer; // index of the very first bintree element from the previous layer
	}
}

// update bintree across all layers:
void update_bintree_layers(double summand, double *bintree, int array_index, int max_layer_needed, int array_length)
{
	int layer;
	int bintree_index;

	for(layer = 0; layer <= max_layer_needed; layer++) {
		if(array_length == max_demes) bintree_index = get_deme_bintree_index(layer, array_index);
		else if(array_length == max_clones_per_deme) bintree_index = get_clone_bintree_index(layer, array_index);
		else bintree_index = get_bintree_index(layer, array_index, array_length);
		bintree[bintree_index] += summand; // add or subtract from the bintree element
	}
}

// update bintree across all layers:
void update_bintree_layers_int(int summand, int *bintree, int array_index, int max_layer_needed, int array_length)
{
	int layer;
	int bintree_index;

	for(layer = 0; layer <= max_layer_needed; layer++) {
		if(array_length == max_demes) bintree_index = get_deme_bintree_index(layer, array_index);
		else if(array_length == max_clones_per_deme) bintree_index = get_clone_bintree_index(layer, array_index);
		else bintree_index = get_bintree_index(layer, array_index, array_length);
		bintree[bintree_index] += summand; // add or subtract from the bintree element
	}
}

// add new bintree layer and/or element if needed:
void extend_bintree(double *bintree, int max_array_index, int max_layer_needed, int array_length)
{
	int i, layer, bintree_index, start_index;
	int previous_array_size = max_array_index;

	// if a new bintree element is needed then initialise it at zero:
	for(layer = 0; layer <= max_layer_needed; layer++) {
		if(array_length == max_demes) bintree_index = get_deme_bintree_index(layer, max_array_index);
		else if(array_length == max_clones_per_deme) bintree_index = get_clone_bintree_index(layer, max_array_index);
		else bintree_index = get_bintree_index(layer, max_array_index, array_length);
		if(max_array_index % powers_list[layer + 1] == 0) bintree[bintree_index] = 0;
	}
	
	// if a new bintree layer is needed then set the value of the first bintree element in this layer:
	if(max_layer_needed > max_layer_needed_array[previous_array_size]) {
		if(array_length == max_demes) {
			bintree_index = min_deme_bintree_index_by_layer[max_layer_needed];
			start_index = min_deme_bintree_index_by_layer[max_layer_needed - 1];
		}
		else if(array_length == max_clones_per_deme) {
			bintree_index = min_clone_bintree_index_by_layer[max_layer_needed];
			start_index = min_clone_bintree_index_by_layer[max_layer_needed - 1];
		}
		else {
			bintree_index = get_bintree_index(max_layer_needed, 0, array_length);
			start_index = get_bintree_index(max_layer_needed - 1, 0, array_length); // index of the very first bintree element from the previous layer
		}
		bintree[bintree_index] = 0;
		// initialise the new bintree element with the sum of sums from the previous layer:
		for(i = start_index; i < start_index + SET_SIZE; i++) bintree[bintree_index] += bintree[i];
	}
}

// add new bintree layer and/or element if needed:
void extend_bintree_int(int *bintree, int max_array_index, int max_layer_needed, int array_length)
{
	int i, layer, bintree_index, start_index;
	int previous_array_size = max_array_index;

	// if a new bintree element is needed then initialise it at zero:
	for(layer = 0; layer <= max_layer_needed; layer++) {
		if(array_length == max_demes) bintree_index = get_deme_bintree_index(layer, max_array_index);
		else if(array_length == max_clones_per_deme) bintree_index = get_clone_bintree_index(layer, max_array_index);
		else bintree_index = get_bintree_index(layer, max_array_index, array_length);
		if(max_array_index % powers_list[layer + 1] == 0) bintree[bintree_index] = 0;
	}
	
	// if a new bintree layer is needed then set the value of the first bintree element in this layer:
	if(max_layer_needed > max_layer_needed_array[previous_array_size]) {
		if(array_length == max_demes) {
			bintree_index = min_deme_bintree_index_by_layer[max_layer_needed];
			start_index = min_deme_bintree_index_by_layer[max_layer_needed - 1];
		}
		else if(array_length == max_clones_per_deme) {
			bintree_index = min_clone_bintree_index_by_layer[max_layer_needed];
			start_index = min_clone_bintree_index_by_layer[max_layer_needed - 1];
		}
		else {
			bintree_index = get_bintree_index(max_layer_needed, 0, array_length);
			start_index = get_bintree_index(max_layer_needed - 1, 0, array_length); // index of the very first bintree element from the previous layer
		}
		bintree[bintree_index] = 0;
		// initialise the new bintree element with the sum of sums from the previous layer:
		for(i = start_index; i < start_index + SET_SIZE; i++) bintree[bintree_index] += bintree[i];
	}
}

// given an array index, get the index of its enclosing bintree element at a particular layer (from scratch):
int get_bintree_index(int layer, int array_index, int array_size)
{
	int res = 0;
	int i;

	for(i = 1; i <= layer; i++) res += (array_size - 1) / powers_list[i] + 1; // integer division
	res = res + array_index / powers_list[layer + 1]; // integer division

	return res;
}

// given an deme index, get the index of its enclosing bintree element at a particular layer (using lookup table):
int get_deme_bintree_index(int layer, int array_index)
{
	return min_deme_bintree_index_by_layer[layer] + array_index / powers_list[layer + 1]; // integer division
}

// given a clone index, get the index of its enclosing bintree element at a particular layer (using lookup table):
int get_clone_bintree_index(int layer, int array_index)
{
	return min_clone_bintree_index_by_layer[layer] + array_index / powers_list[layer + 1]; // integer division
}

// calculate the highest bintree layer necessary, based on the current number of array elements:
int get_max_layer_needed(int n)
{
	int res = -2;
	int out = n - 1;

	do { // n < SET_SIZE^2 => res = 0
		out = out / SET_SIZE; // integer division
		res++;
	} while(out > 0 || res < 0);

	return(res);
}

/////////////////////// calculate diversity metrics (top level):

// calculate diversity metrics based on either genotype or driver genotype:
void get_diversity_metrics(int rao, float *diversity, float *edge_diversity, float *alpha_diversity, float *within_deme_diversity, int num_demes, int num_cells, int num_matrix_cols, int num_clones, 
	int *clone_genotype, int **either_matrix, int dmax, int *populations, long *idum, float gens_elapsed, int* position_in_edge_list, int* demes_at_edge, int* sides_at_edge, int* genotype_edge_pop, int is_drivers_only)
{	
	int i;

	// calculate total diversity:
	*diversity = calculate_diversity(num_matrix_cols, num_cells, populations, either_matrix, dmax, rao);
	
	if(!well_mixed) {
		// calculate edge diversity:
		*edge_diversity = calculate_edge_diversity(num_matrix_cols, num_cells, either_matrix, dmax, num_demes, num_clones, clone_genotype, idum, rao, position_in_edge_list, demes_at_edge, sides_at_edge, genotype_edge_pop);

		// calculate within-deme diversity for each deme:
		for(i=0; i<num_demes; i++) within_deme_diversity[i] = calculate_within_deme_diversity(i, num_matrix_cols, num_clones, clone_genotype, either_matrix, dmax, rao, is_drivers_only, genotype_edge_pop);

		// calculate alpha diversity:
		*alpha_diversity = 0;
		for(i=0; i<num_demes; i++) *alpha_diversity += within_deme_diversity[i] * deme_ints[POPULATION][i];
		*alpha_diversity /= num_cells;
	}
	else {
		*edge_diversity = *diversity;
		*alpha_diversity = *diversity;
	}
}

// calculate diversity and record phylogenies of sampled cells:
void get_biopsy_data(int rao, float *depth_array, int samples, int biopsy_size_per_sample, int num_demes, int num_cells, int num_matrix_cols, int num_clones, int *clone_genotype, int **either_matrix, 
	int dmax, long *idum, FILE *sample_size_log, float gens_elapsed, FILE *output_phylo_of_sample, int calculate_sample_diversities, int record_phylogenies, float centre_X, float centre_Y, 
	int **geno_or_driver_ints, float **geno_or_driver_floats, int *at_edge)
{
	int i, depth, geno_num;
	int genotype_populations_in_sample[num_matrix_cols];
	int sampled_cells;

	// record phylogeny (and optionally calculate diversity) based on biopsy samples at various depths:
	if(!well_mixed) for(depth=0; depth<=10; depth++) {
		find_sample_genotype_pops(genotype_populations_in_sample, &sampled_cells, num_matrix_cols, num_cells, either_matrix, dmax, num_demes, num_clones, depth, centre_X, centre_Y, 
			clone_genotype, samples, biopsy_size_per_sample, idum, sample_size_log, gens_elapsed, rao, at_edge);
		if(exit_code != 0) return;
		if(calculate_sample_diversities) depth_array[depth] = calculate_diversity(num_matrix_cols, sampled_cells, genotype_populations_in_sample, either_matrix, dmax, rao);
		if(record_phylogenies) write_output_phylo(output_phylo_of_sample, num_matrix_cols, gens_elapsed, genotype_populations_in_sample, geno_or_driver_ints, geno_or_driver_floats, 
			samples, biopsy_size_per_sample, depth, num_cells);
	}

	// record phylogeny (and optionally calculate diversity) based on randomly sampled cells:
	sampled_cells = MIN(samples * biopsy_size_per_sample, num_cells);
	for(i = 0; i < num_matrix_cols; i++) genotype_populations_in_sample[i] = 0;
	for(i = 0; i < sampled_cells; i++) {
		geno_num = weighted_random_ints(geno_or_driver_ints[POPULATION], idum, 0, num_matrix_cols); // randomly pick a genotype (probabilities proportional to population sizes)
		genotype_populations_in_sample[geno_num]++;
	}
	if(calculate_sample_diversities) depth_array[11] = calculate_diversity(num_matrix_cols, sampled_cells, genotype_populations_in_sample, either_matrix, dmax, rao);
	if(record_phylogenies) write_output_phylo(output_phylo_of_sample, num_matrix_cols, gens_elapsed, genotype_populations_in_sample, geno_or_driver_ints, geno_or_driver_floats, 
		samples, biopsy_size_per_sample, -1, num_cells);
}

// calculate mean and variance of numbers of mutations:
void calculate_mutation_metrics(float *mean_num_passengers, float *mean_num_drivers, float *var_num_passengers, float *var_num_drivers, int *driver_counts, int num_matrix_cols, int num_cells)
{
	int i;
	int pop1, sum_pop, num_drivers;
	float diff1;

	*mean_num_passengers = 0;
	*mean_num_drivers = 0;
	for(i = 0; i <= MAX_DRIVERS_TO_COUNT; i++) driver_counts[i] = 0;
	sum_pop = 0;
	for(i=0; i<num_matrix_cols; i++) { // loop over all genotypes or driver genotypes
		pop1 = genotype_ints[POPULATION][i];
		num_drivers = genotype_ints[NUM_DRIVER_MUTATIONS][i];
		*mean_num_passengers += pop1 * genotype_ints[NUM_PASSENGER_MUTATIONS][i];
		*mean_num_drivers += pop1 * (num_drivers + genotype_ints[NUM_MIGRATION_MUTATIONS][i]);
		sum_pop += pop1;
		driver_counts[(int)MIN(num_drivers, MAX_DRIVERS_TO_COUNT)] += pop1;
	}
	*mean_num_passengers /= num_cells;
	*mean_num_drivers /= num_cells;

	// error check:
	if(sum_pop != num_cells) {
		printf("Genotype populations don't sum to num_cells (%d versus %d) at iterations %d\n", sum_pop, num_cells, iterations);
		fprintf(error_log, "Genotype populations don't sum to num_cells (%d versus %d) at iterations %d\n", sum_pop, num_cells, iterations);
	}

	*var_num_passengers = 0;
	*var_num_drivers = 0;
	for(i=0; i<num_matrix_cols; i++) { // loop over all genotypes or driver genotypes
		pop1 = genotype_ints[POPULATION][i];
		diff1 = genotype_ints[NUM_PASSENGER_MUTATIONS][i] - *mean_num_passengers;
		*var_num_passengers += pop1 * diff1 * diff1;
		diff1 = genotype_ints[NUM_DRIVER_MUTATIONS][i] + genotype_ints[NUM_MIGRATION_MUTATIONS][i] - *mean_num_drivers;
		*var_num_drivers += pop1 * diff1 * diff1;
	}
	*var_num_passengers /= num_cells;
	*var_num_drivers /= num_cells;
}

// calculate mean and variance of birth or migration rate:
float calculate_variance_of_rate(int num_matrix_cols, int num_cells, double sum_of_rates, float *list_of_rates)
{
	int i;
	int pop1;
	float diff1;
	float var;

	var = 0;
	for(i=0; i<num_matrix_cols; i++) {
		pop1 = genotype_ints[POPULATION][i];
		diff1 = list_of_rates[i] - sum_of_rates / num_cells;
		var += pop1 * diff1 * diff1;
	}

	return var / num_cells;
}

// find and record each genotype's parent, next sister, first daughter and last daughter:
void get_relatives(int num_matrix_cols, int next_genotype_id, int **geno_or_driver_ints)
{
	int i;
	int parent_geno_num, sister_geno_num;
	int id, parent_id;
	int* id_to_index;

	id_to_index = (int *) malloc(next_genotype_id * sizeof *id_to_index);
	if(id_to_index == NULL) {printf("Memory problem!\n"); exit(1);}

	for(i = 0; i < num_matrix_cols; i++) {
		id = geno_or_driver_ints[IDENTITY][i];
		id_to_index[id] = i;
	}

	for(i = 0; i < num_matrix_cols; i++) {
		genotype_relatives[PARENT_INDEX][i] = -1;
		genotype_relatives[NEXT_SISTER][i] = -1;
		genotype_relatives[FIRST_DAUGHTER][i] = -1;
		genotype_relatives[LAST_DAUGHTER][i] = -1;
	}

	for(i = 1; i < num_matrix_cols; i++) { // start at first genotype excluding the original genotype
		parent_id = geno_or_driver_ints[PARENT][i];
		parent_geno_num = id_to_index[parent_id];
		genotype_relatives[PARENT_INDEX][i] = parent_geno_num;
		if(genotype_relatives[FIRST_DAUGHTER][parent_geno_num] < 0) {
			genotype_relatives[FIRST_DAUGHTER][parent_geno_num] = i;
			genotype_relatives[LAST_DAUGHTER][parent_geno_num] = i;
		}
		else {
			sister_geno_num = genotype_relatives[LAST_DAUGHTER][parent_geno_num];
			genotype_relatives[NEXT_SISTER][sister_geno_num] = i;
			genotype_relatives[LAST_DAUGHTER][parent_geno_num] = i;
		}
	}
	free(id_to_index);
}

void get_allele_frequencies(int num_matrix_cols, int *allele_count, int **geno_or_driver_ints)
{
	bool upped = 0; // whether the last move was up
	int now = 0; // now genotype (initially the common ancestor)
	int prev; // previous genotype
	int moves_count = 0; // number of moves made so far (for error check)

	allele_count[0] = geno_or_driver_ints[POPULATION][0]; // record common ancestor's allele count

	while(1) {
		while(!upped) { // if most recent move wasn't up
			prev = now;
			now = move_down(now); // move to first daughter (if there is one)
			moves_count++;
			upped = 0;
			if(now == prev) { // if now genotype has no daughters
				if(now != 0) allele_count[move_up(now)] += allele_count[now]; // add now genotype's population to its parent's allele count (second visit)
				break; // stop trying to move down
			}
			else allele_count[now] = geno_or_driver_ints[POPULATION][now]; // record now genotype's allele count (first visit)
		}
		if(move_right(now) != now) { // if there is a next eldest sister
			prev = now;
			now = move_right(now); // move to next eldest sister
			allele_count[now] = geno_or_driver_ints[POPULATION][now]; // record now genotype's allele count (first visit)
			moves_count++;
			upped = 0;
		}
		else if(now != 0) { // if the now genotype has a parent (i.e. if it's not the common ancestor)
			prev = now;
			now = move_up(now); // move to the parent
			if(now != 0) allele_count[move_up(now)] += allele_count[now]; // add now genotype's population to its parent's allele count (second visit)
			moves_count++;
			upped = 1; // record that the most recent move was up
		}
		if(now == 0) break; // returned to the first genotype (common ancestor)
		if(moves_count > 2 * num_matrix_cols) { // error check: every genotype should be visited exactly twice
			printf("Allele freq finder is stuck in a loop\n");
			fprintf(error_log, "Allele freq finder is stuck in a loop\n");
			return;
		}
	}

	if(moves_count != 2 * num_matrix_cols - 1) {
		printf("Adjacency matrix seems to be bipartite (moves_count = %d != %d)\n", moves_count, 2 * num_matrix_cols - 1);
		fprintf(error_log, "Adjacency matrix seems to be bipartite (moves_count = %d != %d)\n", moves_count, 2 * num_matrix_cols - 1);
	}
}

void get_frequency_table(int input_length, int *count, int **freq_table, int *output_length)
{
	int i, j;
	int matched;

	// set up the first few entries to speed up subsequent lookups:
	for(i = 0; i < 10; i++) {
		freq_table[0][i] = i + 1;
		freq_table[1][i] = 0;
	}
	*output_length = 10;

	for(i = 0; i < input_length; i++) if(count[i] > 0) {
		matched = 0;
		for(j = 0; j < *output_length; j++) {
			if(freq_table[0][j] == count[i]) {
				freq_table[1][j]++;
				matched = 1;
				break;
			}
		}
		if(matched == 0) {
			freq_table[0][*output_length] = count[i];
			freq_table[1][*output_length] = 1;
			*output_length = *output_length + 1;
		}
		if(*output_length > max_distinct_allele_freqs) {
			printf("output_length (%d) > max_distinct_allele_freqs (%d)\n", *output_length, max_distinct_allele_freqs);
			fprintf(error_log, "output_length (%d) > max_distinct_allele_freqs (%d)\n", *output_length, max_distinct_allele_freqs);
			exit_code = 1;
			return;
		}
	}
}

/////////////////////// calculate diversity metrics (lower level):

// find centre of gravity of the tumour, based on deme cell counts:
float centre_of_gravity(int direction, int num_cells)
{
	// direction 1 = horizontal, 2 = vertical
	int x, y;
	int deme_num;
	float sum = 0;
	float coord;

	for(x=0; x<dim_grid; x++) for(y=0; y<dim_grid; y++) {
		deme_num = grid[x][y];
		if(deme_num != EMPTY) {
			if(direction == 1) coord = deme_ints[XCOORD][deme_num];
			else coord = deme_ints[YCOORD][deme_num];
			sum += deme_ints[POPULATION][deme_num] * coord;
		}
	}

	return sum / num_cells;
}

// calculate diversity of the cancer cell population, based on either genotypes or driver genotypes:
float calculate_diversity(int num_matrix_cols, int N, int *populations, int **either_matrix, int dmax, int rao)
{
	int i, j;
	int pop1, pop2;
	double diversity;
	long pops_sum = 0;
	int distance;

	if(rao) for(i=0; i<num_matrix_cols; i++) for(j=0; j<i; j++) { // rao = 1 => use Rao's quadratic diversity
		pop1 = populations[i];
		pop2 = populations[j];
		distance = MIN(either_matrix[j][i], dmax);
		pops_sum += 2 * (long)pop1 * (long)pop2 * (long)distance;
		// factor of two because sum here is over only half of the matrix, whereas it should be over the whole matrix
	}
	else for(i=0; i<num_matrix_cols; i++) { // otherwise use Simpson's index
		pop1 = populations[i];
		pops_sum += (long)pop1 * (long)pop1;
	}

	diversity = (double)((double)pops_sum / N) / N;

	if(rao) return 1.0 / (1.0 - diversity / dmax);
	else return 1.0 / diversity;
}

float calculate_edge_diversity(int num_matrix_cols, int num_cells, int **either_matrix, int dmax, int num_demes, int num_clones, int *clone_genotype, long *idum, int rao, 
	int* position_in_edge_list, int* demes_at_edge, int* sides_at_edge, int* genotype_edge_pop)
{
	int i, j;
	int candidate;
	int num_demes_at_edge = 0;
	int edge_pop = 0;
	int loop;
	int geno_num, deme_num, clone_index;
	int x, y;
	float sample_size_per_deme_edge = sqrt(K);
	int sample_size, clone_sample_size;
	int deme_pop;
	int clone_pop;
	int num_draws;
	int num_already_sampled_from_deme, previous_clone_pops;
	float result;
	int buff;

	// initiate arrays:
	for(i=0; i<num_demes; i++) {
		position_in_edge_list[i] = -1; // yet to be assigned a position in the list
		sides_at_edge[i] = 0; // initiate counters of deme sides on the edge
	}

	// identify the (occupied) demes at the edge and record their indices in the demes_at_edge array;
	// also record sides_at_edge, and position_in_edge_list:
	for(loop=0; loop<4; loop++) for(i=0; i<dim_grid; i++) for(j=0; j<dim_grid; j++) {
		if(loop == 0) {x = i; y = j;} // on the grid images, this is sweeping from left to right
		else if(loop == 1) {x = i; y = dim_grid - j - 1;} // sweeping from right to left
		else if(loop == 2) {x = j; y = i;} // sweeping from bottom to top
		else {x = dim_grid - j - 1; y = i;} // sweeping from top to bottom

		candidate = grid[x][y];

		if(candidate != EMPTY) if(deme_ints[POPULATION][candidate] > 0) {
			if(position_in_edge_list[candidate] == -1) { // first time deme is found to be on the edge
				demes_at_edge[num_demes_at_edge] = candidate; // indices of demes that are on the edge
				sides_at_edge[num_demes_at_edge]++; // number of deme sides that are on the edge
				position_in_edge_list[candidate] = num_demes_at_edge; // position of deme in demes_at_edge and sides_at_edge arrays
				num_demes_at_edge++; // number of demes that are on the edge
			} 
			else sides_at_edge[position_in_edge_list[candidate]]++;
			break;
		}
	}

	// initialise array for recording genotype populations:
	for(i=0; i<num_matrix_cols; i++) genotype_edge_pop[i] = 0;

	for(i=0; i<num_demes_at_edge; i++) { // loop over list of demes at the edge
		deme_num = demes_at_edge[i];
		deme_pop = deme_ints[POPULATION][deme_num]; // deme population
		// sample size for the deme is proportional to the number of deme sides that are on the edge,
		// except that it cannot exceed the deme population size:
		buff = stochastic_round(sides_at_edge[i] * sample_size_per_deme_edge, idum);
		sample_size = MIN(buff, deme_pop);
		num_already_sampled_from_deme = 0;
		previous_clone_pops = 0;
		for(j=0; j<deme_ints[NUM_CLONES_IN_DEME][deme_num]; j++) { // loop over list of clones in deme
			clone_index = clones_list_in_deme[deme_num][j]; // clone index
			clone_pop = clone_ints[POPULATION][clone_index]; // clone population
			geno_num = clone_genotype[clone_index]; // genotype index of the clone
			
			// clone_sample_size is drawn from a hypergeometric distribution (i.e. sampling without replacement):
			num_draws = MAX(sample_size - num_already_sampled_from_deme, 0);
			clone_sample_size = hypergeometric(idum, clone_pop, deme_pop - previous_clone_pops - clone_pop, num_draws);

			num_already_sampled_from_deme += clone_sample_size;
			previous_clone_pops += clone_pop;
			genotype_edge_pop[geno_num] += clone_sample_size; // add clone population to the sample size of the genotype to include in diversity calculation
			edge_pop += clone_sample_size; // add clone population to total edge population
		}
	}

	// calculate diversity, using the genotype edge subpopulations:
	result = calculate_diversity(num_matrix_cols, edge_pop, genotype_edge_pop, either_matrix, dmax, rao);

	return result;
}

float calculate_within_deme_diversity(int deme_index, int num_matrix_cols, int num_clones, int *clone_genotype, int **either_matrix, int dmax, int rao, int is_drivers_only, int* genotype_edge_pop)
{
	int i, j;
	int pop1, pop2;
	double diversity;
	long pops_sum = 0;
	int cells_in_deme = deme_ints[POPULATION][deme_index];
	int distance;
	int geno_num1, geno_num2;
	int clone_num1, clone_num2;
	float result;

	if(cells_in_deme == 0) return 0;

	// calculate the diversity between these clones:
	// rao = 1 => use Rao's quadratic diversity
	if(rao) for(i=0; i<deme_ints[NUM_CLONES_IN_DEME][deme_index]; i++) for(j=0; j<i; j++) {
		clone_num1 = clones_list_in_deme[deme_index][i];
		clone_num2 = clones_list_in_deme[deme_index][j];
		pop1 = clone_ints[POPULATION][clone_num1];
		pop2 = clone_ints[POPULATION][clone_num2];

		// need to use the maxi and mini functions here to ensure that geno_num1 < geno_num2
		// (because only the top half of the matrix is filled)
		geno_num1 = MAX(clone_genotype[clone_num1], clone_genotype[clone_num2]); // either genotype index or driver genotype index
		geno_num2 = MIN(clone_genotype[clone_num1], clone_genotype[clone_num2]);
		distance = MIN(either_matrix[geno_num2][geno_num1], dmax);

		pops_sum += 2 * (long)pop1 * (long)pop2 * (long)distance;
	}
	else if(!is_drivers_only) for(i=0; i<deme_ints[NUM_CLONES_IN_DEME][deme_index]; i++) { // Simpson's index for all genotypes
		clone_num1 = clones_list_in_deme[deme_index][i];
		pop1 = clone_ints[POPULATION][clone_num1];
		pops_sum += (long)pop1 * (long)pop1;
	}
	else { // Simpson's index for driver genotypes (which requires making a list of driver genotype population sizes in the deme)
		// initialise array for recording genotype populations:
		for(i=0; i<num_matrix_cols; i++) genotype_edge_pop[i] = 0;
		
		for(i=0; i<deme_ints[NUM_CLONES_IN_DEME][deme_index]; i++) {
			clone_num1 = clones_list_in_deme[deme_index][i];
			pop1 = clone_ints[POPULATION][clone_num1];
			geno_num1 = clone_genotype[clone_num1]; // genotype index of the clone

			genotype_edge_pop[geno_num1] += pop1; // add clone population to the sample size of the genotype to include in diversity calculation
		}
		// calculate diversity, using the genotype within-deme subpopulations:
		result = calculate_diversity(num_matrix_cols, cells_in_deme, genotype_edge_pop, either_matrix, dmax, 0);
		return(result);
	}

	diversity = (double)((double)pops_sum / cells_in_deme) / cells_in_deme;

	if(rao) return 1.0 / (1.0 - diversity / dmax);
	else return 1.0 / diversity;
}

void find_sample_genotype_pops(int *genotype_populations_in_sample, int *sampled_cells, int num_matrix_cols, int num_cells, int **either_matrix, int dmax, int num_demes, int num_clones, 
	int depth, float centre_X, float centre_Y, int *clone_genotype, int num_directions, int biopsy_size_per_sample, long *idum, FILE *sample_size_log, float gens_elapsed, 
	int rao, int *at_edge)
{
	int i, j;
	int candidate;
	float coordX[num_directions], coordY[num_directions];
	int loop;
	int geno_num;
	int x, y;
	float distX, distY, dist;
	float R = sqrt(biopsy_size_per_sample / (K * PI)); // radius of the sample disc
	float r = 0.5; // deme radius
	float p; // proportion of deme within the sample disc
	int sample_size, clone_sample_size;
	int deme_pop, clone_pop;
	int deme_num, clone_index;
	int num_already_sampled_from_deme, previous_clone_pops, num_draws;
	int limx1, limx2, limy1, limy2, int_buff;

	*sampled_cells = 0;

	for(i=0; i<num_demes; i++) at_edge[i] = 0;

	// identify the coordinates of (occupied) demes at the required depth and record them in the coordX and coordY arrays
	for(loop=0; loop<num_directions; loop++) {

		for(i=0; i<dim_grid; i++) {
			if(loop == 0) {x = ROUNDIT(centre_X); y = i;} // on the grid images, this is sweeping from left to right
			else if(loop == 1) {x = ROUNDIT(centre_X); y = dim_grid - i - 1;} // sweeping from right to left
			else if(loop == 2) {x = i; y = ROUNDIT(centre_Y);} // sweeping from bottom to top
			else {x = dim_grid - i - 1; y = ROUNDIT(centre_Y);} // sweeping from top to bottom

			candidate = grid[x][y];

			if(candidate != EMPTY) if(deme_ints[POPULATION][candidate] > 0) {
				if(loop == 0 || loop == 1) {
					coordX[loop] = x;
					coordY[loop] = y + (centre_Y - y) * (float)depth / 10;
				}
				else {
					coordX[loop] = x + (centre_X - x) * (float)depth / 10;
					coordY[loop] = y;
				}
				break;
			}

			// error check:
			if(i == dim_grid - 1) {
				printf("Warning: couldn't find any non-empty grid squares; centre_X %f, centre_Y %f\n", centre_X, centre_Y);
				fprintf(error_log, "Warning: couldn't find any non-empty grid squares; centre_X %f, centre_Y %f\n", centre_X, centre_Y);
				coordX[loop] = centre_X;
				coordY[loop] = centre_Y;
			}
		}
	}

	for(i=0; i<num_matrix_cols; i++) genotype_populations_in_sample[i] = 0;

	// loop through nearby demes, counting how many cells of each genotype are in each deme:
	for(loop=0; loop<num_directions; loop++) {

		limx1 = ROUNDIT(coordX[loop] - R);
		limx2 = ROUNDIT(coordX[loop] + R);
		for(x = MAX(limx1, 0); x <= MIN(limx2, dim_grid - 1); x++) {
			limy1 = ROUNDIT(coordY[loop] - R);
			limy2 = ROUNDIT(coordY[loop] + R);
			for(y = MAX(limy1, 0); y <= MIN(limy2, dim_grid - 1); y++) {			
				deme_num = grid[x][y];
				if(deme_num != EMPTY) {
					distX = fabs(x - coordX[loop]);
					distY = fabs(y - coordY[loop]);
					dist = sqrt(distX*distX + distY*distY);
					
					// sample from deme proportional to the overlap of the disc and the deme:
					if(dist > R + r) p = 0; // deme is wholly outside neighbourhood
					else if(R >= dist + r) p = 1; // deme is wholly inside neighbourhood
					else if(r >= dist + R) p = (float)biopsy_size_per_sample / K; // neighbourhood is wholly inside deme
					else p = r*r*acos((dist*dist + r*r - R*R)/(2*dist*r)) + R*R*acos((dist*dist + R*R -r*r)/(2*dist*R)) - 1.0/2*sqrt((-dist+r+R)*(dist+r-R)*(dist-r+R)*(dist+r+R));
					// see http://mathworld.wolfram.com/Circle-CircleIntersection.html

					deme_pop = deme_ints[POPULATION][deme_num]; // deme population
					// sample size for the deme is proportional to the overlap,
					// except that it cannot exceed the deme population size:
					int_buff = stochastic_round(p * K, idum);
					sample_size = MIN(int_buff, deme_pop);
					num_already_sampled_from_deme = 0;
					previous_clone_pops = 0;
					for(j=0; j<deme_ints[NUM_CLONES_IN_DEME][deme_num]; j++) { // loop over list of clones in deme
						clone_index = clones_list_in_deme[deme_num][j]; // clone index
						clone_pop = clone_ints[POPULATION][clone_index]; // clone population
						geno_num = clone_genotype[clone_index]; // genotype index of the clone
						
						// clone_sample_size is drawn from a hypergeometric distribution (i.e. sampling without replacement):
						num_draws = MAX(sample_size - num_already_sampled_from_deme, 0);
						clone_sample_size = hypergeometric(idum, clone_pop, deme_pop - previous_clone_pops - clone_pop, num_draws);

						num_already_sampled_from_deme += clone_sample_size;
						previous_clone_pops += clone_pop;
						genotype_populations_in_sample[geno_num] += clone_sample_size; // add clone population to the sample size of the genotype to include in diversity calculation
						*sampled_cells += clone_sample_size; // add clone population to total sampled population
					}
				}
			}
		}
	}

	fprintf(sample_size_log, "%f\t%d\t%d\t%d\t%d\n", gens_elapsed, num_directions, depth, num_directions * biopsy_size_per_sample, *sampled_cells);
}

int move_down(int start_index)
{
	if(genotype_relatives[FIRST_DAUGHTER][start_index] >= 0) return genotype_relatives[FIRST_DAUGHTER][start_index];
	else return start_index;
}

int move_right(int start_index)
{
	if(genotype_relatives[NEXT_SISTER][start_index] >= 0) return genotype_relatives[NEXT_SISTER][start_index];
	else return start_index;
}

int move_up(int start_index)
{
	if(genotype_relatives[PARENT_INDEX][start_index] >= 0) return genotype_relatives[PARENT_INDEX][start_index];
	else return start_index;
}

/////////////////////// generate preamble gnuplot code for images:

char *preamble(char *text, char *buffer_text_long)
{

	char limit[15];
	
	if(text == NULL) {printf("Memory problem!\n"); exit(0);}

	sprintf(limit, "%d", dim_grid - 1);
	strcpy(text, "");
	strcpy(text, concat(text,(char *)"set size square\n", buffer_text_long)); // sets the canvas to be square
	strcpy(text, concat(text,(char *)"unset key\n", buffer_text_long)); // disables a key (or legend) describing plots on a plot
	strcpy(text, concat(text,(char *)"set tic scale 0\n", buffer_text_long)); // control of the major (labelled) tics on all axes at once
	
	// set the range of values to be coloured using the current palette:
	strcpy(text, concat(text,(char *)"set cbrange [0:5]\n", buffer_text_long));
	// define the divisions in the palette
	strcpy(text, concat(text,(char *)"set palette model RGB defined (", buffer_text_long));
	strcpy(text, concat(text,(char *)"0 \"blue\", 1 \"red\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"1 \"red\", 2 \"yellow\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"2 \"yellow\", 3 \"white\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"3 \"grey\", 4 \"grey\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"4 \"black\", 5 \"black\"", buffer_text_long));
	strcpy(text, concat(text,(char *)")\n", buffer_text_long));
	
	strcpy(text, concat(text,(char *)"unset cbtics\n", buffer_text_long)); // removes major (labelled) tics on the color box axis
	
	// x range of the plot:
	strcpy(text, concat(text,(char *)"set xrange [-0.5:", buffer_text_long));
	strcpy(text, concat(text,limit, buffer_text_long));
	strcpy(text, concat(text,(char *)".5]\n", buffer_text_long));
	
	// y range of the plot:
	strcpy(text, concat(text,(char *)"set yrange [-0.5:", buffer_text_long));
	strcpy(text, concat(text,limit, buffer_text_long));
	strcpy(text, concat(text,(char *)".5]\n", buffer_text_long));
	
	strcpy(text, concat(text,(char *)"set view map\n", buffer_text_long)); // converts surface to a two-dimensional 'map' style view
	strcpy(text, concat(text,(char *)"splot '-' matrix with image\n", buffer_text_long)); // makes the plot (splot plots 3-d surfaces and data)

	return text;
}

char *preamble_drivers(char *text, char *buffer_text_long)
{
	// returns a string containing gnuplot code needed to set up each plot

	char limit[15];
	
	if(text == NULL) {printf("Memory problem!\n"); exit(0);}

	sprintf(limit, "%d", dim_grid - 1);
	strcpy(text, "");
	strcpy(text, concat(text,(char *)"set size square\n", buffer_text_long)); // sets the canvas to be square
	strcpy(text, concat(text,(char *)"unset key\n", buffer_text_long)); // disables a key (or legend) describing plots on a plot
	strcpy(text, concat(text,(char *)"set tic scale 0\n", buffer_text_long)); // control of the major (labelled) tics on all axes at once
	
	// set the range of values to be coloured using the current palette:
	strcpy(text, concat(text,(char *)"set cbrange [0:27]\n", buffer_text_long));
	// define the divisions in the palette
	strcpy(text, concat(text,(char *)"set palette model RGB defined (", buffer_text_long));
	strcpy(text, concat(text,(char *)"0 \"#8A7C64\", 1 \"#8A7C64\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"1 \"#599861\", 2 \"#599861\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"2 \"#89C5DA\", 3 \"#89C5DA\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"3 \"#DA5724\", 4 \"#DA5724\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"4 \"#74D944\", 5 \"#74D944\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"5 \"#CE50CA\", 6 \"#CE50CA\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"6 \"#3F4921\", 7 \"#3F4921\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"7 \"#C0717C\", 8 \"#C0717C\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"8 \"#CBD588\", 9 \"#CBD588\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"9 \"#5F7FC7\", 10 \"#5F7FC7\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"10 \"#673770\", 11 \"#673770\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"11 \"#D3D93E\", 12 \"#D3D93E\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"12 \"#38333E\", 13 \"#38333E\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"13 \"#508578\", 14 \"#508578\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"14 \"#D7C1B1\", 15 \"#D7C1B1\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"15 \"#689030\", 16 \"#689030\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"16 \"#AD6F3B\", 17 \"#AD6F3B\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"17 \"#CD9BCD\", 18 \"#CD9BCD\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"18 \"#D14285\", 19 \"#D14285\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"19 \"#6DDE88\", 20 \"#6DDE88\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"20 \"#652926\", 21 \"#652926\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"21 \"#7FDCC0\", 22 \"#7FDCC0\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"22 \"#C84248\", 23 \"#C84248\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"23 \"#8569D5\", 24 \"#8569D5\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"24 \"#5E738F\", 25 \"#5E738F\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"25 \"#D1A33D\", 26 \"#D1A33D\", ", buffer_text_long));
	strcpy(text, concat(text,(char *)"26 \"#000000\", 27 \"#000000\"", buffer_text_long));
	strcpy(text, concat(text,(char *)")\n", buffer_text_long));

	strcpy(text, concat(text,(char *)"unset cbtics\n", buffer_text_long)); // removes major (labelled) tics on the color box axis
	
	// x range of the plot:
	strcpy(text, concat(text,(char *)"set xrange [-0.5:", buffer_text_long));
	strcpy(text, concat(text,limit, buffer_text_long));
	strcpy(text, concat(text,(char *)".5]\n", buffer_text_long));
	
	// y range of the plot:
	strcpy(text, concat(text,(char *)"set yrange [-0.5:", buffer_text_long));
	strcpy(text, concat(text,limit, buffer_text_long));
	strcpy(text, concat(text,(char *)".5]\n", buffer_text_long));
	
	strcpy(text, concat(text,(char *)"set view map\n", buffer_text_long)); // converts surface to a two-dimensional 'map' style view
	strcpy(text, concat(text,(char *)"splot '-' matrix with image\n", buffer_text_long)); // makes the plot (splot plots 3-d surfaces and data)

	return text;
}

/////////////////////// generic functions:

// weighted random choice of an element in an array of floats:
int weighted_random_floats(float *rates, long *idum, int start_index, int num_elements)
{
	if(num_elements == 1) return start_index;

	int i;

	buff_small_float[0] = rates[start_index];
	for(i = 1; i < num_elements; i++) buff_small_float[i] = buff_small_float[i - 1] + rates[start_index + i];

	return start_index + weighted_random_known_sums_floats(buff_small_float, idum, num_elements);
}

// weighted random choice of an element in an array of doubles:
int weighted_random_doubles(double *rates, long *idum, int start_index, int num_elements)
{
	if(num_elements == 1) return start_index;

	int i;

	buff_small_double[0] = rates[start_index];
	for(i = 1; i < num_elements; i++) buff_small_double[i] = buff_small_double[i - 1] + rates[start_index + i];

	return start_index + weighted_random_known_sums_doubles(buff_small_double, idum, num_elements);
}

// weighted random choice of an element in an array of integers:
int weighted_random_ints(int *rates, long *idum, int start_index, int num_elements)
{
	if(num_elements == 1) return start_index;

	int i;

	buff_array_int[0] = rates[start_index];
	for(i = 1; i < num_elements; i++) buff_array_int[i] = buff_array_int[i - 1] + rates[start_index + i];

	return start_index + weighted_random_known_sums_ints(buff_array_int, idum, num_elements);
}

int weighted_random_known_sums_floats(float *cumulative_rates, long *idum, int num_elements)
{
	if(num_elements == 1) return 0;

	double r = ran1(idum) * cumulative_rates[num_elements - 1];

	if(num_elements == 2) {
		if(cumulative_rates[0] < r) return 1;
		else return 0;
	}

	int lowGuess = 0;
	int highGuess = num_elements - 1;
	int guess;

	while(highGuess >= lowGuess) {
		guess = (lowGuess + highGuess) / 2;
		if(cumulative_rates[guess] < r) lowGuess = guess + 1;
		else if(guess > 0 && cumulative_rates[guess - 1] > r) highGuess = guess - 1;
		else return guess;
	}
	return guess;
}

int weighted_random_known_sums_doubles(double *cumulative_rates, long *idum, int num_elements)
{
	if(num_elements == 1) return 0;

	double r = ran1(idum) * cumulative_rates[num_elements - 1];

	if(num_elements == 2) {
		if(cumulative_rates[0] < r) return 1;
		else return 0;
	}

	int lowGuess = 0;
	int highGuess = num_elements - 1;
	int guess;

	while(highGuess >= lowGuess) {
		guess = (lowGuess + highGuess) / 2;
		if(cumulative_rates[guess] < r) lowGuess = guess + 1;
		else if(guess > 0 && cumulative_rates[guess - 1] > r) highGuess = guess - 1;
		else return guess;
	}
	return guess;
}

int weighted_random_known_sums_ints(int *cumulative_rates, long *idum, int num_elements)
{
	if(num_elements == 1) return 0;

	double r = ran1(idum) * cumulative_rates[num_elements - 1];

	if(num_elements == 2) {
		if(cumulative_rates[0] < r) return 1;
		else return 0;
	}

	int lowGuess = 0;
	int highGuess = num_elements - 1;
	int guess;

	while(highGuess >= lowGuess) {
		guess = (lowGuess + highGuess) / 2;
		if(cumulative_rates[guess] < r) lowGuess = guess + 1;
		else if(guess > 0 && cumulative_rates[guess - 1] > r) highGuess = guess - 1;
		else return guess;
	}
	return guess;
}

// create a string of zeroes to prefix a number in a file name:
char *zeroes(int gen, int maxgen, char *buffer_text_short, char *buffer_text_long)
{
	int i;
	
	if(buffer_text_short == NULL) {printf("Memory problem!\n"); exit(0);}
	
	strcpy(buffer_text_short, "");
	for(i = (int)MAX((int)(log(gen + 0.1)/log(10)), 0); i < (int)(log(maxgen + 0.1) / log(10)); i++) buffer_text_short = concat(buffer_text_short, (char *)"0", buffer_text_long);
	// the MAX function is invoked to avoid attempted calculation of log(0), which the compiler equates to some negative number of large magnitude;
	// 0.1 is added to avoid rounding errors (e.g. to avoid log10(1000) = 2.99999)
	
	return buffer_text_short;
}

// return the result of concatenating two strings:
char *concat(char *s1, char *s2, char *buffer_text_long)
{
	strcpy(buffer_text_long, s1);
	strcat(buffer_text_long, s2);
	
	return buffer_text_long;
}

// a random number generator:
double ran1(long *idum)
{
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	double temp;

	if(*idum <= 0 || !iy) {
		if(-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		for(j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if(*idum < 0) *idum += IM;
			if(j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if(*idum < 0) *idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;

	if((temp = AM * iy) > RNMX) return RNMX;
	else return temp;
}

int which_quadrant(int x, int y, float theta, float tan_theta, int l)
{
	if(theta < PI / 2) {
		if(tan_theta < (float)x / y) return 4;
		else return 1;
	}
	else if(theta < PI) {
		if(tan_theta > (float)x / (y - l)) return 2;
		else return 1;
	}
	else if(theta < 3 * PI / 2) {
		if(tan_theta < (float)(x - l) / (y - l)) return 2;
		else return 3;
	}
	else {
		if(tan_theta > (float)(x - l) / y) return 4;
		else return 3;
	}
}

// return an exponentially distributed, positive, random deviate of unit mean:
float expdev(long *idum)
{ 
	float dum;

	do {dum = ran1(idum);} while (dum == 0.0);

	return -log(dum);
}

// sample from a hypergeometric distribution (subpopulation sizes n1 and n2; t samples)
unsigned int hypergeometric(long *idum, unsigned int n1, unsigned int n2, unsigned int t)
{
	const unsigned int n = n1 + n2; // total population size
	unsigned int i;
	unsigned int a = n1;
	unsigned int b = n1 + n2;
	unsigned int k = 0; // number of samples found so far

	if (t > n) t = n; // sample size can't exceed population size

	if (t < n/2) {
		for (i = 0; i < t; i++) {
			double u = ran1(idum);
			if (b * u < a) {
				k++;
				if (k == n1) return k;
				a--;
			}
			b-- ;
		}
		return k;
	}
	else {
		for (i = 0; i < n - t; i++) {
			double u = ran1(idum);
			if (b * u < a) {
				k++;
				if (k == n1) return n1 - k;
				a--;
			}
			b--;
		}
		return n1 - k;
	}
}

// sample from a Poisson distribution
unsigned int poisson(long *idum, double mu)
{
	double emu;
	double prod = 1.0;
	unsigned int k = 0;

	while(mu > 10) { // if the mutation rate is extremely high
		unsigned int m = mu * (7.0 / 8.0);
		double X = gamma(idum, m);

		if (X >= mu) return k + binomial(idum, mu / X, m - 1);
		else {
			k += m;
			mu -= X;
		}
	}

	emu = exp(-mu);

	do {
		prod *= ran1(idum);
		k++;
	} while (prod > emu);

	return k - 1;
}

// sample from a binomial distribution with np < 14
unsigned int binomial(long *idum, double p, unsigned int n)
{
	int ix;
	int flipped = 0;
	double q, s;

	if(n * p >= 14) {
		printf("Warning: in binomial function, np = %f >= 14\n", n * p);
		fprintf(error_log, "Warning: in binomial function, np = %f >= 14\n", n * p);
	}

	if(n == 0) return 0;

	if(p > 0.5) {
		p = 1.0 - p;
		flipped = 1;
	}

	q = 1 - p;
	s = p / q;

	double f0 = power(q, n);

	while(1) {
		double f = f0;
		double u = ran1(idum);

		for (ix = 0; ix <= 110; ++ix) {
			if (u < f) break;
			u -= f;
			f *= s * (n - ix) / (ix + 1);
		}
	}

	return (flipped) ? (n - ix) : (unsigned int)ix;
}

// sample from a gamma distribution with a < 12
double gamma(long *idum, const unsigned int a)
{
	unsigned int i;
	double prod = 1;
	double rnd;

	if(a >= 12) {
		printf("Warning: in gamma function, a = %d >= 12\n", a);
		fprintf(error_log, "Warning: in gamma function, a = %d >= 12\n", a);
	}

	for(i = 0; i < a; i++) {
		do {rnd = ran1(idum);} while(rnd == 0);
		prod *= rnd;
	}

	return -log(prod);
}

// x^n for small integer n
double power(double x, unsigned int n)
{
	double value = 1.0;
	
	do {
		if(n & 1) value *= x;
		n >>= 1;
		x *= x;
	} while(n);

	return value;
}

// returns 1 if the binary representation of "sum" has a 1 at the position specified by "part"
// e.g. component(7, x) = 1 for x = 1,2,4 because 7 = 1+2+4
int component(int sum, int part)
{
    int targetlevel = 0;
    while (part >>= 1) ++targetlevel;

    if((sum>>targetlevel) & 1) return 1;
    else return 0;
}

// stochastically round a positive float either up or down:
int stochastic_round(float a, long *idum)
{ 
	float rnd = ran1(idum);
	if(rnd < a - int(a)) return(int(a) + 1);
	else return(int(a));
}

// malloc an array of integers:
void mallocArray_int(int ***a, int m1, int m2) 
{
	int i;
	*a = (int **)malloc(m1 * sizeof(int *));
	if(*a == NULL) {printf("Memory problem in mallocArray_int (m1 = %d)\n", m1); exit(1);}
	for (i=0; i<m1; i++) {
		(*a)[i] = (int *)malloc(m2 * sizeof(int));
		if((*a)[i] == NULL) {printf("Memory problem in mallocArray_int (m2 = %d)\n", m2); exit(1);}
	}
}

// malloc an array of floats:
void mallocArray_float(float ***a, int m1, int m2) 
{
	int i;
	*a = (float **)malloc(m1 * sizeof(float *));
	if(*a == NULL) {printf("Memory problem in mallocArray_float (m1 = %d)\n", m1); exit(1);}
	for (i=0; i<m1; i++) {
		(*a)[i] = (float *)malloc(m2 * sizeof(float));
		if((*a)[i] == NULL) {printf("Memory problem in mallocArray_float (m2 = %d)\n", m2); exit(1);}
	}
}

// malloc an array of doubles:
void mallocArray_double(double ***a, int m1, int m2) 
{
	int i;
	*a = (double **)malloc(m1 * sizeof(double *));
	if(*a == NULL) {printf("Memory problem in mallocArray_double (m1 = %d)\n", m1); exit(1);}
	for (i=0; i<m1; i++) {
		(*a)[i] = (double *)malloc(m2 * sizeof(double));
		if((*a)[i] == NULL) {printf("Memory problem in mallocArray_double (m2 = %d)\n", m2); exit(1);}
	}
}

// free an array of integers:
void freeArray_int(int **a, int m) 
{
	int i;
	for (i = 0; i < m; ++i) free(a[i]);
	free(a);
}

// free an array of floats:
void freeArray_float(float **a, int m) 
{
	int i;
	for (i = 0; i < m; ++i) free(a[i]);
	free(a);
}

// free an array of doubles:
void freeArray_double(double **a, int m) 
{
	int i;
	for (i = 0; i < m; ++i) free(a[i]);
	free(a);
}