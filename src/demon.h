#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <iostream>
#include <string>
#include <cfloat>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ptree.hpp>

///// mathematical constants:
#define INVROOT2 0.70710678118
#define ROOT2 1.41421356237
#define PI 3.14159265359
#define TRUE 1
#define FALSE 0

// macros:
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define SIGN(X) (0 < X) - (X < 0) // returns -1 if X<0, 0 if X=0, 1 if X>0
#define ROUNDIT(X) floor(X + 0.5)

///// maximum grid width for gnuplot:
#define max_grid_for_output 523

///// options for main_calculations_and_output function:
#define PHYLO_AND_POPS 1
#define DIVERSITIES 2
#define GENOPROPS 4

///// event types (for event_type variable and event_counter array):
#define BIRTH_EVENT 0
#define DEATH_EVENT 1
#define MIGRATION_EVENT 2
#define MUTATION_EVENT 3
#define FISSION_EVENT 4
#define NORMALBIRTH_EVENT 5
#define NORMALDEATH_EVENT 6

///// model-specific constants:
#define EMPTY -9 // value of grid[][] when the site has not yet been occupied by a cancer cell
#define SET_SIZE 3 // max number of elements per bintree (in layer 0, each bintree element contains <=3 demes; in subsequent layers, each bintree element contains <=3 elements from the previous layer)
#define MAX_DRIVERS_TO_COUNT 10 // max value of n when counting number of cells with n drivers
#define MAX_TRIALS 200 // number of times to restart simulation in case of entinction, before giving up;
// Note: probability of success in a single trial, starting from a single cell, is p_fix = (1-r)/(1-r^K) > 1-r, where r = 1/normal_birth_rate; probability of failure in MAX_TRIALS is 1-(1-p_fix)^MAX_TRIALS
// For example, if normal_birth_rate = 0.95 and MAX_TRIALS = 200 then the chance of a successful trial is more than 99.996%

///// clone:
#define POPULATION 0 // number of cancer cells
#define DEME 1 // index of the deme in which the clone is located
#define GENOTYPE 2 // index of the genotype of the clone
#define DRIVER_GENOTYPE 3 // index of the driver genotype of the clone
#define INDEX_IN_DEME 4 // index of the clone within its deme

#define NUM_CLONE_INT_PROPS 5

///// genotype and driver genotype:
//POPULATION already defined as 0
#define PARENT 1 // parent's unique ID
#define IDENTITY 2 // unique ID
#define DRIVER_IDENTITY 3 // unique ID of corresponding driver genotype
#define NUM_DRIVER_MUTATIONS 4 // number of driver mutations
#define NUM_MIGRATION_MUTATIONS 5 // number of migration rate mutations
#define IMMORTAL 6 // whether genotype record can be overwritten
#define NUM_PASSENGER_MUTATIONS 7 // number of passenger mutations

#define NUM_GENOTYPE_INT_PROPS 8
//
#define BIRTH_RATE 0 // birth rate conferred by the genotype
#define MIGRATION_RATE 1 // migration rate conferred by the genotype
#define ORIGIN_TIME 2 // generation at which genotype originated

#define NUM_GENOTYPE_FLOAT_PROPS 3

///// genotype relatives:
#define PARENT_INDEX 0
#define FIRST_DAUGHTER 1
#define LAST_DAUGHTER 2
#define NEXT_SISTER 3

///// deme:
//POPULATION already defined as 0
#define XCOORD 1 // x coordinate of the deme
#define YCOORD 2 // y coordinate of the deme
#define NORMAL_CELLS 3 // number of normal cells in the deme
#define NUM_CLONES_IN_DEME 4 // number of clones in the deme

#define NUM_DEME_INT_PROPS 5
//
#define DEATH_RATE 0 // death rate of the deme
#define MIGRATION_MODIFIER 1 // factor by which migration rate is multiplied (depends on what cells are in the deme)

#define NUM_DEME_FLOAT_PROPS 2
//
#define SUM_BIRTH_RATES 0 // sum of birth rates of cancer cells in the deme
#define SUM_MIGRATION_RATES 1 // sum of migration rates of cancer cells in the deme
#define SUM_RATES 2 // sum of all rates of cancer cells and normal cells in the deme

#define NUM_DEME_DOUBLE_PROPS 3

///// parameters of random number generator:
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

//////////////

// set up:
char* get_input_path(int argc, char *argv[]);
void read_parameters(boost::property_tree::ptree pt);
float set_init_migration_rate(int K, float init_migration_rate, float A, float B, float C);
void initialise(int *num_cells, int *num_clones, int *num_demes, int *num_matrix_cols, int *num_empty_cols, int init_driver_birth_mutations, 
	int init_driver_mig_mutations, int init_passengers, int *num_empty_driver_cols, int *num_driver_matrix_cols, int *next_driver_genotype_id, 
	int *next_genotype_id, long *idum, int *num_extinct_genotypes, int *num_empty_demes, int *num_extinct_driver_genotypes);

// run the simulation:
void run_sim(char *input_and_output_path, char *config_file_with_path);

// choose update type:
int choose_bintree_doubles(double *bintree, int max_layer_needed, int max_array_index, int array_size, long *idum);
int choose_bintree_ints(int *bintree, int max_layer_needed, int max_array_index, int array_size, long *idum);
int choose_clone_without_bintree(int chosen_deme, float *buff_array, float deme_sum_cancer, long *idum);
int choose_clone_with_bintree_by_population(int bintree_index, int chosen_deme, int num_clones_in_deme, long *idum);
int choose_clone_with_bintree_by_rates(int bintree_index, int chosen_deme, int num_clones_in_deme, long *idum);
int choose_event_for_clone(bool include_death, int chosen_deme, int chosen_clone, float *buff_array, long *idum);
int choose_event_for_deme(int chosen_deme, float *buff_array, long *idum);

// cell events (top level):
void cell_division(int *event_counter, int *num_cells, int parent_deme_num, int *new_passengers, int *new_mig_mutations, int *new_birth_mutations, 
	int *new_mutations, long *idum, int chosen_clone, int parent_geno_num, int *daughter_clone_nums, int *num_empty_cols, int *num_matrix_cols, 
	int *empty_cols, int *num_clones, int parent_driver_geno_num, int *num_empty_driver_cols, int *num_driver_matrix_cols, int *empty_driver_cols, 
	int *next_driver_genotype_id, int num_demes, int *next_genotype_id, int *num_extinct_genotypes, int *num_empty_demes, 
	int *num_extinct_driver_genotypes, float gens_elapsed);
void cell_death(int *event_counter, int *num_cells, int parent_deme_num, int parent_geno_num, int *empty_cols, int *num_empty_cols,
	int chosen_clone, int *num_clones, int parent_driver_geno_num, int *num_empty_driver_cols, int *empty_driver_cols, int num_demes, 
	int *num_extinct_genotypes, int *num_empty_demes, int *num_extinct_driver_genotypes);
void cell_migration(int *event_counter, int parent_deme_num, long *idum, int *num_demes, int *num_clones, int parent_clone, int *num_cells,
	int daughter_geno_num, int daughter_driver_geno_num, int *num_empty_demes, int *empty_cols, int *num_empty_cols, int *num_empty_driver_cols, int *empty_driver_cols, 
	int *num_extinct_genotypes, int *num_extinct_driver_genotypes);
void deme_fission(int *event_counter, int origin_deme_num, long *idum, int *num_demes, int *num_clones, int *num_cells, int *num_empty_demes, int *num_empty_cols, int *num_empty_driver_cols, 
	int *empty_cols, int *empty_driver_cols, int *num_extinct_genotypes, int *num_extinct_driver_genotypes, int num_matrix_cols);

// genotype and driver genotype events (lower level):
void choose_number_mutations(int *new_passengers, int *new_mig_mutations, int *new_birth_mutations, int *new_mutations, long *idum);
int select_genotype_index(int *num_empty_cols, int *num_matrix_cols, int *empty_cols);
void create_genotype(int **geno_or_driver_ints, float **geno_or_driver_floats, int *num_matrix_cols, int daughter_geno_num, int parent_geno_num, int *next_genotype_id, int daughter_driver_id, 
	float new_birth_rate, float new_migration_rate, int new_passengers, int new_birth_mutations, int new_mig_mutations, float gens_elapsed);
void increment_or_decrement_genotype(int **geno_or_driver_ints, int parent_geno_num, int *empty_cols, int *num_empty_cols, int change, int *num_extinct_genotypes);
void create_column(int **either_matrix, int num_matrix_cols, int parent_geno_num, int daughter_geno_num, int num_mutations);

// clone events (lower level):
int select_clone_in_deme(int deme_index, int genotype_index);
void create_clone(int daughter_geno_num, int *num_clones, int daughter_deme_num, int daughter_driver_geno_num, int num_demes, int initial_clone_pop);
void increment_or_decrement_clone(int chosen_clone, int deme_index, int *num_clones, int change, int num_demes);
void remove_clone(int chosen_clone, int deme_index, int *num_clones, int num_clones_in_deme);

// deme events (lower level):
void move_cells(long *idum, int origin_deme_num, int dividing_beyond_the_edge, int new_deme_index, int *num_cells, int *num_demes, int *num_empty_demes, int *num_clones, 
	int *num_empty_cols, int *num_empty_driver_cols, int *num_extinct_genotypes, int *num_extinct_driver_genotypes, int *event_counter);
void budge_demes(int old_x, int old_y, int *x_to_fill, int *y_to_fill);
void get_deme_coordinates(int *x_to_fill, int *y_to_fill, int old_x, int old_y, long *idum);
void choose_grid_square(int old_deme_index, long *idum, int *new_x, int *new_y);
void create_deme(int new_x, int new_y, int *num_demes, int num_cells);
void remove_deme(int deme_index, int *num_cells, int *num_clones, int *num_demes, int *num_empty_demes, int *empty_cols, int *num_empty_cols, int *num_empty_driver_cols, int *empty_driver_cols, 
	int *num_extinct_genotypes, int *num_extinct_driver_genotypes);
void increment_or_decrement_deme(int change, int deme_index, int num_cells, int num_demes, int *num_empty_demes);
void add_or_remove_normal_cell(int change, int deme_index, int num_cells, int *event_counter, int num_demes);

// memory handling:
void assign_memory();
void free_memory();

// file handling (top level):
void open_files(char *input_and_output_path);
void initiate_files(int *num_samples_list);
void main_calculations_and_output(long *idum, int num_demes, int num_matrix_cols, int num_driver_matrix_cols, float gens_elapsed, int *driver_counts, int num_cells, int num_clones, 
	int record_phylogenies, int *num_samples_list, int next_genotype_id, int next_driver_genotype_id, int *event_counter, int num_empty_cols, int num_empty_driver_cols, int num_extinct_genotypes, 
	int num_extinct_driver_genotypes, int num_empty_demes, long t1, int print_output, int which_parts);
void grids_output(char *preamble_text, char *preamble_drivers_text, char *preamble_passengers_text, float gens_elapsed, int num_clones, int num_demes, char *input_and_output_path, 
	char *buffer_text_short, char *buffer_text_long, bool to_file);
void end_of_loop_output(int num_cells, float gens_elapsed, long t1);
void final_output(int trial_num, int num_cells, long t1);
void close_files();

// write system state to screen or file (lower level):
void print_to_screen(float gens_elapsed, int num_cells, long t1, int num_demes, int num_matrix_cols, int num_clones, float mean_num_passengers, float mean_num_drivers, float diversity, float driver_diversity, 
	double sum_birth_rates, double sum_death_rates, double sum_migration_rates, int num_empty_cols, float alpha_diversity, float alpha_driver_diversity, int *event_counter, int num_driver_matrix_cols, 
	double sum_normal_death_rates, double sum_normal_birth_rates, float edge_diversity, int num_empty_driver_cols, int num_empty_demes, int num_extinct_genotypes, int num_extinct_driver_genotypes, int next_genotype_id);
void write_output_pops(FILE *output_pops, int num_cells, int *event_counter, float mean_num_passengers, float mean_num_drivers, double sum_birth_rates, double sum_migration_rates, 
	float var_num_passengers, float var_num_drivers, float variance_birth_rate, float variance_mig_rate, int num_matrix_cols, int num_empty_cols, int num_driver_matrix_cols, int num_empty_driver_cols, 
	int num_extinct_genotypes, int num_extinct_driver_genotypes, int num_demes, int num_empty_demes, float gens_elapsed, int num_clones, int *driver_counts, int next_genotype_id);
void write_diversities(FILE *output_diversities, float diversity, float alpha_diversity, float edge_diversity, float driver_diversity, float alpha_driver_diversity, float edge_driver_diversity, 
	float **depth_diversity, float **depth_diversity_bigsample, float gens_elapsed, int num_cells);
void write_other_files(FILE *output_demes, FILE *output_clones, FILE *output_genotype_counts, FILE *output_driver_genotype_counts, FILE *output_phylo, float gens_elapsed, 
	int num_demes, int num_matrix_cols, int num_driver_matrix_cols, int num_clones, float *within_deme_diversity, float *within_deme_driver_diversity, int num_cells);
void write_output_phylo(FILE* output, int num_cols, float gens_elapsed, int *populations, int **genotype_or_driver_ints, float **genotype_or_driver_floats, 
	int samples, int biopsy_size_per_sample, int depth, int num_cells);
void write_frequency_table(FILE *output_allele_counts, int **freq_table, int num_cells, float gens_elapsed, int num_freqs);
void write_genotypes(FILE *output_genotype_properties, int num_matrix_cols, int *allele_count, int **genotype_or_driver_ints, float **genotype_or_driver_floats);
void write_matrix_to_file(FILE *output_matrix, int **either_matrix, int num_matrix_cols);
void write_pops_grid_to_file(FILE *output_popgrid);
void write_passengers_grid_to_file(FILE *output_passengersgrid, int num_clones, int num_demes);
void write_normalcells_grid_to_file(FILE *output_normalcellsgrid);
void write_deathrates_grid_to_file(FILE *output_deathrates_grid);
void write_rates_grid_to_file(FILE *output, double *deme_rates, int *deme_pops);
void write_drivers_grid_to_file(FILE *output);
void plot_birth_rate_grid(FILE *gp, char *preamble_text, float gens_elapsed, int num_clones, char *input_and_output_path, char *buffer_text_short, char *buffer_text_long);
void plot_drivers_grid(FILE *gp, char *preamble_drivers_text, float gens_elapsed, char *input_and_output_path, char *buffer_text_short, char *buffer_text_long);
void plot_pops_grid(FILE *gp, char *preamble_text, float gens_elapsed, char *input_and_output_path, char *buffer_text_short, char * buffer_text_long);
void plot_passengers_grid(FILE *gp, char *preamble_text, float gens_elapsed, int num_clones, int num_demes, char *input_and_output_path, char *buffer_text_short, char *buffer_text_long);
void plot_migration_grid(FILE *gp, char *preamble_text, float gens_elapsed, int num_clones, char *input_and_output_path, char *buffer_text_short, char *buffer_text_long);

// set rates:
float set_birth_rate(int new_birth_mutations, int new_passengers, float parent_birth_rate, long *idum);
float set_death_rate(int deme_index, int num_cells);
float set_migration_rate(int new_mig_mutations, float parent_mig_rate, long *idum);
float set_migration_modifier(int normal_cells, int population);

// calculate or reset sums of rates:
void update_deme_bintree_environmental_rates(int change_pop, int change_normal_cells, int num_demes, int deme_index, int num_cells);
void update_deme_bintree_phenotype_rates(int chosen_clone, int deme_index, int change, int num_demes);
void set_clone_in_deme(int deme_num, int geno_num, int index_in_deme, int clone_pop);
void reset_deme_and_bintree_sums(int num_demes, int num_clones, int print_always);
void reset_clone_bintree_sums(int num_demes, int num_clones);
void calculate_sums_of_rates(double *sum_death_rates, double *sum_birth_rates, double *sum_migration_rates, double *sum_normal_birth_rates, double *sum_normal_death_rates, int num_demes, 
	int num_matrix_cols, float gens_elapsed);
float sum_of_all_rates(int num_demes);

// checks:
void check_clone_bintree_sums(bool full);
void check_deme_sums(int chosen_deme, float deme_sum_cancer);
void check_rates_sum(float b0, float b1, int chosen_deme, int cell_type);
void check_geno_populations(int chosen_clone, int event_type);
void check_matrix_cols1(int num_matrix_cols, int num_driver_matrix_cols, int num_cells);
void check_matrix_cols2(int num_matrix_cols, int num_driver_matrix_cols, int num_empty_cols, int num_empty_driver_cols, int num_extinct_genotypes, int num_extinct_driver_genotypes);
void check_clones_in_deme1(int origin_deme_num, int clone_num);
void check_clones_in_deme2(int num_clones_in_deme);
void check_chosen_clone(int chosen_clone, int num_clones_in_deme);
void check_clone_populations(int chosen_clone, int event_type, int parent_deme_num);
void check_normal_pops(int chosen_deme);
void check_genotype_counts(int num_matrix_cols, int num_empty_cols, int num_extinct_genotypes, int num_driver_matrix_cols, int num_empty_driver_cols, int num_extinct_driver_genotypes, 
	int cell_type, int event_type, int chosen_deme);

// generic bintree functions:
void set_bintree_sums_layer0(double *bintree, double *array, int num_array_elements);
void set_bintree_sums_layer0_int(int *bintree, int *array, int num_array_elements);
void set_bintree_sums_subsequent_layers(double *bintree, int num_demes, int max_layer_needed, int array_length);
void set_bintree_sums_subsequent_layers_int(int *bintree, int num_demes, int max_layer_needed, int array_length);
void update_bintree_layers(double summand, double *bintree, int array_index, int max_layer_needed, int array_length);
void update_bintree_layers_int(int summand, int *bintree, int array_index, int max_layer_needed, int array_length);
void extend_bintree(double *bintree, int array_index, int max_layer_needed, int array_length);
void extend_bintree_int(int *bintree, int array_index, int max_layer_needed, int array_length);
int get_bintree_index(int layer, int array_index, int array_size);
int get_deme_bintree_index(int layer, int array_index);
int get_clone_bintree_index(int layer, int array_index);
int get_max_layer_needed(int n);

// calculate diversity metrics (top level):
void get_diversity_metrics(int rao, float *diversity, float *edge_diversity, float *alpha_diversity, float *within_deme_diversity, int num_demes, int num_cells, int num_matrix_cols, int num_clones, 
	int *clone_genotype, int **either_matrix, int dmax, int *populations, long *idum, float gens_elapsed, int* position_in_edge_list, int* demes_at_edge, int* sides_at_edge, int* genotype_edge_pop, int is_drivers_only);
void get_biopsy_data(int rao, float *depth_array, int samples, int biopsy_size_per_sample, int num_demes, int num_cells, int num_matrix_cols, int num_clones, int *clone_genotype, int **either_matrix, 
	int dmax, long *idum, FILE *sample_size_log, float gens_elapsed, FILE *output_phylo_of_sample, int calculate_sample_diversities, int record_phylogenies, float centre_X, float centre_Y, 
	int **genotype_or_driver_ints, float **genotype_or_driver_floats, int *at_edge);
void calculate_mutation_metrics(float *mean_num_passengers, float *mean_num_drivers, float *var_num_passengers, float *var_num_drivers, int *driver_counts, int num_matrix_cols, int num_cells);
float calculate_variance_of_rate(int num_matrix_cols, int num_cells, double sum_of_rates, float *list_of_rates);
void get_relatives(int num_matrix_cols, int next_genotype_id, int **geno_or_driver_ints);
void get_allele_frequencies(int num_matrix_cols, int *allele_count, int **geno_or_driver_ints);
void get_frequency_table(int input_length, int *count, int **freq_table, int *output_length);

// calculate diversity metrics (lower level):
float centre_of_gravity(int direction, int num_cells);
float calculate_diversity(int num_matrix_cols, int N, int *populations, int **either_matrix, int dmax, int rao);
float calculate_edge_diversity(int num_matrix_cols, int num_cells, int **either_matrix, int dmax, int num_demes, int num_clones, int *clone_genotype, long *idum, int rao, int* position_in_edge_list, 
	int* demes_at_edge, int* sides_at_edge, int* genotype_edge_pop);
float calculate_within_deme_diversity(int deme_index, int num_matrix_cols, int num_clones, int *clone_genotype, int **either_matrix, int dmax, int rao, int is_drivers_only, int* genotype_edge_pop);
void find_sample_genotype_pops(int *genotype_populations_in_sample, int *sampled_cells, int num_matrix_cols, int num_cells, int **either_matrix, int dmax, int num_demes, int num_clones, 
	int depth, float centre_X, float centre_Y, int *clone_genotype, int num_directions, int biopsy_size_per_sample, long *idum, FILE *sample_size_log, float gens_elapsed, int rao, int *at_edge);
int move_down(int start_index);
int move_right(int start_index);
int move_up(int start_index);

// generate preamble gnuplot code for images:
char *preamble(char *text, char *buffer_text_long);
char *preamble_drivers(char *text, char *buffer_text_long);

// generic functions:
int weighted_random_floats(float *rates, long *idum, int start_index, int num_elements);
int weighted_random_doubles(double *rates, long *idum, int start_index, int num_elements);
int weighted_random_ints(int *rates, long *idum, int start_index, int num_elements);
int weighted_random_known_sums_floats(float *cumulative_rates, long *idum, int num_elements);
int weighted_random_known_sums_doubles(double *cumulative_rates, long *idum, int num_elements);
int weighted_random_known_sums_ints(int *cumulative_rates, long *idum, int num_elements);
char *zeroes(int, int, char *buffer_text_short, char *buffer_text_long);
char *concat(char *, char *, char *buffer_text_long);
double ran1(long *);
int which_quadrant(int x, int y, float theta, float tan_theta, int l);
float expdev(long *);
unsigned int hypergeometric(long *, unsigned int, unsigned int, unsigned int);
unsigned int poisson(long *idum, double mu);
unsigned int binomial(long *idum, double p, unsigned int n);
double gamma(long *idum, const unsigned int a);
double power(double x, unsigned int n);
int component(int sum, int part);
int stochastic_round(float, long *);
void mallocArray_int(int ***a, int m1, int m2);
void mallocArray_float(float ***a, int m1, int m2);
void mallocArray_double(double ***a, int m1, int m2);
void freeArray_int(int **a, int m);
void freeArray_float(float **a, int m);
void freeArray_double(double **a, int m);