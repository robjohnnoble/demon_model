#include "configuration_structures.h"

#include <iostream>
#include <string>

//initializes parameters
void Parameters::GetFromFile(std::string & file_name) {

	  boost::property_tree::ptree pt;
    boost::property_tree::info_parser::read_info(file_name, pt);

    spatial_structure.log2_deme_carrying_capacity =
			pt.get <int> ("parameters.spatial_structure.log2_deme_carrying_capacity");

    dispersal.migration_type =
			pt.get <int> ("parameters.dispersal.migration_type");
    dispersal.init_migration_rate =
			pt.get <float> ("parameters.dispersal.init_migration_rate");
    dispersal.migration_edge_only  =
			pt.get <int> ("parameters.dispersal.migration_edge_only");
    dispersal.migration_rate_scales_with_K  =
			pt.get <int> ("parameters.dispersal.migration_rate_scales_with_K");

		mutation_rates.mu_driver_birth =
			pt.get <float> ("parameters.mutation_rates.mu_driver_birth");
		mutation_rates.mu_passenger =
			pt.get <float> ("parameters.mutation_rates.mu_passenger");
		mutation_rates.mu_driver_migration =
			pt.get <float> ("parameters.mutation_rates.mu_driver_migration");
		mutation_rates.passenger_pop_threshold =
			pt.get <int> ("parameters.mutation_rates.passenger_pop_threshold");

		fitness_effects.normal_birth_rate =
			pt.get <float> ("parameters.fitness_effects.normal_birth_rate");
        fitness_effects.baseline_death_rate =
			pt.get <float> ("parameters.fitness_effects.baseline_death_rate");
		fitness_effects.s_driver_birth =
			pt.get <float> ("parameters.fitness_effects.s_driver_birth");
		fitness_effects.s_passenger =
			pt.get <float> ("parameters.fitness_effects.s_passenger");
		fitness_effects.s_driver_migration  =
			pt.get <float> ("parameters.fitness_effects.s_driver_migration");
		fitness_effects.max_relative_birth_rate  =
			pt.get <float> ("parameters.fitness_effects.max_relative_birth_rate");
		fitness_effects.max_relative_migration_rate  =
			pt.get <float> ("parameters.fitness_effects.max_relative_migration_rate");

		non_biological_parameters.init_pop =
			pt.get <int> ("parameters.non_biological_parameters.init_pop");
		non_biological_parameters.matrix_max =
			pt.get <int> ("parameters.non_biological_parameters.matrix_max");
		non_biological_parameters.max_pop =
			pt.get <int> ("parameters.non_biological_parameters.max_pop");
		non_biological_parameters.max_time =
			pt.get <int> ("parameters.non_biological_parameters.max_time");
		non_biological_parameters.max_generations =
			pt.get <int> ("parameters.non_biological_parameters.max_generations");
		non_biological_parameters.seed =
			pt.get <int> ("parameters.non_biological_parameters.seed");
		non_biological_parameters.write_grid =
			pt.get <int> ("parameters.non_biological_parameters.write_grid");
		non_biological_parameters.write_clones_file =
			pt.get <int> ("parameters.non_biological_parameters.write_clones_file");
		non_biological_parameters.write_demes_file =
			pt.get <int> ("parameters.non_biological_parameters.write_demes_file");
		non_biological_parameters.record_matrix =
			pt.get <int> ("parameters.non_biological_parameters.record_matrix");
		non_biological_parameters.write_phylo =
			pt.get <int> ("parameters.non_biological_parameters.write_phylo");
		non_biological_parameters.calculate_total_diversity =
			pt.get <int> ("parameters.non_biological_parameters.calculate_total_diversity");
		non_biological_parameters.biopsy_size_per_sample =
			pt.get <int> ("parameters.non_biological_parameters.biopsy_size_per_sample");

		p_tree = pt;
};

void Parameters::SetToFile(std::string & file_name) {
  boost::property_tree::ptree pt;

  pt.add("parameters.spatial_structure.log2_deme_carrying_capacity" ,
					spatial_structure.log2_deme_carrying_capacity);

  pt.add("parameters.dispersal.migration_type",
					dispersal.migration_type);
  pt.add("parameters.dispersal.init_migration_rate",
					dispersal.init_migration_rate);
  pt.add("parameters.dispersal.migration_edge_only",
					dispersal.migration_edge_only);
  pt.add("parameters.dispersal.migration_rate_scales_with_K",
					dispersal.migration_rate_scales_with_K);

  pt.add("parameters.mutation_rates.mu_driver_birth",
					mutation_rates.mu_driver_birth);
  pt.add("parameters.mutation_rates.mu_passenger",
				  mutation_rates.mu_passenger);
  pt.add("parameters.mutation_rates.mu_driver_migration",
				  mutation_rates.mu_driver_migration);
  pt.add("parameters.mutation_rates.passenger_pop_threshold",
				  mutation_rates.passenger_pop_threshold);

  pt.add("parameters.fitness_effects.normal_birth_rate",
				  fitness_effects.normal_birth_rate);
  pt.add("parameters.fitness_effects.baseline_death_rate",
				  fitness_effects.baseline_death_rate);
  pt.add("parameters.fitness_effects.s_driver_birth",
					fitness_effects.s_driver_birth);
  pt.add("parameters.fitness_effects.s_passenger",
					fitness_effects.s_passenger);
  pt.add("parameters.fitness_effects.s_driver_migration",
					fitness_effects.s_driver_migration);
  pt.add("parameters.fitness_effects.max_relative_birth_rate",
					fitness_effects.max_relative_birth_rate);
  pt.add("parameters.fitness_effects.max_relative_migration_rate",
					fitness_effects.max_relative_migration_rate);

  pt.add("parameters.non_biological_parameters.init_pop",
				  non_biological_parameters.init_pop);
  pt.add("parameters.non_biological_parameters.matrix_max",
					non_biological_parameters.matrix_max);
  pt.add("parameters.non_biological_parameters.max_pop",
					non_biological_parameters.max_pop);
  pt.add("parameters.non_biological_parameters.max_time",
					non_biological_parameters.max_time);
  pt.add("parameters.non_biological_parameters.max_generations",
					non_biological_parameters.max_generations);
  pt.add("parameters.non_biological_parameters.seed",
					non_biological_parameters.seed);
  pt.add("parameters.non_biological_parameters.write_grid",
				  non_biological_parameters.write_grid);
  pt.add("parameters.non_biological_parameters.write_clones_file",
				  non_biological_parameters.write_clones_file);
  pt.add("parameters.non_biological_parameters.write_demes_file",
				  non_biological_parameters.write_demes_file);
  pt.add("parameters.non_biological_parameters.record_matrix",
				  non_biological_parameters.record_matrix);
  pt.add("parameters.non_biological_parameters.write_phylo",
				  non_biological_parameters.write_phylo);
  pt.add("parameters.non_biological_parameters.calculate_total_diversity",
				  non_biological_parameters.calculate_total_diversity);
  pt.add("parameters.non_biological_parameters.biopsy_size_per_sample",
				  non_biological_parameters.biopsy_size_per_sample);

  boost::property_tree::info_parser::write_info(file_name, pt);


};

boost::property_tree::ptree Parameters::getPTree() {
	boost::property_tree::ptree pt;

  pt.add("parameters.spatial_structure.log2_deme_carrying_capacity" ,
					spatial_structure.log2_deme_carrying_capacity);

  pt.add("parameters.dispersal.migration_type",
					dispersal.migration_type);
  pt.add("parameters.dispersal.init_migration_rate",
					dispersal.init_migration_rate);
  pt.add("parameters.dispersal.migration_edge_only",
					dispersal.migration_edge_only);
  pt.add("parameters.dispersal.migration_rate_scales_with_K",
					dispersal.migration_rate_scales_with_K);

  pt.add("parameters.mutation_rates.mu_driver_birth",
					mutation_rates.mu_driver_birth);
  pt.add("parameters.mutation_rates.mu_passenger",
				  mutation_rates.mu_passenger);
  pt.add("parameters.mutation_rates.mu_driver_migration",
				  mutation_rates.mu_driver_migration);
  pt.add("parameters.mutation_rates.passenger_pop_threshold",
				  mutation_rates.passenger_pop_threshold);

  pt.add("parameters.fitness_effects.normal_birth_rate",
				  fitness_effects.normal_birth_rate);
  pt.add("parameters.fitness_effects.baseline_death_rate",
					fitness_effects.baseline_death_rate);
  pt.add("parameters.fitness_effects.s_driver_birth",
					fitness_effects.s_driver_birth);
  pt.add("parameters.fitness_effects.s_passenger",
					fitness_effects.s_passenger);
  pt.add("parameters.fitness_effects.s_driver_migration",
					fitness_effects.s_driver_migration);
  pt.add("parameters.fitness_effects.max_relative_birth_rate",
					fitness_effects.max_relative_birth_rate);
  pt.add("parameters.fitness_effects.max_relative_migration_rate",
					fitness_effects.max_relative_migration_rate);

  pt.add("parameters.non_biological_parameters.init_pop",
				  non_biological_parameters.init_pop);
  pt.add("parameters.non_biological_parameters.matrix_max",
					non_biological_parameters.matrix_max);
  pt.add("parameters.non_biological_parameters.max_pop",
					non_biological_parameters.max_pop);
  pt.add("parameters.non_biological_parameters.max_time",
					non_biological_parameters.max_time);
  pt.add("parameters.non_biological_parameters.max_generations",
					non_biological_parameters.max_generations);
  pt.add("parameters.non_biological_parameters.seed",
					non_biological_parameters.seed);
  pt.add("parameters.non_biological_parameters.write_grid",
				  non_biological_parameters.write_grid);
  pt.add("parameters.non_biological_parameters.write_clones_file",
				  non_biological_parameters.write_clones_file);
  pt.add("parameters.non_biological_parameters.write_demes_file",
				  non_biological_parameters.write_demes_file);
  pt.add("parameters.non_biological_parameters.record_matrix",
				  non_biological_parameters.record_matrix);
  pt.add("parameters.non_biological_parameters.write_phylo",
				  non_biological_parameters.write_phylo);
  pt.add("parameters.non_biological_parameters.calculate_total_diversity",
				  non_biological_parameters.calculate_total_diversity);
  pt.add("parameters.non_biological_parameters.biopsy_size_per_sample",
				  non_biological_parameters.biopsy_size_per_sample);

	return pt;
};

void Parameters::GetFromPTree(boost::property_tree::ptree pt){

    spatial_structure.log2_deme_carrying_capacity =
			pt.get <int> ("parameters.spatial_structure.log2_deme_carrying_capacity");

    dispersal.migration_type =
			pt.get <int> ("parameters.dispersal.migration_type");
    dispersal.init_migration_rate =
			pt.get <float> ("parameters.dispersal.init_migration_rate");
    dispersal.migration_edge_only  =
			pt.get <int> ("parameters.dispersal.migration_edge_only");
    dispersal.migration_rate_scales_with_K  =
			pt.get <int> ("parameters.dispersal.migration_rate_scales_with_K");

		mutation_rates.mu_driver_birth =
			pt.get <float> ("parameters.mutation_rates.mu_driver_birth");
		mutation_rates.mu_passenger =
			pt.get <float> ("parameters.mutation_rates.mu_passenger");
		mutation_rates.mu_driver_migration =
			pt.get <float> ("parameters.mutation_rates.mu_driver_migration");
		mutation_rates.passenger_pop_threshold =
			pt.get <int> ("parameters.mutation_rates.passenger_pop_threshold");

		fitness_effects.normal_birth_rate =
			pt.get <float> ("parameters.fitness_effects.normal_birth_rate");
		fitness_effects.baseline_death_rate =
			pt.get <float> ("parameters.fitness_effects.baseline_death_rate");
		fitness_effects.s_driver_birth =
			pt.get <float> ("parameters.fitness_effects.s_driver_birth");
		fitness_effects.s_passenger =
			pt.get <float> ("parameters.fitness_effects.s_passenger");
		fitness_effects.s_driver_migration  =
			pt.get <float> ("parameters.fitness_effects.s_driver_migration");
		fitness_effects.max_relative_birth_rate  =
			pt.get <float> ("parameters.fitness_effects.max_relative_birth_rate");
		fitness_effects.max_relative_migration_rate  =
			pt.get <float> ("parameters.fitness_effects.max_relative_migration_rate");

		non_biological_parameters.init_pop =
			pt.get <int> ("parameters.non_biological_parameters.init_pop");
		non_biological_parameters.matrix_max =
			pt.get <int> ("parameters.non_biological_parameters.matrix_max");
		non_biological_parameters.max_pop =
			pt.get <int> ("parameters.non_biological_parameters.max_pop");
		non_biological_parameters.max_time =
			pt.get <int> ("parameters.non_biological_parameters.max_time");
		non_biological_parameters.max_generations =
			pt.get <int> ("parameters.non_biological_parameters.max_generations");
		non_biological_parameters.seed =
			pt.get <int> ("parameters.non_biological_parameters.seed");
		non_biological_parameters.write_grid =
			pt.get <int> ("parameters.non_biological_parameters.write_grid");
		non_biological_parameters.write_clones_file =
			pt.get <int> ("parameters.non_biological_parameters.write_clones_file");
		non_biological_parameters.write_demes_file =
			pt.get <int> ("parameters.non_biological_parameters.write_demes_file");
		non_biological_parameters.record_matrix =
			pt.get <int> ("parameters.non_biological_parameters.record_matrix");
		non_biological_parameters.write_phylo =
			pt.get <int> ("parameters.non_biological_parameters.write_phylo");
		non_biological_parameters.calculate_total_diversity =
			pt.get <int> ("parameters.non_biological_parameters.calculate_total_diversity");
		non_biological_parameters.biopsy_size_per_sample =
			pt.get <int> ("parameters.non_biological_parameters.biopsy_size_per_sample");

			p_tree = pt;
};

void Parameters::PrintToConsole() {
		std::cout<<"parameters.spatial_structure.log2_deme_carrying_capacity = "<<
								spatial_structure.log2_deme_carrying_capacity <<std::endl;

		std::cout<<"parameters.dispersal.migration_type = "<<
								dispersal.migration_type <<std::endl;
		std::cout<<"parameters.dispersal.init_migration_rate = "<<
								dispersal.init_migration_rate <<std::endl;
		std::cout<<"parameters.dispersal.migration_edge_only = "<<
								dispersal.migration_edge_only <<std::endl;
		std::cout<<"parameters.dispersal.migration_rate_scales_with_K = "<<
								dispersal.migration_rate_scales_with_K <<std::endl;


		std::cout<< "parameters.mutation_rates.mu_driver_birth = " <<
								 mutation_rates.mu_driver_birth <<std::endl;
		std::cout<<	"parameters.mutation_rates.mu_passenger = " <<
								 mutation_rates.mu_passenger << std::endl;
		std::cout<<	"parameters.mutation_rates.mu_driver_migration = " <<
								 mutation_rates.mu_driver_migration << std::endl;
		std::cout<<	"parameters.mutation_rates.passenger_pop_threshold = " <<
								 mutation_rates.passenger_pop_threshold << std::endl;

		std::cout<< "parameters.fitness_effects.normal_birth_rate = " <<
								 fitness_effects.normal_birth_rate<< std::endl;
		std::cout<< "parameters.fitness_effects.baseline_death_rate = " <<
								 fitness_effects.baseline_death_rate << std::endl;
		std::cout<< "parameters.fitness_effects.s_driver_birth = " <<
								 fitness_effects.s_driver_birth << std::endl;
		std::cout<< "parameters.fitness_effects.s_passenger = " <<
								 fitness_effects.s_passenger << std::endl;
		std::cout<< "parameters.fitness_effects.s_driver_migration = " <<
								 fitness_effects.s_driver_migration<< std::endl;
		std::cout<< "parameters.fitness_effects.max_relative_birth_rate = " <<
								 fitness_effects.max_relative_birth_rate<< std::endl;
		std::cout<< "parameters.fitness_effects.max_relative_migration_rate = " <<
								 fitness_effects.max_relative_migration_rate<< std::endl;

		std::cout<< "parameters.non_biological_parameters.init_pop = " <<
							   non_biological_parameters.init_pop << std::endl;
		std::cout<< "parameters.non_biological_parameters.matrix_max = " <<
								 non_biological_parameters.matrix_max << std::endl;
		std::cout<< "parameters.non_biological_parameters.max_pop = " <<
								 non_biological_parameters.max_pop << std::endl;
		std::cout<< "parameters.non_biological_parameters.max_time = " <<
								 non_biological_parameters.max_pop << std::endl;
		std::cout<< "parameters.non_biological_parameters.max_generations = " <<
								 non_biological_parameters.max_generations << std::endl;
		std::cout<< "parameters.non_biological_parameters.seed =" <<
								 non_biological_parameters.seed << std::endl;
		std::cout<< "parameters.non_biological_parameters.write_grid = " <<
								 non_biological_parameters.write_grid << std::endl;
		std::cout<< "parameters.non_biological_parameters.write_clones_file = " <<
								 non_biological_parameters.write_clones_file << std::endl;
		std::cout<< "parameters.non_biological_parameters.write_demes_file = " <<
								 non_biological_parameters.write_demes_file << std::endl;
		std::cout<< "parameters.non_biological_parameters.record_matrix = " <<
								 non_biological_parameters.record_matrix << std::endl;
		std::cout<< "parameters.non_biological_parameters.write_phylo = " <<
								 non_biological_parameters.write_phylo << std::endl;
		std::cout<< "parameters.non_biological_parameters.calculate_total_diversity = " <<
								 non_biological_parameters.calculate_total_diversity << std::endl;
		std::cout<< "parameters.non_biological_parameters.biopsy_size_per_sample = " <<
								 non_biological_parameters.biopsy_size_per_sample << std::endl;
};

