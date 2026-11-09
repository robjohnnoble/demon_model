# methdemon
MethDemon (fluctuating **meth**ylation clock **dem**e-based **on**cology model in 1D) is a framework for simulating fluctuating CpG loci in a cancer population expanding via gland fission, of which only a predetermined number of glands are tracked and sampled.

## Prerequisites

The program is written in C++; apart from standard libraries, it requires the Boost C++ library (e.g. `sudo apt install libboost-dev`, `sudo pacman -S boost` or `brew install boost`).

## Installing and compiling

```
git clone https://github.com/vesmanojlovic/methdemon
cd methdemon
make all
```
If `boost` is not installed in a standard location, add `-I <path_to_boost_include>` to `CXXFLAGS` in the `Makefile`.

Create a directory where you would like your outputs stored and make a copy of the `config.dat` file modified with your required parameters. From the repo parent directory run
```
bin/methdemon <output_dir_path> <config_file_name>
```
For example, to run the first example config
```
bin/methdemon examples/eg1 config.dat
```
The output is written to `final_demes.csv` in the output directory, with one row per gland every few generations and a final snapshot at the end of the simulation.

To clear logfiles and binaries run
```
make clean
```

## Config parameters
`capacity`
- `deme_carrying_capacity`: number of stem cells per gland (K), a gland can fission once it is full

`dispersal`
- `init_migration_rate`: gland fission rate per cell
- `left_demes`, `right_demes`: number of sampled glands on each side of the tumour (`-1` for no constraint, in which case 8 glands are sampled)

`mutation`
- `mu_driver_birth`: driver mutation rate per cell division
- `mu_driver_migration`: mutation rate of drivers that increase the fission rate

`fitness`
- `baseline_death_rate`: cell death rate
- `s_driver_birth`: mean selective advantage of a driver
- `s_driver_migration`: mean relative increase in fission rate of a fission driver
- `max_relative_birth_rate`, `max_relative_migration_rate`: caps on the birth and fission rates of driver genotypes

`methylation`
- `meth_rate`, `demeth_rate`: methylation and demethylation rates per fCpG allele per cell division
- `fCpG_loci_per_cell`: number of fCpG loci (each cell has two alleles per locus)
- `manual_array`: path to a file with the initial methylation array (`2 * fCpG_loci_per_cell` values of 0 or 1, one per line), or `-1` for a random initial array

`stopping_conditions`
- `max_fissions`: expansion stops once the tracked glands have undergone this many fissions on average, so the tumour has roughly `2^max_fissions` glands; this also sets the probability that a fission is tracked
- `turnover`: length of the turnover phase after expansion, as a fraction of the expansion time (`0` for growth only)

`output_indicators`
- `write_demes_file`: write `final_demes.csv`

`rng_seed`
- `seed`: random seed

`normal_birth_rate`, `migration_rate_scales_with_K`, `max_time`, `max_generations`, `init_pop`, `fission_config` and `write_clones_file` are read from the config but not used in the current version.
