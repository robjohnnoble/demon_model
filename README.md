# methdemon
MethDemon (fluctuating **meth**ylation clock **dem**e-based **on**cology model in 1D) is a proprietary framework for simulating fluctuating CpG loci in a cancer population expanding via gland fission, of which only up to 8 glands are tracked and sampled.

## Prerequisites

The program is written in C++; apart from standard libraries, it requires the Boost C++ library.

## Installing and compiling

```
git clone -b methdemon1d https://github.com/vesmanojlovic/methdemon
cd methdemon1d
```
Then modify the `Makefile` to include the path to your copy of `boost` and run
```
make all
```

Create a directory where you would like your outputs stored and make a copy of the `config.dat` file modified with your required parameters. From the repo parent directory run
```
bin/methdemon <output_dir_path> <config_file_name>
```

To clear logfiles and binaries run
```
make clean
```
