# owl

BGV key switching [[2023/1642](https://eprint.iacr.org/2023/1642)]


## files

- `build.sh`: build script
- `ring.h`: polynomial ring interface
- `ring_hexl.cpp`: polynomial ring using HEXL
- `swk.c`: key switching variants
- `test.c`: tests
- `bench.c`: benchmarks
- `params.py`: benchmarking parameter generation
- `params.c`: generated benchmarking parameters
- `latex.c`: benchmarking data to LaTeX conversion
- `data/`: benchmarking data used in the paper


## dependencies

- HEXL (included as submodule)
- C/C++ compiler (defaults to gcc/g++)
- CMake for building HEXL
- ssh/scp for remote building/testing


## build & test

```
# local
./build.sh build tests

# remote via ssh
./build.sh -r <remote> rcopy build tests
```

The remote build copies the files to `${HOME}/owl`, then builds and runs the tests on the remote.
If building HEXL fails, either run `./build.sh clean` or remove `build/hexl` before trying again.


## benchmarks

Running the HEXL benchmark `build/bench_hexl` prints the benchmarking results.
It also creates a `build/bench_hexl.dat`; `build/latex` converts this data file to LaTeX.
