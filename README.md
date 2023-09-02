# A Faster Software Implementation of SQISign
## This code is an attachment of the manuscript "A Faster Software Implementation of SQISign", utilizing several techniques to optimize the implementation of SQISign. 

The code is based on the [code](https://github.com/SQISign/sqisign-ec23) proposed by Luca De Feo et al. For more details of their implementation see the following original readme.

The default compilation option also adaptes the precomputation technique as we proposed in Section 5.3 of our manuscript. It may enlarge the cost of key generation, but reduces the signing cost. 
To compile the code without the precomputation technique do
```
make NOPRECOMP=1
```

**The following is the original readme.**


# Optimized implementation of SQISign
## Accompanying implementation of the paper "New algorithms for the Deuring correspondence - Towards practical and secure SQISign signatures" (EUROCRYPT'23)

This code implements an efficient variant of the isogeny-based signature scheme SQISign.
For details see this [paper](https://eprint.iacr.org/2022/234).

(C) MIT license.

## Dependencies

The code depends on the latest stable version of the [PARI/GP
library](http://pari.math.u-bordeaux.fr/), 2.11.4.

The code has an optional dependency on [GMP](https://gmplib.org/),
which is also an optional dependency of PARI/GP and is typically
installed along with it.

## Supported platforms

The code compiles and runs on Linux and MacOS.

It contains an implementation of the low-level arithmetic functions using handwritten assembly for the x86-64 platform,
  starting from Broadwell architectures.

## Compile

To compile, test and benchmark, run

```
make
make check
make benchmark
```

## Changing base field

There are two supported primes pXXXX: `p3923` and `p6983`.
The default prime is `p3923`. 
To switch prime do

```
make PRIME=p6983
make check PRIME=p6983
make benchmark PRIME=p6983
```

Object files are created in separate folder, there is no need to
recompile from scratch with `make -B` when you change prime.


## Running individual tests and benchmarks

All tests are found in the folder `build/pXXXX/test/`, all benchmarks in the
folder `build/pXXXX/bench/`. To run them type with the corresponding prime

```
./build/pXXXX/test/sqisign
./build/pXXXX/bench/keygen
./build/pXXXX/bench/sign
./build/pXXXX/bench/verif
```