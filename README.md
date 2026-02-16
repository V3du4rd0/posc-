# Posc++
Library for Power-Series Composition

## Requirements and compilation

This library compiles with modern versions of `g++`. In particularly, it has been tested with gcc-14.2.0 and gcc-15.2.1.

The provided Makefile should avoid any linking issues. Perform a test compilation by running: 

```
make
```

A successfull compilation will create a directory called `bin` containing the binary files. This process creates an test example in the `bin` directory.

This library supports single, double and quadruple precission. Before compilation, make sure to select the desired precission in the `src/config.h` file.

