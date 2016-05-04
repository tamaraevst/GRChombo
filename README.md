# GRChombo modernisation
## C++ification of RHS
`@saran` : I have rewritten the CCZ4 RHS in C++ and also attempted to vectorise it using SIMD templates.

Here's a quick guide for code navigation:
- `CCZ4` contains actual RHS, differentiation, and a driver that loops over `FArrayBox`. This code has explicit vectorisation.
- `CCZ4Test` contains the test case that compares the new RHS with the old `ChF` version. You can build an executable binary by issuing `make CCZ4Test` in here. (Also performs a rough-and-ready timing comparison between the two versions)
- `simd` contains the templated SIMD library. Supports AVX512, AVX and SSE right now.
- `utils` contains various useful, but as yet unorganised, `#define` and `struct`.
