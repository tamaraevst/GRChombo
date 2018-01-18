#ifndef MISCUTILS_H_
#define MISCUTILS_H_

#define FOR1(IDX) for (int IDX = 0; IDX < DEFAULT_TENSOR_DIM; ++IDX)
#define FOR2(IDX1, IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1, IDX2, IDX3) FOR2(IDX1, IDX2) FOR1(IDX3)
#define FOR4(IDX1, IDX2, IDX3, IDX4) FOR2(IDX1, IDX2) FOR2(IDX3, IDX4)

#define MIN_CUT_OFF(VAR, MIN_VAL)                                              \
    VAR = simd_conditional(simd_compare_lt(VAR, MIN_VAL), MIN_VAL, VAR)

#endif /* MISCUTILS_H_*/
