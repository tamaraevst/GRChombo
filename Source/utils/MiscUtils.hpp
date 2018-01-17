#ifndef MISCUTILS_H_
#define MISCUTILS_H_

#define MIN_CUT_OFF(VAR, MIN_VAL)                                              \
    VAR = simd_conditional(simd_compare_lt(VAR, MIN_VAL), MIN_VAL, VAR)

#endif /* MISCUTILS_H_*/
