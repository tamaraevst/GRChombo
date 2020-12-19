/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRUTILS_HPP_
#define GRUTILS_HPP_

#ifndef GR_SPACEDIM
#define GR_SPACEDIM 3
#endif

#ifndef DEFAULT_TENSOR_DIM
#define DEFAULT_TENSOR_DIM 3
#endif

#define FOR1(IDX) for (int IDX = 0; IDX < DEFAULT_TENSOR_DIM; ++IDX)
#define FOR2(IDX1, IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1, IDX2, IDX3) FOR2(IDX1, IDX2) FOR1(IDX3)
#define FOR4(IDX1, IDX2, IDX3, IDX4) FOR2(IDX1, IDX2) FOR2(IDX3, IDX4)
#define FOR5(IDX1, IDX2, IDX3, IDX4, IDX5) FOR2(IDX1, IDX2) FOR3(IDX3, IDX4, IDX5)
#define FOR6(IDX1, IDX2, IDX3, IDX4, IDX5, IDX6) FOR3(IDX1, IDX2, IDX3)
                                                         FOR3(IDX4, IDX5, IDX6)
#define FOR7(IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7) FOR4(IDX1, IDX2, IDX3, IDX4)
                                                          FOR3(IDX5, IDX6, IDX7)
#define FOR8(IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8)
                       FOR4(IDX1, IDX2, IDX3, IDX4) FOR4(IDX5, IDX6, IDX7, IDX8)

#endif /* GRUTILS_HPP_*/
