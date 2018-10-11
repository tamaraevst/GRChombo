/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DERIVATIVESETUP_HPP_
#define DERIVATIVESETUP_HPP_

inline const Derivative Derivative::LOCAL;

inline const Derivative Derivative::dx(0);
inline const Derivative Derivative::dy(1);
inline const Derivative Derivative::dz(2);

inline const Derivative Derivative::dxdx(0, 0);
inline const Derivative Derivative::dydy(1, 1);
inline const Derivative Derivative::dzdz(2, 2);

inline const Derivative Derivative::dxdy(0, 1);
inline const Derivative Derivative::dxdz(0, 2);
inline const Derivative Derivative::dydz(1, 2);

#endif /* DERIVATIVESETUP_HPP_ */
