#ifndef _QUINTICCONVOLUTION_HPP_
#define _QUINTICCONVOLUTION_HPP_

#include "InterpSource.hpp"
#include <utility>

class QuinticConvolution
{
    const InterpSource& m_source;
    bool m_verbosity;

    vector<IntVect> m_interp_points;
    vector<double> m_interp_weights;

public:
    QuinticConvolution(const InterpSource& source, bool verbosity = false);

    void setup(const Array<int, CH_SPACEDIM>& deriv, const Array<double, CH_SPACEDIM>& dx, const Array<double, CH_SPACEDIM>& evalCoord, const IntVect& nearest);
    double interpData(const FArrayBox& fab, int comp);

    const static string TAG;
};

#include "QuinticConvolution.impl.hpp"

#endif /* _QUINTICCONVOLUTION_HPP_ */
