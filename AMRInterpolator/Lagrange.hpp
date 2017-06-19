#ifndef _LAGRANGE_HPP_
#define _LAGRANGE_HPP_

#include <utility>
#include "InterpSource.hpp"

template <int Order>
class Lagrange
{
    const InterpSource& m_source;
    bool m_verbosity;

    struct Stencil;

    struct Stencil
    {
        int m_width;
        int m_deriv;
        double m_dx;
        double m_point_offset;

        double *m_weights;

        Stencil(int width, int deriv, double point_offset, double dx);
        inline bool operator==(const Stencil& rhs) const;
        inline bool isSameAs(int width, int deriv, double point_offset, double dx) const;

        inline const double& operator[](unsigned int i) const;
    };

    typedef vector<Stencil> stencil_collection_t;
    stencil_collection_t m_memoized_stencils;

    Stencil getStencil(int width, int deriv, double point_offset, double dx);

    // Helper function to generate tensor product weights
    // Argument 'dim' is used for recursion over dimensions.
    pair<vector<IntVect>, vector<double> >
    generateStencil(const std::array<int, CH_SPACEDIM>& deriv, const std::array<double, CH_SPACEDIM>& dx, const std::array<double, CH_SPACEDIM>& evalCoord, const IntVect& nearest, int dim = CH_SPACEDIM - 1);

    vector<IntVect> m_interp_points;
    vector<double> m_interp_weights;

    // We are adding 216+ numbers at roughly the same magnitudes but alternating signs.
    // Let's keep track of positive and negative terms separately to make sure we don't run into trouble.
    multiset<double> m_interp_neg;
    multiset<double> m_interp_pos;

public:
    Lagrange(const InterpSource& source, bool verbosity = false);

    void setup(const std::array<int, CH_SPACEDIM>& deriv, const std::array<double, CH_SPACEDIM>& dx, const std::array<double, CH_SPACEDIM>& evalCoord, const IntVect& nearest);
    double interpData(const FArrayBox& fab, int comp);

    const static string TAG;
};

#include "Lagrange.impl.hpp"

#endif /* _LAGRANGE_HPP_ */
