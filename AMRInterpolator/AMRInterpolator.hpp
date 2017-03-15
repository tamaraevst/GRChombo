#ifndef _AMRINTERPOLATOR_HPP_
#define _AMRINTERPOLATOR_HPP_

// Chombo includes

#include "AMRLevel.hpp"
#include "AMR.hpp"

#include "UsingNamespace.hpp"

// Our includes

#include "Array.hpp"
#include "InterpolationLayout.hpp"
#include "InterpSource.hpp"
#include "InterpolationAlgorithm.hpp"
#include "MPIContext.hpp"
#include "InterpolationQuery.hpp"

// End include

template <typename InterpAlgo>
class AMRInterpolator
{
public:
    AMRInterpolator(const AMR& amr, const Array<double, CH_SPACEDIM>& coarsest_origin, const Array<double, CH_SPACEDIM>& coarsest_dx, int verbosity = 0);
    void refresh();
    void limit_num_levels(unsigned int num_levels);
    void interp(InterpolationQuery& query);

private:
    void computeLevelLayouts();
    InterpolationLayout findBoxes(InterpolationQuery& query);
    void prepareMPI(InterpolationQuery& query, const InterpolationLayout layout);
    void exchangeMPIQuery();
    void calculateAnswers(InterpolationQuery& query);
    void exchangeMPIAnswer();

    const AMR& m_amr;

    // Coordinates of the point represented by IntVect::Zero in coarsest grid
    const Array<double, CH_SPACEDIM> m_coarsest_origin;

    // Grid spacing in each direction
    const Array<double, CH_SPACEDIM> m_coarsest_dx;

    int m_num_levels;
    const int m_verbosity;

    vector<Array<double, CH_SPACEDIM> > m_origin;
    vector<Array<double, CH_SPACEDIM> > m_dx;

    MPIContext m_mpi;
    vector<int> m_mpi_mapping;

    // Memoisation of boxes previously found
    vector<int> m_mem_level;
    vector<int> m_mem_box;

    vector<int> m_query_level;
    vector<int> m_query_box;
    vector<double> m_query_coords[CH_SPACEDIM];
    vector<vector<double> > m_query_data;

    vector<int> m_answer_level;
    vector<int> m_answer_box;
    vector<double> m_answer_coords[CH_SPACEDIM];
    vector<vector<double> > m_answer_data;

    // A bit of Android-ism here, but it's really useful!
    // Identifies the printout as originating from this class.
    const static string TAG;

};

#include "AMRInterpolator.impl.hpp"

#endif /* _AMRINTERPOLATOR_HPP_ */
