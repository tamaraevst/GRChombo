/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "computeSum.H"
#include "computeNorm.H"
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>

Vector<LevelData<FArrayBox>* > GRAMR::getLevelDataPtrs()
{
    // First get a std::vector of AMRLevel pointers
    const std::vector<AMRLevel*> level_ptrs = getAMRLevels().stdVector();

    // Instatiate a Vector of LevelData<FArrayBox> pointers as this is the
    // the format that Chombo's computeSum requires
    Vector<LevelData<FArrayBox>*> level_data_ptrs(level_ptrs.size());

    // Now get the level data pointers
    std::transform(level_ptrs.begin(), level_ptrs.end(),
                   level_data_ptrs.stdVector().begin(),
                   [](AMRLevel* level_ptr)
                   {
                       return const_cast<GRLevelData*>
                           (&(GRAMRLevel::gr_cast(level_ptr)->getLevelData()));
                   });

    return level_data_ptrs;
}

Real GRAMR::compute_sum(const int a_comp, const Real a_dx_coarse)
{
    const Vector<LevelData<FArrayBox>* > level_data_ptrs {getLevelDataPtrs()};
    return computeSum(level_data_ptrs, m_ref_ratios, a_dx_coarse,
                      Interval(a_comp, a_comp), 0);
}

Real GRAMR::compute_norm(const int a_comp, const double a_p,
                        const Real a_dx_coarse)
{
    const Vector<LevelData<FArrayBox>* > level_data_ptrs {getLevelDataPtrs()};
    return computeNorm(level_data_ptrs, m_ref_ratios, a_dx_coarse,
                       Interval(a_comp, a_comp), a_p, 0);
}

Real GRAMR::compute_max(const Interval a_comps)
{
    const Vector<LevelData<FArrayBox>* > level_data_ptrs {getLevelDataPtrs()};
    return computeMax(level_data_ptrs, m_ref_ratios, a_comps, 0);
}

Real GRAMR::compute_min(const Interval a_comps)
{
    const Vector<LevelData<FArrayBox>* > level_data_ptrs {getLevelDataPtrs()};
    return computeMin(level_data_ptrs, m_ref_ratios, a_comps, 0);
}

Real GRAMR::compute_inf_norm(const Interval a_comps)
{
    return std::max(std::abs(compute_max(a_comps)),
                    std::abs(compute_min(a_comps)));
}
