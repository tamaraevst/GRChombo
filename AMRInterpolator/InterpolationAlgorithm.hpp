#ifndef _INTERPOLATIONALGORITHM_HPP_
#define _INTERPOLATIONALGORITHM_HPP_

class InterpolationAlgorithm {
public:
    virtual ~InterpolationAlgorithm() = 0;
};

class NearestNeighbour : public InterpolationAlgorithm {
public:
    static inline double interpPoint(const std::array<double, CH_SPACEDIM>& gridCoord, const FArrayBox& fab, int comps, const IntVect& nearest)
    {
        return fab.get(nearest, comps);
    }
};

#endif
