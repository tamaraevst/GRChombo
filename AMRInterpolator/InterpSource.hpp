#ifndef _INTERPSOURCE_HPP_
#define _INTERPSOURCE_HPP_

#include "Array.hpp"

// Abstrace base class to get the FABs out of an AMRLevel
class InterpSource {
public:
    virtual const LevelData<FArrayBox>& getLevelData() const = 0;
    virtual bool contains(const Array<double, CH_SPACEDIM>& point) const = 0;
    virtual void refresh() = 0;
};

#endif /* _INTERPSOURCE_HPP_ */
