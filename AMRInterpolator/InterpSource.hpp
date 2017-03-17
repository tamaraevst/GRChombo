#ifndef _INTERPSOURCE_H_
#define _INTERPSOURCE_H_

#include "Array.hpp"

// Abstrace base class to get the FABs out of an AMRLevel
class InterpSource {
public:
    virtual const LevelData<FArrayBox>& getLevelData() const = 0;
    virtual bool contains(const Array<double, CH_SPACEDIM>& point) const = 0;
    virtual void fillAllGhosts() = 0;
};

#endif /* _INTERPSOURCE_H_ */
