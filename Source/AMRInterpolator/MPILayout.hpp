/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to Copyright.txt in GRChombo's root directory.
 */

#ifndef MPILAYOUT_HPP_
#define MPILAYOUT_HPP_

class MPILayout
{

  public:
    // Getters
    inline int count(int rank) const;
    inline int totalCount() const;
    inline int displ(int rank) const;

    // Setters
    inline void setCount(int rank, int count);
    inline void incrementCount(int rank);
    inline void clearCounts();

  private:
    friend class MPIContext;

    MPILayout(int num_process);

    const int m_num_process;
    std::vector<int> m_counts;
    mutable std::vector<int> m_displs;

    mutable int m_total_count;
    mutable bool m_dirty;

    inline void updateDirty() const;
    inline int *countsPtr();
    inline int *displsPtr();
};

#include "MPILayout.impl.hpp"

#endif /* MPILAYOUT_HPP_ */
