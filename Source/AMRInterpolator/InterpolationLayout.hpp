#ifndef _INTERPOLATIONLAYOUT_HPP_
#define _INTERPOLATIONLAYOUT_HPP_

class InterpolationLayout
{
  private:
    template <typename InterpAlgo> friend class AMRInterpolator;

    vector<int> rank;
    vector<int> level_idx;
    vector<int> box_idx;

    InterpolationLayout(int num_points)
        : rank(num_points, -1), level_idx(num_points, -1),
          box_idx(num_points, -1)
    {
    }
};

#endif /* _INTERPOLATIONLAYOUT_HPP_ */
