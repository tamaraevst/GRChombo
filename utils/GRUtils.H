#ifndef _GRUTILS_H_
#define _GRUTILS_H_

#define GR_SPACEDIM 3
#define IDX_SPACEDIM 3

#define FOR1(IDX) for (int IDX = 0; IDX < IDX_SPACEDIM; ++IDX)
#define FOR2(IDX1,IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1,IDX2,IDX3) FOR2(IDX1,IDX2) FOR1(IDX3)
#define FOR4(IDX1,IDX2,IDX3,IDX4) FOR2(IDX1,IDX2) FOR2(IDX3,IDX4)

// C++ standard, 12.8 Copying class objects:
// Each subobject is assigned in the manner appropriate to its type:
// - if the subobject is an array, EACH ELEMENT IS ASSIGNED, in the manner appropriate to the element type;
template <int N, class data_t>
class tensor
{
	template <int, class> friend class tensor;
	typedef typename tensor<N-1, data_t>::arr_t arr_t[IDX_SPACEDIM];
	arr_t arr;

public:
	operator arr_t& ()
	{
		return arr;
	}

	operator const arr_t& () const
	{
		return arr;
	}
};

template <class data_t>
class tensor<0, data_t>
{
	template <int, class> friend class tensor;
	typedef data_t arr_t;
	arr_t arr;

public:
	operator arr_t& ()
	{
		return arr;
	}

	operator const arr_t& () const
	{
		return arr;
	}
};

template <int stencil_size>
struct stencil
{
	int offset[stencil_size];
	double weight[stencil_size];
};

#endif /* _GRUTILS_H_*/
