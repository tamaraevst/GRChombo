/**
  * BINARY BLACK HOLES
  * Assuming that initial A_{ij} is approximately linear
  */

template <typename T>
class BinaryBH
{

public:
	const T& bh1;
	const T& bh2;

	BinaryBH(T& bh1, T& bh2);
	Real psi(Real x, Real y, Real z) const;
	void Aij(Real x, Real y, Real z, Real out[3][3]) const;

};

template <typename T>
BinaryBH<T>::BinaryBH(T& bh1, T& bh2) :
	bh1 (bh1),
	bh2 (bh2)
{

}

template <typename T>
Real
BinaryBH<T>::psi(Real x, Real y, Real z) const
{
	return 1 + bh1.psi_minus_one(x, y, z) + bh2.psi_minus_one(x, y, z);
}

template <typename T>
void
BinaryBH<T>::Aij(Real x, Real y, Real z, Real out[3][3]) const
{
	Real Aij1[3][3] = {{0}};
	bh1.Aij(x, y, z, Aij1);

	Real Aij2[3][3] = {{0}};
	bh2.Aij(x, y, z, Aij2);

	for (int i = 0; i < 3; ++i)
	{
		for (int j = i; j < 3; ++j)
		{
			out[i][j] = Aij1[i][j] + Aij2[i][j];
		}
	}
}
