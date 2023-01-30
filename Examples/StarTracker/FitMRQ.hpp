#ifndef FITMRQ_HPP_
#define FITMRQ_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

class Fitmrq 
{
	public:

	static const int NDONE=4, ITMAX=1000;
	int ndat, ma, mfit;
	std::vector<double> &x,&y,&sig;
	double tol;
	void (*funcs)(const double, std::vector<double> &, double &, std::vector<double> &);
	std::vector<bool> ia;
	std::vector<double> a;
	double chisq;

	Fitmrq(std::vector<double> &xx, std::vector<double> &yy, std::vector<double> &ssig, std::vector<double> &aa,
	void funks(const double, std::vector<double> &, double &, std::vector<double> &), const double
	TOL=1.e-3) : ndat(xx.size()), ma(aa.size()), x(xx), y(yy), sig(ssig),
	tol(TOL), funcs(funks), ia(ma), a(aa) {
		for (int i=0;i<ma;i++) ia[i] = true;
	}

	template <class data_t> void fit() {
		int j,k,l,iter,done=0;
		double alamda=.001,ochisq;
		std::vector<double> atry(ma),beta(ma),da(ma);
		Tensor<2, data_t> alpha, covar;
		mfit=0;
		for (j=0;j<ma;j++) if (ia[j]) mfit++;
		Tensor<2, data_t> oneda, temp;
		mrqcof(a,alpha,beta);
		for (j=0;j<ma;j++) atry[j]=a[j];
		ochisq=chisq;
		for (iter=0;iter<ITMAX;iter++) {
			if (done==NDONE) alamda=0.;
			for (j=0;j<mfit;j++) {
				for (k=0;k<mfit;k++) covar[j][k]=alpha[j][k];
				covar[j][j]=alpha[j][j]*(1.0+alamda);
				for (k=0;k<mfit;k++) temp[j][k]=covar[j][k];
				oneda[j][0]=beta[j];
			}
			gaussj(temp,oneda);
			for (j=0;j<mfit;j++) {
				for (k=0;k<mfit;k++) covar[j][k]=temp[j][k];
				da[j]=oneda[j][0];
			}
			if (done==NDONE) {
				covsrt(covar);
				covsrt(alpha);
				return;
			}
			for (j=0,l=0;l<ma;l++)
				if (ia[l]) atry[l]=a[l]+da[j++];
			mrqcof(atry,covar,da);
			if (abs(chisq-ochisq) < MAX(tol,tol*chisq)) done++;
			if (chisq < ochisq) {
				alamda *= 0.1;
				ochisq=chisq;
				for (j=0;j<mfit;j++) {
					for (k=0;k<mfit;k++) alpha[j][k]=covar[j][k];
						beta[j]=da[j];
				}
				for (l=0;l<ma;l++) a[l]=atry[l];
			} else {
				alamda *= 10.0;
				chisq=ochisq;
			}
		}
		throw("Fitmrq too many iterations");
	}

	void hold(const int i, const double val) {ia[i]=false; a[i]=val;}
	void free(const int i) {ia[i]=true;}

	template <class data_t> void mrqcof(std::vector<double> &a, Tensor<2, data_t> &alpha, std::vector<double> &beta) {
		int i,j,k,l,m;
		double ymod,wt,sig2i,dy;
		std::vector<double> dyda(ma);

		for (j=0;j<mfit;j++) {
			for (k=0;k<=j;k++) alpha[j][k]=0.0;
			beta[j]=0.;
		}
		chisq=0.;
		for (i=0;i<ndat;i++) {
			funcs(x[i],a,ymod,dyda);
			sig2i=1.0/(sig[i]*sig[i]);
			dy=y[i]-ymod;
			for (j=0,l=0;l<ma;l++) {
				if (ia[l]) {
					wt=dyda[l]*sig2i;
					for (k=0,m=0;m<l+1;m++)
						if (ia[m]) alpha[j][k++] += wt*dyda[m];
					beta[j++] += dy*wt;
				}
			}
			chisq += dy*dy*sig2i;
		}
		for (j=1;j<mfit;j++)
			for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];
	}

	template <class data_t> void covsrt(Tensor<2, data_t> &covar) {
		int i,j,k;
		for (i=mfit;i<ma;i++)
			for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
		k=mfit-1;
		for (j=ma-1;j>=0;j--) {
			if (ia[j]) {
				for (i=0;i<ma;i++) SWAP(covar[i][k],covar[i][j]);
				for (i=0;i<ma;i++) SWAP(covar[k][i],covar[j][i]);
				k--;
			}
		}
	}

};

#endif /* FITMRQ_HPP_ */