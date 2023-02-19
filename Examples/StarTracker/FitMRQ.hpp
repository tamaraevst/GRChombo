#ifndef FITMRQ_HPP_
#define FITMRQ_HPP_

#include "Matrix.hpp"
#include <stdexcept>

class Fitmrq 
{
	public:

	static const int NDONE=4, ITMAX=5000;
	int ndat, ma, mfit;
	std::vector<double> &x,&y,&sig;
	double tol;
	void (*funcs)(const double, std::vector<double> &, double &, std::vector<double> &);
	std::vector<bool> ia;
	std::vector<double> a;
	double chisq;

	Fitmrq(std::vector<double> &xx, std::vector<double> &yy, std::vector<double> &ssig, std::vector<double> &aa,
	void funks(const double, std::vector<double> &, double &, std::vector<double> &), const double
	TOL=1.e-6) : ndat(xx.size()), ma(aa.size()), x(xx), y(yy), sig(ssig),
	tol(TOL), funcs(funks), ia(ma), a(aa) {
		for (int i=0;i<ma;i++) ia[i] = true;
	}

	int fit() {
		//std::cout<<"Entering the fit function" << std::endl;
		int j,k,l,iter,done=0;
		double alamda=.001,ochisq;
		//std::cout<<"Creating the vectors"<<std::endl;
		std::vector<double> atry(ma),beta(ma),da(ma);
		Matrix alpha(ma, ma), covar(ma, ma);
		mfit=0;
		//std::cout <<"Filling ia" <<std::endl;
		for (j=0;j<ma;j++) if (ia[j]) mfit++;
		Matrix oneda(mfit, mfit), temp(mfit, mfit);
		//std::cout<<"Calling MRQCOF" <<std::endl;
		mrqcof(a,alpha,beta);
		//std::cout<<"Filling atry-vector !"<<std::endl;
		for (j=0;j<ma;j++) atry[j]=a[j];
		ochisq=chisq;
		//std::cout<<"Filling Covar Mat" <<std::endl;
		for (iter=0;iter<ITMAX;iter++) {
			if (done==NDONE) alamda=0.;
			for (j=0;j<mfit;j++) {
				for (k=0;k<mfit;k++) covar.At(j,k)=alpha.At(j,k);
				covar.At(j,j)=alpha.At(j,j)*(1.0+alamda);
				for (k=0;k<mfit;k++) temp.At(j,k)=covar.At(j,k);
				oneda.At(j,0)=beta[j];
			}
			//std::cout<<"Calling Gauss"<<std::endl;
			gaussj(temp,oneda);
			for (j=0;j<mfit;j++) {
				for (k=0;k<mfit;k++) covar.At(j,k)=temp.At(j,k);
				da[j]=oneda.At(j,0);
			}
			if (done==NDONE) {
				covsrt(covar);
				covsrt(alpha);
				return;
			}
			for (j=0,l=0;l<ma;l++)
				if (ia[l]) atry[l]=a[l]+da[j++];
			//std::cout<<"Calling MRQCOF" <<std::endl;
			mrqcof(atry,covar,da);
			if (abs(chisq-ochisq) < MAX(tol,tol*chisq)) done++;
			if (chisq < ochisq) {
				alamda *= 0.1;
				ochisq=chisq;
				for (j=0;j<mfit;j++) {
					for (k=0;k<mfit;k++) alpha.At(j,k)=covar.At(j,k);
						beta[j]=da[j];
				}
				for (l=0;l<ma;l++) a[l]=atry[l];
			} else {
				alamda *= 10.0;
				chisq=ochisq;
			}
		}
		//throw std::runtime_error("Fitmrq too many iterations");
		return 0;
	}

	void hold(const int i, const double val) {ia[i]=false; a[i]=val;}
	void free(const int i) {ia[i]=true;}

	void mrqcof(std::vector<double> &a, Matrix &alpha, std::vector<double> &beta) {
		//std::cout<<"Enetering MRQCOF"<<std::endl;
		int i,j,k,l,m;
		double ymod,wt,sig2i,dy;
		std::vector<double> dyda(ma);
		//std::cout<<"Generating alpha"<<std::endl;
		for (j=0;j<mfit;j++) {
			for (k=0;k<=j;k++) alpha.At(j,k)=0.0;
			beta[j]=0.;
		}
		chisq=0.;
		//                std::cout<<"Generating sigma"<<std::endl;

		for (i=0;i<ndat;i++) {
			funcs(x[i],a,ymod,dyda);
			sig2i=1.0/(sig[i]*sig[i]);
			dy=y[i]-ymod;
			for (j=0,l=0;l<ma;l++) {
				if (ia[l]) {
					wt=dyda[l]*sig2i;
					for (k=0,m=0;m<l+1;m++)
						if (ia[m]) alpha.At(j,k++) += wt*dyda[m];
					beta[j++] += dy*wt;
				}
			}
			chisq += dy*dy*sig2i;
		}
		for (j=1;j<mfit;j++)
			for (k=0;k<j;k++) alpha.At(k,j)=alpha.At(j,k);
	}

	void covsrt(Matrix &covar) {
		//std::cout<<"Entering COVRST"<<std::endl;
		int i,j,k;
		for (i=mfit;i<ma;i++)
			for (j=0;j<i+1;j++) covar.At(i,j)=covar.At(j,i)=0.0;
		k=mfit-1;
		for (j=ma-1;j>=0;j--) {
			if (ia[j]) {
				for (i=0;i<ma;i++) SWAP(covar.At(i,k),covar.At(i,j));
				for (i=0;i<ma;i++) SWAP(covar.At(k,i),covar.At(j,i));
				k--;
			}
		}
	}

	void gaussj(Matrix &a, Matrix &b)
	{
		int i,icol,irow,j,k,l,ll,n=a.nrows(),m=b.ncols();
		double big,dum,pivinv;
		std::vector<double> indxc(n),indxr(n),ipiv(n);
		for (j=0;j<n;j++) ipiv[j]=0;
		for (i=0;i<n;i++) {
			big=0.0;
			for (j=0;j<n;j++)
				if (ipiv[j] != 1)
					for (k=0;k<n;k++) {
						if (ipiv[k] == 0) {
							if (abs(a.At(j,k)) >= big) {
								big=abs(a.At(j,k));
								irow=j;
								icol=k;
							}
						}
					}
			++(ipiv[icol]);
			if (irow != icol) {
				for (l=0;l<n;l++) SWAP(a.At(irow,l),a.At(icol,l));
				for (l=0;l<m;l++) SWAP(b.At(irow,l),b.At(icol,l));
			}
			indxr[i]=irow;
			indxc[i]=icol;
			if (a.At(icol,icol) == 0.0) throw("gaussj: Singular Matrix");
			pivinv=1.0/a.At(icol,icol);
			a.At(icol,icol)=1.0;
			for (l=0;l<n;l++) a.At(icol,l) *= pivinv;
			for (l=0;l<m;l++) b.At(icol,l) *= pivinv;
			for (ll=0;ll<n;ll++)
				if (ll != icol) {
					dum=a.At(ll,icol);
					a.At(ll,icol)=0.0;
					for (l=0;l<n;l++) a.At(ll,l) -= a.At(icol,l)*dum;
					for (l=0;l<m;l++) b.At(ll,l) -= b.At(icol,l)*dum;
				}
		}
		for (l=n-1;l>=0;l--) {
			if (indxr[l] != indxc[l])
				for (k=0;k<n;k++)
					SWAP(a.At(k,indxr[l]),a.At(k,indxc[l]));
		}
	}

};

#endif /* FITMRQ_HPP_ */
