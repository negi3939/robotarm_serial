#include<iostream>
#include<string>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<math.h>
#include<unistd.h>
#include<pthread.h>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/StdVector>
#ifndef GCC3p3
#include <Eigen/LU>
#endif
#include <functional> //-std=c++11をつける
#include "mymath.h"

using namespace Eigen;

#ifndef GCC3p3
namespace Mymath{

	Funcvec::Funcvec(){}
	Funcvec::~Funcvec(){}
	VectorXd Funcvec::function(VectorXd x){return x;}

	double sign(double A){
    	if(A>0.0d) return 1.0d;
    	else if(A<0.0d) return -1.0d;
    	else return 0.0d;
	}

	double step(double A){
		if(A>0.0d) return 1.0d;
		else return 0.0d;
	}

	double sigmoid(double x,double alp){
		return 0.5d*(std::tanh(0.5d*alp*x)+1.0d);
	}
	VectorXd sigmoid(VectorXd x,double alp){
		VectorXd ans = VectorXd::Zero(x.size());
		for(int ii=0;ii<x.size();ii++){
			ans(ii) = 0.5d*(std::tanh(0.5d*alp*x(ii))+1.0d);
		}
		return ans;
	}

	#if defined(CALC_EIGEN)
	void eig(MatrixXd aA,MatrixXd bB,MatrixXd &eigenV,MatrixXd &eigenD){
		if((aA.rows()!=aA.cols())||(bB.rows()!=bB.cols())||(aA.rows()!=bB.rows())){std::cout << "Matrix size is not same."<<std::endl;}
		int ii,jj,size = aA.rows();
		GeneralizedEigenSolver<MatrixXd> solve(aA,bB);
   		MatrixXcd Dc = solve.eigenvalues();
   		VectorXd Dv = Dc.real();
   		eigenD = Dv.asDiagonal();
   		MatrixXcd Vc = solve.eigenvectors();
   		MatrixXd Vbuf = Vc.real();
		eigenV = MatrixXd::Zero(size,size);
		for(jj=0;jj<size;jj++){
			eigenV(size-1,jj) = pow(-1.0d,(double)(jj+1));
			for(ii=0;ii<size-1;ii++){
				eigenV(ii,jj) = Vbuf(ii,jj)/(Vbuf(size-1,jj)/eigenV(size-1,jj));
			}
		}
		for(ii=0;ii<size;ii++){
			VectorXd buff = aA*eigenV.block(0,ii,size,1) - eigenD(ii,ii)*bB*eigenV.block(0,ii,size,1);
			for(int jj=0;jj<size;jj++){
				if(abs(buff(jj))>0.1){
				std::cout << "eigen val or eigen vector is not correct." << std::endl; 
			}
			}
		}	
	}
	
	void eig(MatrixXd aA,MatrixXd &eigenV,MatrixXd &eigenD){
		if(aA.rows()!=aA.cols()){std::cout << "Matrix size is not same."<<std::endl;}
		int size = aA.rows();
		EigenSolver<MatrixXd> solve(aA);
   		MatrixXcd Dc = solve.eigenvalues();
   		VectorXd Dv = Dc.real();
   		eigenD = Dv.asDiagonal();
   		MatrixXcd Vc = solve.eigenvectors();
   		eigenV = Vc.real();
		for(int ii=0;ii<size;ii++){
			VectorXd buff = aA*eigenV.block(0,ii,size,1) - eigenD(ii,ii)*eigenV.block(0,ii,size,1);
			for(int jj=0;jj<size;jj++){
				if(abs(buff(jj))>0.1){
				std::cout << "eigen val or eigen vector is not correct." << std::endl; 
			}
			}
		}	
	}
	#endif

	MatrixXd inv(MatrixXd aA){
		FullPivLU< MatrixXd > invM(aA);
		return invM.inverse();
	}
	
	MatrixXd pseudo_inv(MatrixXd aA){
		return inv(aA.transpose()*aA)*aA.transpose();	
	}

	MatrixXd absmat(MatrixXd aA){
		int ii,jj;
		MatrixXd ans = MatrixXd::Zero(aA.rows(),aA.cols()); 
		for(ii=0;ii<aA.rows();ii++){
	   		for(jj=0;jj<aA.cols();jj++){
		  		ans(ii,jj) = std::abs(aA(ii,jj));  
   			}
   		}
		return ans;
	}

	MatrixXd diffvec(VectorXd x,Funcvec *func){
    	int ii;
    	VectorXd fx = func->function(x);
    	MatrixXd ans(fx.size(),x.size());
    	MatrixXd bef(fx.size(),x.size());
    	MatrixXd aft(fx.size(),x.size());
    	double delta = 0.000001d;
    	VectorXd deltax(x.size());
    	for(ii=0;ii<x.size();ii++){
        	deltax = VectorXd::Zero(x.size());
        	deltax(ii) =  delta; 
        	bef.block(0,ii,fx.size(),1) = (fx-func->function(x-deltax))/delta;
        	aft.block(0,ii,fx.size(),1) = (func->function(x+deltax)-fx)/delta;
    	}
    	ans = (bef+aft)/2.0;
    	return ans;
	}

	MatrixXd diffdiffvec(VectorXd x,Funcvec *func){
    	VectorXd fxv = func->function(x);
		double fx = fxv(0);
    	MatrixXd ans(x.size(),x.size());
    	double delta = 0.00001d;
		double deltainv = 1.0/delta;
    	VectorXd deltaxpp(x.size());
		VectorXd deltaxpm(x.size());
		VectorXd deltaxmp(x.size());
		VectorXd deltaxmm(x.size());
    	for(int ii=0;ii<x.size();ii++){
			for(int jj=0;jj<x.size();jj++){
				deltaxpp = VectorXd::Zero(x.size());
				deltaxpm = VectorXd::Zero(x.size());
				deltaxmp = VectorXd::Zero(x.size());
				deltaxmm = VectorXd::Zero(x.size());
				if(ii==jj){
					deltaxpp(ii) = delta;
					deltaxmm(ii) = -delta;
					ans(ii,jj) = 0.5*deltainv*deltainv*(2.0*func->function(x+deltaxpp)(0) + 2.0*func->function(x+deltaxmm)(0) -4.0*fx);
				}else{
					deltaxpp(ii) = delta;
					deltaxpp(jj) = delta;
					deltaxmm(ii) = -delta;
					deltaxmm(jj) = -delta;
					deltaxpm(ii) = delta;
					deltaxpm(jj) = -delta;
					deltaxmp(ii) = -delta;
					deltaxmp(jj) = delta;
					ans(ii,jj) = 0.25d*deltainv*deltainv*(func->function(x+deltaxpp)(0) + func->function(x+deltaxmm)(0) - func->function(x+deltaxpm)(0) - func->function(x+deltaxmp)(0));
				} 
			}
    	}
    	return ans;
	}

	MatrixXd diffnorm(VectorXd x,Funcvec *func){
    	int ii;
    	double fx = func->function(x).norm();
    	MatrixXd ans(1,x.size());
    	MatrixXd bef(1,x.size());
    	MatrixXd aft(1,x.size());
    	double delta = 0.000001d;
    	VectorXd deltax(x.size());
    	for(ii=0;ii<x.size();ii++){
        	deltax = VectorXd::Zero(x.size());
        	deltax(ii) =  delta; 
        	bef(0,ii) = (fx-func->function(x-deltax).norm())/delta;
        	aft(0,ii) = (func->function(x+deltax).norm()-fx)/delta;
    	}
    	ans = (bef+aft)/2.0;
    	return ans;
	}

	void showvec(VectorXd &a){
		for(int ii=0;ii<a.size()-1;ii++){
			std::cout << a(ii) << " , " << std::flush;
		}
		std::cout << a(a.size()-1) << std::endl;
	}
	void showvec(VectorXf &a){
		for(int ii=0;ii<a.size()-1;ii++){
			std::cout << a(ii) << " , " << std::flush;
		}
		std::cout << a(a.size()-1) << std::endl;
	}
	void showvec(Vector3d &a){
		for(int ii=0;ii<2;ii++){
			std::cout << a(ii) << " , " << std::flush;
		}
		std::cout << a(a.size()-1) << std::endl;
	}
	void showvec(Vector3f &a){
		for(int ii=0;ii<2;ii++){
			std::cout << a(ii) << " , " << std::flush;
		}
		std::cout << a(a.size()-1) << std::endl;
	}
}

#else
//GCC 3.3 用
namespace Mymath{

	Funcvec::Funcvec(){}
	Funcvec::~Funcvec(){}
	VectorXd Funcvec::function(VectorXd x){return x;}

	double sign(double A){
    	if(A>0.0) return 1.0;
    	else if(A<0.0) return -1.0;
    	else return 0.0;
	}

	double step(double A){
		if(A>0.0) return 1.0;
		else return 0.0;
	}

	double sigmoid(double x,double alp){
		return 0.5*(std::tanh(0.5*alp*x)+1.0);
	}
	VectorXd sigmoid(VectorXd x,double alp){
		VectorXd ans = VectorXd::Zero(x.size());
		for(int ii=0;ii<x.size();ii++){
			ans(ii) = 0.5*(std::tanh(0.5*alp*x(ii))+1.0);
		}
		return ans;
	}

	#if defined(CALC_EIGEN)
	void eig(MatrixXd aA,MatrixXd bB,MatrixXd &eigenV,MatrixXd &eigenD){
		if((aA.rows()!=aA.cols())||(bB.rows()!=bB.cols())||(aA.rows()!=bB.rows())){std::cout << "Matrix size is not same."<<std::endl;}
		int ii,jj,size = aA.rows();
		GeneralizedEigenSolver<MatrixXd> solve(aA,bB);
   		MatrixXcd Dc = solve.eigenvalues();
   		VectorXd Dv = Dc.real();
   		eigenD = Dv.asDiagonal();
   		MatrixXcd Vc = solve.eigenvectors();
   		MatrixXd Vbuf = Vc.real();
		eigenV = MatrixXd::Zero(size,size);
		for(jj=0;jj<size;jj++){
			eigenV(size-1,jj) = pow(-1.0,(double)(jj+1));
			for(ii=0;ii<size-1;ii++){
				eigenV(ii,jj) = Vbuf(ii,jj)/(Vbuf(size-1,jj)/eigenV(size-1,jj));
			}
		}
		for(ii=0;ii<size;ii++){
			VectorXd buff = aA*eigenV.block(0,ii,size,1) - eigenD(ii,ii)*bB*eigenV.block(0,ii,size,1);
			for(int jj=0;jj<size;jj++){
				if(abs(buff(jj))>0.1){
				std::cout << "eigen val or eigen vector is not correct." << std::endl; 
			}
			}
		}	
	}
	
	void eig(MatrixXd aA,MatrixXd &eigenV,MatrixXd &eigenD){
		if(aA.rows()!=aA.cols()){std::cout << "Matrix size is not same."<<std::endl;}
		int size = aA.rows();
		EigenSolver<MatrixXd> solve(aA);
   		MatrixXcd Dc = solve.eigenvalues();
   		VectorXd Dv = Dc.real();
   		eigenD = Dv.asDiagonal();
   		MatrixXcd Vc = solve.eigenvectors();
   		eigenV = Vc.real();
		for(int ii=0;ii<size;ii++){
			VectorXd buff = aA*eigenV.block(0,ii,size,1) - eigenD(ii,ii)*eigenV.block(0,ii,size,1);
			for(int jj=0;jj<size;jj++){
				if(abs(buff(jj))>0.1){
				std::cout << "eigen val or eigen vector is not correct." << std::endl; 
			}
			}
		}	
	}
	#endif

	MatrixXd inv(MatrixXd aA){
		return aA.inverse();
	}
	
	MatrixXd pseudo_inv(MatrixXd aA){
		return inv(aA.transpose()*aA)*aA.transpose();	
	}

	MatrixXd absmat(MatrixXd aA){
		int ii,jj;
		MatrixXd ans = MatrixXd::Zero(aA.rows(),aA.cols()); 
		for(ii=0;ii<aA.rows();ii++){
	   		for(jj=0;jj<aA.cols();jj++){
		  		ans(ii,jj) = std::abs(aA(ii,jj));  
   			}
   		}
		return ans;
	}

	MatrixXd diffvec(VectorXd x,Funcvec *func){
    	int ii;
    	VectorXd fx = func->function(x);
    	MatrixXd ans(fx.size(),x.size());
    	MatrixXd bef(fx.size(),x.size());
    	MatrixXd aft(fx.size(),x.size());
    	double delta = 0.000001;
    	VectorXd deltax(x.size());
    	for(ii=0;ii<x.size();ii++){
        	deltax = VectorXd::Zero(x.size());
        	deltax(ii) =  delta; 
        	bef.block(0,ii,fx.size(),1) = (fx-func->function(x-deltax))/delta;
        	aft.block(0,ii,fx.size(),1) = (func->function(x+deltax)-fx)/delta;
    	}
    	ans = (bef+aft)/2.0;
    	return ans;
	}

	MatrixXd diffnorm(VectorXd x,Funcvec *func){
    	int ii;
    	double fx = func->function(x).norm();
    	MatrixXd ans(1,x.size());
    	MatrixXd bef(1,x.size());
    	MatrixXd aft(1,x.size());
    	double delta = 0.000001;
    	VectorXd deltax(x.size());
    	for(ii=0;ii<x.size();ii++){
        	deltax = VectorXd::Zero(x.size());
        	deltax(ii) =  delta; 
        	bef(0,ii) = (fx-func->function(x-deltax).norm())/delta;
        	aft(0,ii) = (func->function(x+deltax).norm()-fx)/delta;
    	}
    	ans = (bef+aft)/2.0;
    	return ans;
	}

	void showvec(VectorXd &a){
		for(int ii=0;ii<a.size()-1;ii++){
			std::cout << a(ii) << " , " << std::flush;
		}
		std::cout << a(a.size()-1) << std::endl;
	}
	void showvec(VectorXf &a){
		for(int ii=0;ii<a.size()-1;ii++){
			std::cout << a(ii) << " , " << std::flush;
		}
		std::cout << a(a.size()-1) << std::endl;
	}
	void showvec(Vector3d &a){
		for(int ii=0;ii<2;ii++){
			std::cout << a(ii) << " , " << std::flush;
		}
		std::cout << a(a.size()-1) << std::endl;
	}
	void showvec(Vector3f &a){
		for(int ii=0;ii<2;ii++){
			std::cout << a(ii) << " , " << std::flush;
		}
		std::cout << a(a.size()-1) << std::endl;
	}
}
#endif