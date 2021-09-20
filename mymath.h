#ifndef MYMATH_H
#define MYMATH_H

#define PRINT_MAT(X) std::cout << #X << ":\n" << X << std::endl << std::endl

using namespace Eigen;

namespace Mymath{
	class Funcvec{
		public:
			Funcvec();
			virtual VectorXd function(VectorXd x);
			virtual ~Funcvec();
	};
	
	double sign(double A);
	double step(double A);
	double sigmoid(double x,double alp);
	VectorXd sigmoid(VectorXd x,double alp);
	/*
	void eig(MatrixXd aA,MatrixXd bB,MatrixXd &eigenV,MatrixXd &eigenD);
	void eig(MatrixXd aA,MatrixXd &eigenV,MatrixXd &eigenD);
	*/
	MatrixXd inv(MatrixXd aA);
	MatrixXd pseudo_inv(MatrixXd aA);
	MatrixXd absmat(MatrixXd aA);
	MatrixXd diffvec(VectorXd x,Funcvec *func);
	MatrixXd diffdiffvec(VectorXd x,Funcvec *func);
	MatrixXd diffnorm(VectorXd x,Funcvec *func);
	void showvec(VectorXd &a);
	void showvec(VectorXf &a);
	void showvec(Vector3d &a);
	void showvec(Vector3f &a);
}

#endif