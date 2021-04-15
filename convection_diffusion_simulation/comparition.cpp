#include"Variable.h"
#include<iostream>
#include"Eigen/Eigen"
#include<fstream>

using namespace std;
using namespace Eigen;


MatrixXd A = MatrixXd::Zero(n, n);
VectorXd x = VectorXd::Zero(n), b = VectorXd::Zero(n);

void setBoundary(MatrixXd& M, VectorXd& v) {
	//四边west、east、south、north
	for (int i = 0; i <= ny - 1; i++) {
		M(i, i) = 1;
	}
	for (int i = (nx - 1) * ny; i <= nx * ny - 1; i++) {
		M(i, i) = 1;
	}
	for (int i = ny; i < (nx - 1) * ny;) {
		M(i, i) = 1;
		i += ny;
	}
	for (int i = 2 * ny - 1; i < nx * ny - 1;) {
		M(i, i) = 1;
		i += ny;
	}
    //西侧恒定100摄氏度，其他侧0摄氏度
	for (int i = 0; i <= ny - 1; i++) {
		v[i] = 100;
	}
	for (int i = (nx - 1) * ny + 1; i <= nx * ny - 1; i++) {
		v[i] = 0;
	}
	for (int i = ny; i < (nx - 1) * ny;) {
		v[i] = 0;
		i += ny;
	}
	for (int i = 2 * ny ; i <= nx * ny - 1;) {
		v[i] = 0;
		i += ny;
	}
}
void setBoundary1(MatrixXd& M, VectorXd& v) {
	//四边west、east、south、north
	for (int i = 0; i <= ny - 1; i++) {
		M(i, i) = 1;
	}
	for (int i = (nx - 1) * ny; i <= nx * ny - 1; i++) {
		M(i, i) = 1;
	}
	for (int i = ny; i < (nx - 1) * ny;) {
		M(i, i) = 1;
		i += ny;
	}
	for (int i = 2 * ny - 1; i < nx * ny - 1;) {
		M(i, i) = 1;
		i += ny;
	}
	//西侧恒定100摄氏度，其他侧0摄氏度
	for (int i = 0; i <= ny - 1; i++) {
		v[i] = 0;
	}
	for (int i = (nx - 1) * ny + 1; i <= nx * ny - 1; i++) {
		v[i] = 0;
	}
	for (int i = ny; i < (nx - 1) * ny;) {
		v[i] = 0;
		i += ny;
	}
	for (int i = 2 * ny; i <= nx * ny - 1;) {
		v[i] = 0;
		i += ny;
	}
}
void diffusion(MatrixXd& M, VectorXd& v){

	//内部点
	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			int num = i + j;
			M(num, num) = 4 * D + ap0;
			M(num, num - 1) = -D;
			M(num, num + 1) = -D;
			M(num, num - ny) = -D ;
			M(num, num + ny) = -D ;
		}
		i += ny;
	}

	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			v[i + j] = 0;
		}
		i += ny;
	}
}
void convection(MatrixXd& M, VectorXd& v) {
	//内部点
	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			int num = i + j;
			M(num, num) =  F + ap0;                 //test------------------------------
			M(num, num - 1) = 0;
			M(num, num + 1) = 0;
			M(num, num - ny) = 0 - F;
			M(num, num + ny) = 0;
		}
		i += ny;
	}

	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			v[i + j] = 0;
		}
		i += ny;
	}
}
void FUDrenew(MatrixXd& M, VectorXd& v) {

	//内部点
	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			int num = i + j;
			M(num, num) = 4 * D + F + ap0;                 //test------------------------------
			M(num, num - 1) = -D;
			M(num, num + 1) = -D;
			M(num, num - ny) = -D - F ;
			M(num, num + ny) = -D ;
		}
		i += ny;
	}

	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			v[i + j] = 0;
		}
		i += ny;
	}

}
void CDrenew(MatrixXd& M, VectorXd& v) {

	//内部点
	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			int num = i + j;
			M(num, num) = 4 * D + ap0;
			M(num, num - 1) = -D;
			M(num, num + 1) = -D;
			M(num, num - ny) = -D - F / 2;
			M(num, num + ny) = -D + F / 2;
		}
		i += ny;
	}

	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			v[i + j] = 0;
		}
		i += ny;
	}

}
//未完成-----------------------------------------------
void QUICKrenew(MatrixXd& M, VectorXd& v) {

	//内部点
	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			int num = i + j;
			M(num, num) = 4 * D + ap0;
			M(num, num - 1) = -D;
			M(num, num + 1) = -D;
			M(num, num - ny) = -D - F / 2;
			M(num, num + ny) = -D + F / 2;
		}
		i += ny;
	}

	for (int i = ny; i < (nx - 1) * ny; ) {
		for (int j = 1; j <= ny - 2; j++) {
			v[i + j] = 0;
		}
		i += ny;
	}

}

void main() {
	cout << "D/F: " << D / F << endl;

	VectorXd b0 = VectorXd::Zero(n);
	b0[nx * ny / 2 + ny / 2] = 100;
	for (int i = 0; i < 5; i++) {
		diffusion(A, b);
		b += ap0 * b0;
		setBoundary1(A, b);
		x = A.fullPivHouseholderQr().solve(b);
		b0 = x;

	}
	for (int i = 0; i < 100; i++) {
		convection(A, b);
		b += ap0 * b0;
		setBoundary1(A, b);
		x = A.fullPivHouseholderQr().solve(b);
		b0 = x;

		ostringstream name;
		name << i << "data.dat";
		ofstream out(name.str().c_str(), ios::app);
		out << x;
		out << endl;
		out.close();
	}
}