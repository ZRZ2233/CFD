#include<algorithm>
//#include"function.h"
#include<iostream>

using namespace std;

//二维五点矩阵ILU分解，并将因子应用于更新方程的一次迭代。
//a,b,c,d,e代表5对角矩阵从下到上五条对角线
void ilu5(int m, int n, double* a, double* b, double* c, double* d, double* e, double* r, double* x)
{

    //a[0] = a[0];
    //b[0] = b[0];
    //c[0] = c[0];
    d[0] = d[0] / c[0];
    e[0] = e[0] / c[0];
    x[0] = r[0] / c[0];

    for (int i = 1; i < m - 1; i++)
    {
        //a[i] = a[i];
        //b[i] = b[i];
        c[i] = c[i] - b[i] * d[i - 1];
        d[i] = d[i] / c[i];
        e[i] = e[i] / c[i];
        x[i] = (r[i] - b[i] * x[i - 1]) / c[i];
    }
    for (int i = m; i < m * n - 1; i++)
    {
        //a[i] = a[i];
        //b[i] = b[i];
        c[i] = c[i] - a[i] * e[i - m] - b[i] * d[i - 1];
        d[i] = d[i] / c[i];
        e[i] = e[i] / c[i];
        x[i] = (r[i] - b[i] * x[i - 1] - a[i] * x[i - m]) / c[i];
    }

    x[m * n - 1] = x[m * n - 1];

    for (int i = m * n - 2; i > m * (n - 1); i--)
    {
        x[i] = x[i] - d[i] * x[i + 1];
    }
    for (int i = m * (n - 1); i > 0; i--)
    {
        x[i] = x[i] - d[i] * x[i + 1] - e[i] * x[i + m];
    }
}

//二维五点矩阵SIP-ILU分解，并将因子应用于更新方程的一次迭代。
void sip5(int m, int n, double alf, double* a, double* b, double* c, double* d, double* e, double* r, double* x)
{
    if (alf == 0 || (m <= 3 || n <= 3))
        ilu5(m, n, a, b, c, d, e, r, x);
    else
    {
        a[0] = a[0];
        b[0] = b[0];
        c[0] = c[0];
        d[0] = d[0] / c[0];
        e[0] = e[0] / c[0];
        x[0] = r[0] / c[0];

        for (int i = 1; i <= m - 1; i++)
        {
            a[i] = a[i];
            b[i] = b[i] / (1 + alf * e[i - 1]);
            c[i] = c[i] - b[i] * d[i - 1] + alf * b[i] * e[i - 1];
            d[i] = d[i] / c[i];
            e[i] = (e[i] - alf * b[i] * e[i - 1]) / e[i - 1];
            x[i] = (r[i] - b[i] * x[i - 1]) / c[i];
        }
        for (int i = m; i <= m * n - 1; i++)
        {
            a[i] = a[i] / (1 + alf * d[i - m]);
            b[i] = b[i] / (1 + alf * e[i - m]);
            c[i] = c[i] - a[i] * e[i - m] - b[i] * d[i - 1] + alf * (a[i] * d[i - m] + b[i] * e[i - 1]);
            d[i] = (d[i] - alf * a[i] * d[i - m]) / c[i];
            e[i] = (e[i] - alf * b[i] * e[i - 1]) / c[i];
            x[i] = (r[i] - b[i] * x[i - 1] - a[i] * x[i - m]) / c[i];
        }


		x[m * n - 1] = x[m * n - 1];

		for (int i = m * n - 2; i >= m * (n - 1); i--)
		{
			x[i] = x[i] - d[i] * x[i + 1];
		}
		for (int i = m * (n - 1)-1; i >= 0; i--)
		{
			x[i] = x[i] - d[i] * x[i + 1] - e[i] * x[i + m];
		}

              
    }
}

//调整五点矩阵的对角线，以满足与对角线优势相关的属性
void adjustM(int m, int n, double* matrix, double rf)
{
    int a = n * 5;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            matrix[i * a + j * 5 + 2] = max(matrix[i * a + j * 5 + 2], -rf * (matrix[i * a + j * 5 + 0] + matrix[i * a + j * 5 + 1] + matrix[i * a + j * 5 + 3] + matrix[i * a + j * 5 + 4]));
        }
    }
}

//m,n注意----------------------来自网上-----误差在2%左右----------
int SIPsolver(double* ALL, double* AL, double* AC, double* AR, double* ARR, double* x, double* r, int nn, int mm)
{
	//int m, int n
	int n = nn * mm - 1;
	int m = mm-1;
	//求解如下形式矩阵：
	//AC0	AR0			ARR0						
	//AL1	AC1	AR1			ARR1					
	//	    AL2	AC2	AR2			ARR2				
	//		    AL3	AC3	AR3			ARR3			
	//ALLm			AL4	AC4	AR4			ARR4		
	//	ALL5			AL5	AC5	AR5			ARRn-m	
	//		ALL6			AL6	AC6	AR6			
	//			ALL7			AL7	AC7	AR7		
	//				ALL8			AL8	AC8	AR8	
	//					ALLn			ALn	ACn	
	double a = 0.99;
	double P1 = 0;
	double P2 = 0;
	double* LLL, * LL, * LC;
	double* UR, * URR;
	//double* RES, * y;      
	double* RES;                //delta
	double RESN = 0;            //收敛判据
	int MAXITER = 1000;         //最大循环次数
	int ITER;                  //计数器
	LLL = new double[n + 1];
	LL = new double[n + 1];
	LC = new double[n + 1];
	UR = new double[n + 1];
	URR = new double[n + 1];
	RES = new double[n + 1];
	for (int i = 0; i <= n; i++) {
		LLL[i] = 0;
		LL[i] = 0;
		LC[i] = 0;
		UR[i] = 0;
		URR[i] = 0;
		RES[i] = 0;
	}
	for (int i = 0; i <= n; i++)//相比solps内的没有对x的改变了
	{
		if (i >= m) {
			LLL[i] = ALL[i] / (1 + a * UR[i - m]);
		}
		if (i >= 1) {
			LL[i] = AL[i] / (1 + a * URR[i - 1]);
		}
		if (i >= m) {
			P1 = a * LLL[i] * UR[i - m];
		}
		else {
			P1 = 0;
		}
		if (i >= 1) {
			P2 = a * LL[i] * URR[i - 1];
		}
		else {
			P2 = 0;
		}
		if (i < 1) {
			LC[i] = AC[i] + (P1 + P2);
		}
		else if (i >= 1 && i < m) {
			LC[i] = AC[i] + (P1 + P2) - LL[i] * UR[i - 1];
		}
		else {
			LC[i] = AC[i] + (P1 + P2) - LLL[i] * URR[i - m] - LL[i] * UR[i - 1];
		}
		if (i <= n - 1) {
			UR[i] = (AR[i] - P1) / LC[i];
		}
		else {
			UR[i] = 0;
		}
		if (i <= n - m) {
			URR[i] = (ARR[i] - P2) / LC[i];
		}
		else {
			URR[i] = 0;
		}

	}
	for (int i = 0; i <= n ; i++) {
		cout << LLL[i] << "---" << LL[i] << "---" << LC[i] << "---" << UR[i] << "---" << URR[i] << endl;
	}
	
	ITER = 0;
	do {
		RESN = 0;
		for (int i = 0; i <= n; i++)
		{
			RES[i] = r[i] - (ALL[i] * (i - m >= 0 ? x[i - m] : 0) + AL[i] * (i - 1 >= 0 ? x[i - 1] : 0) + AC[i] * x[i] + AR[i] * (i + 1 <= n ? x[i + 1] : 0) + ARR[i] * (i + m <= n ? x[i + m] : 0));
			RESN += abs(RES[i]);
		}
		


		for (int i = 0; i <= n; i++)
		{
			RES[i] = (RES[i] - (i - 1 >= 0 ? RES[i - 1] : 0) * LL[i] - (i - m >= 0 ? RES[i - m] : 0) * LLL[i]) / LC[i];
		}

		for (int i = n; i >= 0; i--)
		{
			x[i] = x[i] + (RES[i] - (i + 1 <= n ? RES[i + 1] : 0) * UR[i] - (i + m <= n ? RES[i + m] : 0) * URR[i]);
		}
		ITER++;
		//bool bb = abs(RESN) > 1E-6 && ITER < MAXITER;
		//cout << "num: "<<ITER <<" RESN: " << abs(RESN) << "---" << bb << endl;
		//for (int i = 0; i < 16; i++) {
		//	//cout <<a[i]<<"----------"<<b[i]<<"-------"<<c[i]<<"------------"<<d[i]<<"----------"<< e[i] << endl;
		//	cout << x[i] << endl;
		//}
	} while (abs(RESN) > 1E-6 && ITER < MAXITER);

	delete[] LLL;
	delete[] LL;
	delete[] LC;
	delete[] UR;
	delete[] URR;
	delete[] RES;
	return ITER > MAXITER ? -1 : 1;
}

//结合网上那个和solps内置的，和Eigen算的结果一样
//void sip(int m, int n, double alf, double* a, double* b, double* c, double* d, double* e, double* r, double* x)
//{
//	double p1 = 0;
//	double p2 = 0;
//	double* acopy = new double[n * m];
//	double* bcopy = new double[n * m];
//	double* ccopy = new double[n * m];
//	double* dcopy = new double[n * m];
//	double* ecopy = new double[n * m];
//	double* res = new double[n * m];
//	int count = 0;      //循环计数
//	double con_jud = 0; //收敛判断
//	for (int i = 0; i < n * m; i++)
//	{
//		acopy[i] = a[i];
//		bcopy[i] = b[i];
//		ccopy[i] = c[i];
//		dcopy[i] = d[i];
//		ecopy[i] = e[i];
//		res[i] = 0;
//	}
//
//	acopy[0] = acopy[0];
//	bcopy[0] = bcopy[0];
//	ccopy[0] = ccopy[0];
//	dcopy[0] = dcopy[0] / ccopy[0];
//	ecopy[0] = ecopy[0] / ccopy[0];
//
//	for (int i = 1; i <= m - 1; i++)
//	{
//		acopy[i] = acopy[i];
//		bcopy[i] = bcopy[i] / (1 + alf * ecopy[i - 1]);
//		ccopy[i] = ccopy[i] - bcopy[i] * dcopy[i - 1] + alf * bcopy[i] * ecopy[i - 1];
//		dcopy[i] = dcopy[i] / ccopy[i];
//		ecopy[i] = (ecopy[i] - alf * bcopy[i] * ecopy[i - 1]) / ecopy[i - 1];
//	}
//	for (int i = m; i <= m * n - 1; i++)
//	{
//		acopy[i] = acopy[i] / (1 + alf * dcopy[i - m]);
//		bcopy[i] = bcopy[i] / (1 + alf * ecopy[i - m]);
//		ccopy[i] = ccopy[i] - acopy[i] * ecopy[i - m] - bcopy[i] * dcopy[i - 1] + alf * (acopy[i] * d[i - m] + bcopy[i] * e[i - 1]);
//		dcopy[i] = (dcopy[i] - alf * acopy[i] * d[i - m]) / ccopy[i];
//		ecopy[i] = (ecopy[i] - alf * bcopy[i] * e[i - 1]) / ccopy[i];
//		/*ccopy[i] = ccopy[i] - acopy[i] * ecopy[i - m] - bcopy[i] * dcopy[i - 1] + alf * (acopy[i] * dcopy[i - m] + bcopy[i] * ecopy[i - 1]);
//		dcopy[i] = (dcopy[i] - alf * acopy[i] * dcopy[i - m]) / ccopy[i];
//		ecopy[i] = (ecopy[i] - alf * bcopy[i] * ecopy[i - 1]) / ccopy[i];*/
//	}
//	
//	do
//	{
//		con_jud = 0;
//		for (int i = 0; i < n * m; i++) //收敛判断
//		{
//			res[i] = r[i] - (a[i] * (i - m >= 0 ? x[i - m] : 0) + b[i] * (i - 1 >= 0 ? x[i - 1] : 0) + c[i] * x[i] + d[i] * (i + 1 <= n * m - 1 ? x[i + 1] : 0) + e[i] * (i + m <= n * m - 1 ? x[i + m] : 0));
//			con_jud += abs(res[i]);
//		}
//
//		for (int i = 0; i < n * m; i++)  //更新差值
//		{ 
//			res[i] = (res[i] - (i - 1 >= 0 ? res[i - 1] : 0) * bcopy[i] - (i - m >= 0 ? res[i - m] : 0) * acopy[i]) / ccopy[i];
//		}
//
//		for (int i = n * m - 1; i >= 0; i--)  //根据差值，更新解
//		{
//			x[i] = x[i] + (res[i] - (i + 1 <= n * m - 1 ? res[i + 1] : 0) * dcopy[i] - (i + m <= n * m - 1 ? res[i + m] : 0) * ecopy[i]);
//		}
//	} while (abs(con_jud) > 1E-6 && count < 1000);
//}

//根据上一个进一步更改，复制改为赋0，加快速度
void sip(int m, int n, double alf, double* a, double* b, double* c, double* d, double* e, double* r, double* x)
{
	double p1 = 0;
	double p2 = 0;
	double* acopy = new double[n * m];
	double* bcopy = new double[n * m];
	double* ccopy = new double[n * m];
	double* dcopy = new double[n * m];
	double* ecopy = new double[n * m];
	double* res = new double[n * m];
	int count = 0;      //循环计数
	double con_jud = 0; //收敛判断
	for (int i = 0; i < n * m; i++)
	{
		acopy[i] = 0;
		bcopy[i] = 0;
		ccopy[i] = 0;
		dcopy[i] = 0;
		ecopy[i] = 0;
		res[i] = 0;
	}

	acopy[0] = a[0];
	bcopy[0] = b[0];
	ccopy[0] = c[0];
	dcopy[0] = d[0] / ccopy[0];
	ecopy[0] = e[0] / ccopy[0];

	for (int i = 1; i <= m - 1; i++)
	{
		acopy[i] = a[i];
		bcopy[i] = b[i] / (1 + alf * ecopy[i - 1]);
		ccopy[i] = c[i] - bcopy[i] * dcopy[i - 1] + alf * bcopy[i] * ecopy[i - 1];
		dcopy[i] = d[i] / ccopy[i];
		ecopy[i] = (e[i] - alf * bcopy[i] * ecopy[i - 1]) / ccopy[i - 1];
	}
	for (int i = m; i <= m * n - 1; i++)
	{
		acopy[i] = a[i] / (1 + alf * dcopy[i - m]);
		bcopy[i] = b[i] / (1 + alf * ecopy[i - m]);
		ccopy[i] = c[i] - acopy[i] * ecopy[i - m] - bcopy[i] * dcopy[i - 1] + alf * (acopy[i] * d[i - m] + bcopy[i] * e[i - 1]);
		dcopy[i] = (d[i] - alf * acopy[i] * d[i - m]) / ccopy[i];
		ecopy[i] = (e[i] - alf * bcopy[i] * e[i - 1]) / ccopy[i];
	}

	do
	{
		con_jud = 0;
		for (int i = 0; i < n * m; i++) //收敛判断
		{
			res[i] = r[i] - (a[i] * (i - m >= 0 ? x[i - m] : 0) + b[i] * (i - 1 >= 0 ? x[i - 1] : 0) + c[i] * x[i] + d[i] * (i + 1 <= n * m - 1 ? x[i + 1] : 0) + e[i] * (i + m <= n * m - 1 ? x[i + m] : 0));
			con_jud += abs(res[i]);
		}

		for (int i = 0; i < n * m; i++)  //更新差值
		{
			res[i] = (res[i] - (i - 1 >= 0 ? res[i - 1] : 0) * bcopy[i] - (i - m >= 0 ? res[i - m] : 0) * acopy[i]) / ccopy[i];
		}

		for (int i = n * m - 1; i >= 0; i--)  //根据差值，更新解
		{
			x[i] = x[i] + (res[i] - (i + 1 <= n * m - 1 ? res[i + 1] : 0) * dcopy[i] - (i + m <= n * m - 1 ? res[i + m] : 0) * ecopy[i]);
		}
	} while (abs(con_jud) > 1E-6 && count < 1000);
}
