// pattern0.cpp : 定义控制台应用程序的入口点。
//
//#define _CRT_SECURE_NO_WARNINGS
//#define _CRT_SECURE_NO_DEPRECATE

//#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <string>
#include <string.h> 
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
using namespace std;

static int n1=2;
static int n2=2;
static int n3=2;
static int n4=2;
static int n=20;
static int T=200;
static double sigma=0.02;


static int randomtime=5;
static int nnet=81;
static int nset=100;


static double tmin=0.1;
static double tmax=10.0;
static double kmin=0.01;
static double kmax=1.0;
static double alphamin=0.0;
static double alphamax=2.0;
static double ndim=6;
static double deltamax=0.0000001;

const char c='\t';

struct initial
{
	double data[800];
};

list<initial> iinitial;

struct input_sample//输入参数
{
	double t_a;
	double t_b;
	double k_1;
	double k_2;
	double k_3;
	double k_4;
	double alpha;
};

list<input_sample> isample;//输入参数

struct solve_data
{
	int inet;
	input_sample params;
	double data[5][800];
};

list<solve_data> idata;

int load_data()
{

	ifstream inputfile("data");
	string linshi;
	//cout<<"OK0";
	input_sample s;
	char *kk=new char[500];
	//cout<<"OK1";
	while (getline(inputfile,linshi))
	{
		if(linshi.length()<5)continue;
		strcpy(kk,linshi.c_str());
		s.t_a=atof(strtok(kk,&c));
		s.t_b=atof(strtok(NULL,&c));
		s.k_1=atof(strtok(NULL,&c));
		s.k_2=atof(strtok(NULL,&c));
		s.k_3=atof(strtok(NULL,&c));
		s.k_4=atof(strtok(NULL,&c));
		s.alpha=atof(strtok(kk,&c));
		isample.push_back(s);
	}
	//cout<<"OK2";
	inputfile.close();
	inputfile.open("initials");
	linshi.clear();
	initial k;
	char *pp=new char[50000];
	cout<<"loaddata"<<endl;
	while (getline(inputfile,linshi))
	{
		if(linshi.length()<5)continue;
		strcpy(pp,linshi.c_str());
		for (int i = 0; i < 2*n*n; ++i)
		{
			//cout<<i<<"	";
			if(i==0)k.data[i]=atof(strtok(pp,&c));
			else k.data[i]=atof(strtok(NULL,&c));
		}
		iinitial.push_back(k);
	}
	cout<<"OK"<<endl;
	delete kk;
	delete pp;
	return 0;
}
double mypow(double a,double b)
{
	if(a<0)a=0.000000001;
	if(a>1)a=1.;
	return pow(a,b);
}
double mystepfunction(double a){
	return a>=0.?1.:0.;
}

int func (double t, const double y[], double f[],void *params)
{	
	///////////////
	/*说明
		paramas是传入的所有参数
		y是自变量，共有2*n*n个800个，其中y[2*(x*n+y)]代表A[x][y](从0-19共20行20列)
		y[2*(x*n+y)+1]代表B[x][y]。
		下面过程中的A_sum矩阵是指B对A的影响。同理B_sum
		f是方程，共有800个每个点一个方程，排序和y保持一致。
	*/
	///////////////
	double abin = *(double *)params;
	double about = *((double *)params+1);
	double bain = *((double *)params+2);
	double baout = *((double *)params+3);
	double fk1 = *((double *)params+4);
	double fk2 = *((double *)params+5);
	double fk3 = *((double *)params+6);
	double fk4 = *((double *)params+7);
	double falpha = *((double *)params+8);

	double A_sum[20][20];
	double B_sum[20][20];
	///changed!lxr:A->2n,B->2n+1
	A_sum[0][0]=(y[2*n]+y[2])/2.0;
	B_sum[0][0]=(y[2*n+1]+y[3])/2.0;
	A_sum[0][n-1]=(y[2*(2*n-1)]+y[2*(n-2)])/2.0;
	B_sum[0][n-1]=(y[2*(2*n-1)+1]+y[2*(n-2)+1])/2.0;
	A_sum[n-1][0]=(y[2*n*(n-2)]+y[2*(n*(n-1)+1)])/2.0;
	B_sum[n-1][0]=(y[2*n*(n-2)+1]+y[2*(n*(n-1)+1)+1])/2.0;
	A_sum[n-1][n-1]=(y[2*(n*(n-2)+n-1)]+y[2*(n*(n-1)+n-2)])/2.0;
	B_sum[n-1][n-1]=(y[2*(n*(n-2)+n-1)+1]+y[2*(n*(n-1)+n-2)+1])/2.0;//////////

	/*A_sum[0][0]=(y[2*n+1]+y[3])/2.0;
	B_sum[0][0]=(y[2*n]+y[2])/2.0;
	A_sum[0][n-1]=(y[2*n+1]*y[2*(n-2)+1])/2.0;
	B_sum[0][n-1]=(y[2*n]*y[2*(n-2)])/2.0;
	A_sum[n-1][0]=(y[2*n*(n-2)+1]+y[2*(n*(n-1)+1)+1])/2.0;
	B_sum[n-1][0]=(y[2*n*(n-2)]+y[2*(n*(n-1)+1)])/2.0;
	A_sum[n-1][n-1]=(y[2*(n*(n-2)+n-1)+1]+y[2*(n*(n-1)+n-2)])/2.0*/
	
	for (int i = 1; i < n-1; ++i)
	{
		B_sum[i][0]=(y[2*n*(i-1)+1]+y[2*(n*i+1)+1]+y[2*n*(i+1)+1])/3.0;
		A_sum[i][0]=(y[2*n*(i-1)]+y[2*(n*i+1)]+y[2*n*(i+1)])/3.0;
		B_sum[0][i]=(y[2*(i-1)+1]+y[2*(n+i)+1]+y[2*(i+1)+1])/3.0;
		A_sum[0][i]=(y[2*(i-1)]+y[2*(n+i)]+y[2*(i+1)])/3.0;
		B_sum[i][n-1]=(y[2*(n*(i-1)+n-1)+1]+y[2*(n*i+n-2)+1]+y[2*(n*(i+1)+n-1)+1])/3.0;
		A_sum[i][n-1]=(y[2*(n*(i-1)+n-1)]+y[2*(n*i+n-2)]+y[2*(n*(i+1)+n-1)])/3.0;
		B_sum[n-1][i]=(y[2*(n*(n-1)+i-1)+1]+y[2*(n*(n-2)+i)+1]+y[2*(n*(n-1)+i+1)+1])/3.0;
		A_sum[n-1][i]=(y[2*(n*(n-1)+i-1)]+y[2*(n*(n-2)+i)]+y[2*(n*(n-1)+i+1)])/3.0;

		//A_sum[i][0]=(y[2*n*(i-1)+1]+y[2*(n*i+1)+1]+y[2*n*(i+1)+1])/3.0;
		//B_sum[i][0]=(y[2*n*(i-1)]+y[2*(n*i+1)]+y[2*n*(i+1)])/3.0;
		//A_sum[0][i]=(y[2*(i-1)+1]+y[2*(n+i)+1]+y[2*(i+1)+1])/3.0;
		//B_sum[0][i]=(y[2*(i-1)]+y[2*(n+i)]+y[2*(i+1)])/3.0;

		for (int j = 1; j < n-1; ++j)
		{
			B_sum[i][j]=(y[2*(n*(i-1)+j)+1]+y[2*(n*(i+1)+j)+1]+y[2*(n*i+j-1)+1]+y[2*(n*i+j+1)+1])/4.0;
			A_sum[i][j]=(y[2*(n*(i-1)+j)]+y[2*(n*(i+1)+j)]+y[2*(n*i+j-1)]+y[2*(n*i+j+1)])/4.0;
			//A_sum[i][j]=(y[2*(n*(i-1)+j)+1]+y[2*(n*(i+1)+j)+1]+y[2*(n*i+j-1)+1]+y[2*(n*i+j+1)+1])/4.0;
			//B_sum[i][j]=(y[2*(n*(i-1)+j)]+y[2*(n*(i+1)+j)]+y[2*(n*i+j-1)]+y[2*(n*i+j+1)])/4.0;
		}
	}

	///以下是方程
	for (int i = 0; i < n*n; ++i)
	{
		f[2*i]=(mystepfunction(bain)*mypow(y[2*i+1],n1)+mystepfunction(-bain)*mypow(fk1,n1))/(mypow(y[2*i+1],n1)+mypow(fk1,n1))
			*(mystepfunction(baout)*mypow(falpha*B_sum[i/20][i%20],n3)+mystepfunction(-baout)*mypow(fk3,n3))/(mypow(falpha*B_sum[i/20][i%20],n3)+mypow(fk3,n3))-y[2*i];/////-y[2*i]
		f[2*i+1]=(mystepfunction(abin)*mypow(y[2*i],n2)+mystepfunction(-abin)*mypow(fk2,n2))/(mypow(y[2*i],n2)+mypow(fk2,n2))
			*(mystepfunction(about)*mypow(falpha*A_sum[i/20][i%20],n4)+mystepfunction(-about)*mypow(fk4,n4))/(mypow(falpha*A_sum[i/20][i%20],n4)+mypow(fk4,n4))-y[2*i+1];/////-y[2*i+1]
		//f[2*i]=(mystepfunction(bain)*mypow(y[2*i+1],n1)+mystepfunction(-bain)*mypow(fk1,n1))/(mypow(y[2*i+1],n1)+mypow(fk1,n1))
		//	*(mystepfunction(baout)*falpha*mypow(A_sum[i/20][i%20],n3)+mystepfunction(-baout)*mypow(fk3,n3))/(mypow(A_sum[i/20][i%20],n3)+mypow(fk3,n3))-y[2*i];/////-y[2*i]
		//f[2*i+1]=(mystepfunction(abin)*mypow(y[2*i],n2)+mystepfunction(-abin)*mypow(fk2,n2))/(mypow(y[2*i],n2)+mypow(fk2,n2))
		//	*(mystepfunction(about)*falpha*mypow(B_sum[i/20][i%20],n4)+mystepfunction(-about)*mypow(fk4,n4))/(mypow(B_sum[i/20][i%20],n4)+mypow(fk4,n4))-y[2*i+1];/////-y[2*i+1]
	}
	return GSL_SUCCESS;
}


int caculate(int inet)//传入网络参数
{
	///文件名///
	char buffer[20];
	sprintf(buffer,"%d",inet);
	string name=buffer;
	const string tname=name+".txt";
	const string dname=name+".dat";	
	ofstream outputtext (tname.c_str());
	int nn=0;
	idata.clear();
	///网络参数计算///
	int abin=inet%3-1;
	int about=(inet/3)%3-1;
	int bain=(inet/27)%3-1;
	int baout=(inet/81)%3-1;
	outputtext<<"ab    "<<abin<<"	"<<about<<"	"<<"ba    "<<bain<<"	"<<baout<<endl;
	if(about==0&&baout==0)//没有反馈的网络
	{
		return 1;
	}
	list<input_sample>::iterator plist;
	list<initial>::iterator ilist;

	double params[9]={0};//定义参数列表
	double y[800]={0};//解的值；
	double y_1[800]={0};//前面一步的解得值；
	double delta=0;
	solve_data thedata;
	thedata.inet=inet;
	///一下gsl ode解法定义
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 800);
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-6, 0.0);
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (800);
	gsl_odeiv_system sys = {func, NULL, 800, params};
	double t = 0.0, t1 = 100.0;
	double h = 1e-6;

	///参数的循环
	for(plist=isample.begin();plist!=isample.end();plist++)
	{
		outputtext<<"t_a	"<<(*plist).t_a<<"	t_b	"<<(*plist).t_b<<"	k_1	"<<(*plist).k_1<<"	k_2	"<<(*plist).k_2<<"	k_3	"<<(*plist).k_3<<"	k_4	"<<(*plist).k_4<<"	alpha	"<<(*plist).alpha<<endl;
		params[0]=(double)abin;
		params[1]=(double)about;
		params[2]=(double)bain;
		params[3]=(double)baout;
		params[4]=(*plist).k_1;
		params[5]=(*plist).k_2;
		params[6]=(*plist).k_3;
		params[7]=(*plist).k_4;
		params[8]=(*plist).alpha;
		///二进制储存params
		thedata.params.alpha=(*plist).alpha;
		thedata.params.k_1=(*plist).k_1;
		thedata.params.k_2=(*plist).k_2;
		thedata.params.k_3=(*plist).k_3;
		thedata.params.k_4=(*plist).k_4;
		thedata.params.t_a=(*plist).t_a;
		thedata.params.t_b=(*plist).t_b;
		///初始值的循环
		nn=0;
		cout<<"caculate"<<endl;
		for(ilist=iinitial.begin();ilist!=iinitial.end();ilist++)
		{
			///载入初值
			for (int i = 0; i < 2*n*n; ++i)
			{
				y[i]=(*ilist).data[i];
				y_1[i]=(*ilist).data[i];
			}
			double delta=0;
			int time=0;
			double told=0;
			while (t < t1)
			{
				time++;
				int status = gsl_odeiv_evolve_apply (e, c, s,&sys,&t, t1,&h, y);
				if (status != GSL_SUCCESS) break;
				delta =0;
				for (int i = 0; i < 2*n*n; ++i)
				{
					delta+=(y[i]-y_1[i])*(y[i]-y_1[i]);
				}
				//if(time==100){
				//	cout<<"误差"<<delta/(t-told)<<"	"<<t<<endl;
				//	time=0;
				//}
				//cout<<y[0]<<"  "<<y[1]<<endl;
				//getchar();
				told=t;
				if (delta<deltamax) ;//break;
				else
				{
					memcpy(y_1, y,sizeof(y));//复制y到y_1;
				}
			}
			outputtext<<"初始值:	"<<nn<<endl;
			outputtext<<"A:"<<endl;
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					outputtext<<y[2*(i*n+j)]<<"	";
				}
				outputtext<<endl;
			}
			outputtext<<"B:"<<endl;
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					outputtext<<y[2*(i*n+j)+1]<<"	";
				}
				outputtext<<endl;
			}
			for (int i = 0; i < 2*n*n; ++i)
			{
				thedata.data[nn][i]=y[i];
			}
			nn++;
		}
		idata.push_back(thedata);
	}
	outputtext.close();
	return 0;
}


int main()
{
	load_data();
	for(int i=0;i<81;i++)
	{
		cout<<i<<endl;
		caculate(i);
	}
	//caculate(2);
//	getchar();
//	cout<<"OK";
	getchar();
	return 0;
}

