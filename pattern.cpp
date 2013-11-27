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
#include <TTree.h>
#include <TObject.h>
#include <TFile.h>
#include "pattern.h"
using namespace std;

ClassImp(solve_data);
ClassImp(input_sample);

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


list<input_sample> isample;//输入参数

//list<solve_data> idata;

int load_data()
{

	ifstream inputfile("parameter");
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
		for(int i=0;i<8;i++){
			s.k[i]=atof(strtok(NULL,&c));
		}
		for(int i=0;i<8;i++){
			s.n[i]=atof(strtok(NULL,&c));
		}
		s.alpha=atof(strtok(NULL,&c));
		s.beta=atof(strtok(NULL,&c));
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
double msf(double a){
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
	double abin = *(double *)params;
	double about = *((double *)params+1);
	double bain = *((double *)params+2);
	double baout = *((double *)params+3);
	double aain = *((double *)params+4);
	double aaout = *((double *)params+5);
	double bbin = *((double *)params+6);
	double bbout = *((double *)params+7);
/*
	double fkabin = *((double *)params+8);
	double fkabout = *((double *)params+9);
	double fkbain = *((double *)params+10);
	double fkbaout = *((double *)params+11);
	double fkaain = *((double *)params+12);
	double fkaaout = *((double *)params+13);
	double fkbbin = *((double *)params+14);
	double fkbbout = *((double *)params+15);
*/
	double fnabin = *((double *)params+16);
	double fnabout = *((double *)params+17);
	double fnbain = *((double *)params+18);
	double fnbaout = *((double *)params+19);
	double fnaain = *((double *)params+20);
	double fnaaout = *((double *)params+21);
	double fnbbin = *((double *)params+22);
	double fnbbout = *((double *)params+23);

	double falpha = *((double *)params+24);
	double fbeta = *((double *)params+25);

	double ft_ainv = 1/(*((double *)params+26));
	double ft_binv = 1/(*((double *)params+27));

	double knabin=mypow(*((double *)params+8),*((double *)params+16));
	double knabout=mypow(*((double *)params+9),*((double *)params+17));
	double knbain=mypow(*((double *)params+10),*((double *)params+18));
	double knbaout=mypow(*((double *)params+11),*((double *)params+19));
	double knaain=mypow(*((double *)params+12),*((double *)params+20));
	double knaaout=mypow(*((double *)params+13),*((double *)params+21));
	double knbbin=mypow(*((double *)params+14),*((double *)params+22));
	double knbbout=mypow(*((double *)params+15),*((double *)params+23));


	double A_sum[20][20];
	double B_sum[20][20];
	A_sum[0][0]=(y[2*n]+y[2])/2.0;
	B_sum[0][0]=(y[2*n+1]+y[3])/2.0;
	A_sum[0][n-1]=(y[2*(2*n-1)]+y[2*(n-2)])/2.0;
	B_sum[0][n-1]=(y[2*(2*n-1)+1]+y[2*(n-2)+1])/2.0;
	A_sum[n-1][0]=(y[2*n*(n-2)]+y[2*(n*(n-1)+1)])/2.0;
	B_sum[n-1][0]=(y[2*n*(n-2)+1]+y[2*(n*(n-1)+1)+1])/2.0;
	A_sum[n-1][n-1]=(y[2*(n*(n-2)+n-1)]+y[2*(n*(n-1)+n-2)])/2.0;
	B_sum[n-1][n-1]=(y[2*(n*(n-2)+n-1)+1]+y[2*(n*(n-1)+n-2)+1])/2.0;

	
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

		for (int j = 1; j < n-1; ++j)
		{
			B_sum[i][j]=(y[2*(n*(i-1)+j)+1]+y[2*(n*(i+1)+j)+1]+y[2*(n*i+j-1)+1]+y[2*(n*i+j+1)+1])/4.0;
			A_sum[i][j]=(y[2*(n*(i-1)+j)]+y[2*(n*(i+1)+j)]+y[2*(n*i+j-1)]+y[2*(n*i+j+1)])/4.0;

		}
	}

	///equation below
	for (int i = 0; i < n*n; ++i)
	{

		f[2*i]=ft_ainv*(
			 (bain==0?1.:bain==1?(1.-knbain/(mypow(y[2*i+1],fnbain)+knbain)):(knbain/(mypow(y[2*i+1],fnbain)+knbain)))
			*(baout==0?1.:baout==1?(1.-knbaout/(mypow(fbeta*B_sum[i/20][i%20],fnbaout)+knbaout)):(knbaout/(mypow(fbeta*B_sum[i/20][i%20],fnbaout)+knbaout)))
			*(aain==0?1.:aain==1?(1.-knaain/(mypow(y[2*i],fnaain)+knaain)):(knaain/(mypow(y[2*i],fnaain)+knaain)))
			*(aaout==0?1.:aaout==1?(1.-knaaout/(mypow(falpha*A_sum[i/20][i%20],fnaaout)+knaaout)):(knaaout/(mypow(falpha*A_sum[i/20][i%20],fnaaout)+knaaout)))-y[2*i]);/////-y[2*i]
		f[2*i+1]=ft_binv*(
			 (abin==0?1.:abin==1?(1.-knabin/(mypow(y[2*i],fnabin)+knabin)):(knabin/(mypow(y[2*i],fnabin)+knabin)))
			*(about==0?1.:about==1?(1.-knabout/(mypow(falpha*A_sum[i/20][i%20],fnabout)+knabout)):(knabout/(mypow(falpha*A_sum[i/20][i%20],fnabout)+knabout)))
			*(bbin==0?1.:bbin==1?(1.-knbbin/(mypow(y[2*i+1],fnbbin)+knbbin)):(knbbin/(mypow(y[2*i+1],fnbbin)+knbbin)))
			*(bbout==0?1.:bbout==1?(1.-knbbout/(mypow(fbeta*B_sum[i/20][i%20],fnbbout)+knbbout)):(knbbout/(mypow(fbeta*B_sum[i/20][i%20],fnbbout)+knbbout)))-y[2*i+1]);/////-y[2*i]
	/*	f[2*i+1]=ft_binv*((msf(abin)*mypow(y[2*i],fnabin)+msf(-abin)*mypow(fkabin,fnabin))/(mypow(y[2*i],fnabin)+mypow(fkabin,fnabin))
			*(msf(about)*mypow(falpha*A_sum[i/20][i%20],fnabout)+msf(-about)*mypow(fkabout,fnabout))/(mypow(falpha*A_sum[i/20][i%20],fnabout)+mypow(fkabout,fnabout))
			*(msf(bbin)*mypow(y[2*i+1],fnbbin)+msf(-bbin)*mypow(fkbbin,fnbbin))/(mypow(y[2*i+1],fnbbin)+mypow(fkbbin,fnbbin))
			*(msf(bbout)*mypow(fbeta*B_sum[i/20][i%20],fnbbout)+msf(-bbout)*mypow(fkbbout,fnbbout))/(mypow(fbeta*B_sum[i/20][i%20],fnbbout)+mypow(fkbbout,fnbbout))-y[2*i+1]);/////-y[2*i+1]
*/	}
	return GSL_SUCCESS;
}


int caculate(int inet)//传入网络参数
{
	///文件名///
	char buffer[20];
	sprintf(buffer,"%d",inet);
	string name=buffer;
	const string tname="data/"+name+".txt";
	const string dname="data/"+name+".root";	
	ofstream outputtext (tname.c_str());
	int nn=0;
	//idata.clear();
	TFile hfile(dname.c_str(),"RECREATE","Program Data");
	///网络参数计算///
	int abin=inet%3-1;
	int about=(inet/3)%3-1;
	int bain=(inet/9)%3-1;
	int baout=(inet/27)%3-1;
	int aain=(inet/81)%3-1;
	int aaout=(inet/243)%3-1;
	int bbin=(inet/729)%3-1;
	int bbout=(inet/2187)%3-1;
	outputtext<<"ab    "<<abin<<"	"<<about<<"	"<<"ba    "<<bain<<"	"<<baout<<"aa    "<<aain<<"	"<<aaout<<"	"<<"bb    "<<bbin<<"	"<<bbout<<endl;
	if(about==0&&baout==0&&aaout==0&&bbout==0)//没有反馈的网络
	{
		return 1;
	}
	list<input_sample>::iterator plist;
	list<initial>::iterator ilist;

	double params[28]={0};//定义参数列表
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
	double t = 0.0, t1 = 100000.0;
	double h = 1e-6;

	TTree *t1=new TTree("t1","Program");
	t1 ->Branch("data",&thedata);
	///参数的循环
	for(plist=isample.begin();plist!=isample.end();plist++)
	{
		outputtext<<"t_a	"<<(*plist).t_a<<"	t_b	"<<(*plist).t_b<<endl;
		outputtext<<"kabin  "<<(*plist).k[0]<<"  kabout	"<<(*plist).k[1]<<"	 kbain	"<<(*plist).k[2]<<"	kbaout	"<<(*plist).k[3];
		outputtext<<"	kaain  "<<(*plist).k[4]<<"	 kaaout	"<<(*plist).k[5]<<"	 kbbin	"<<(*plist).k[6]<<"	kbbout	"<<(*plist).k[7];
		outputtext<<"	nabin  "<<(*plist).n[0]<<"  nabout	"<<(*plist).n[1]<<"	 nbain	"<<(*plist).n[2]<<"	nbaout	"<<(*plist).n[3];
		outputtext<<"	naain  "<<(*plist).n[4]<<"	 naaout	"<<(*plist).n[5]<<"	 nbbin	"<<(*plist).n[6]<<"	nbbout	"<<(*plist).n[7];
		outputtext<<"	alpha	"<<(*plist).alpha<<"beta  "<<(*plist).beta<<endl;
		params[0]=(double)abin;
		params[1]=(double)about;
		params[2]=(double)bain;
		params[3]=(double)baout;
		params[4]=(double)aain;
		params[5]=(double)aaout;
		params[6]=(double)bbin;
		params[7]=(double)bbout;
		for(int i=0;i<8;i++){
			params[i+8]=(*plist).k[i];
			params[i+16]=(*plist).n[i];
		}
		params[24]=(*plist).alpha;
		params[25]=(*plist).beta;
		params[26]=(*plist).t_a;
		params[27]=(*plist).t_b;
		///二进制储存params
		thedata.params.alpha=(*plist).alpha;
		thedata.params.beta=(*plist).beta;
		thedata.params.t_a=(*plist).t_a;
		thedata.params.t_b=(*plist).t_b;
		for(int i=0;i<8;i++){
			thedata.params.k[i]=(*plist).k[i];
			thedata.params.n[i]=(*plist).n[i];
		}
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
			//while (t < t1)
			while(1)
			{
				time++;
				int status = gsl_odeiv_evolve_apply (e, c, s,&sys,&t, t1,&h, y);
				if (status != GSL_SUCCESS) break;
				if (time==100){
					delta=0;
					for (int i = 0; i < 2*n*n; ++i)
					{
						delta+=(y[i]-y_1[i])*(y[i]-y_1[i]);
					}
					cout<<"误差"<<delta/(t-told)<<"	"<<t<<endl;
					told=t;
					if (delta<deltamax) break;//break;
					else
					{
						memcpy(y_1, y,sizeof(y));//复制y到y_1;
					}
					time=0;
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
		t1->Fill();
	}
	outputtext.close();
	hfile.Write();
	hfile.Close();
	return 0;
}


int main()
{
	load_data();
	for(int i=3240;i<3321;i++)
	{
		cout<<i<<endl;
		caculate(i);
	}
	return 0;
}

