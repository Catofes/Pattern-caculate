#include "TObject.h"
class input_sample : public TObject//输入参数
{
	public:
	double t_a;
	double t_b;
	double k[8];
	double n[8];
	double alpha,beta;
	input_sample(){}
	ClassDef(input_sample,1);
};


class solve_data : public TObject
{
	public:
		int inet;
		input_sample params;
		double data[5][800];
		solve_data(){}

		ClassDef(solve_data,2);
};

