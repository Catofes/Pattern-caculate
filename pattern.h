#include "TObject.h"
class input_sample : public TObject//输入参数
{
	public:
	double t_a;
	double t_b;
	double k_1;
	double k_2;
	double k_3;
	double k_4;
	double alpha;
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

		ClassDef(solve_data,1);
};

