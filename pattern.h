#include "TObject.h"
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


class solve_data : public TObject
{
	public:
		int inet;
		input_sample params;
		double data[5][800];
		solve_data(){}

		ClassDef(solve_data,1);
};

