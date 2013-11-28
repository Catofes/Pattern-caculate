#include"iostream"
#include"ftttest.h"
#include"TFile.h"
#include"TTree.h"
using namespace std;
ClassImp(solve_data);
ClassImp(input_sample);
int main(){
	TFile *f1=new TFile("data/3243.root");
	TTree *ptree=(TTree*)f1->Get("t1");
	solve_data *thedata;
	cout<<&thedata<<endl;
	getchar();
	ptree->SetBranchAddress("data",&thedata);
	cout<<&thedata<<endl;
	getchar();
	ptree->GetEntry(1);
	cout<<&thedata<<endl;
	getchar();
	cout<<thedata->inet<<"  "<<thedata->params.t_a<<"  "<<thedata->data[2][3];
	
//data
//	t1->Show(0);
}
