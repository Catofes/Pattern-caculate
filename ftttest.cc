#include"iostream"
#include"ftttest.h"
#include"TFile.h"
#include"TTree.h"
using namespace std;
ClassImp(solve_data);
ClassImp(input_sample);
int main(){
	TFile *f1=new TFile("data/3243.root");
	TTree *t1=(TTree*)f1->Get("t1");
	solve_data thedata;
	t1->SetBranchAddress("data",&thedata);
	t1->GetEntry(0);
	cout<<thedata.data[0][1];
	
//data
//	t1->Show(0);
}
