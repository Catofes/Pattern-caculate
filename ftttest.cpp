#include"iostream"
#include"ftttest.h"
#include"TFile.h"
#include"TTree.h"
using namespace std;
ClassImp(solve_data);
ClassImp(input_sample);
int ftttest() {
	TFile *f1=new TFile("data/3243.root");
	TTree *tree=(TTree*)f1->Get("tree");
	solve_data *thedata=0;
	tree->SetBranchAddress("thedata",&thedata);
	tree->GetEntry(0);
	tree->Show();
	//solve_data * thedata=0;
	//tree->SetBranchAddress("thedata",&thedata);
	//TBranch *branch  = t1->GetBranch("data");
	//branch->SetAddress(&thedata);
	//tree->GetEntry(0);
	//tree->Print(0);
	//TBranch *b_destep = t2->GetBranch("inet");
	//cout<<thedata<<"	";
	//tree->StartViewer();
	//cout<<thedata->inet<<endl;
	//for (i = 0; i < 20; i++) {
	//	for (j = 0; j < 20; j++) {
	//		cout<<thedata->data[0][20*i+j]<<"	";
	//	}
	//	cout<<endl;
	//}
	return 0;
}

int main() {
	ftttest();
	return 0;
}
