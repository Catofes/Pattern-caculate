#include"iostream"
#include"ftttest.h"
#include"TFile.h"
#include"TTree.h"
#include"fftw++/Array.h"
#include"fftw++/fftw++.h"
#include"TCanvas.h"
#include"TGraph2D.h"
#include"TColor.h"
#include"TStyle.h"


using namespace std;
using namespace Array;
using namespace fftwpp;

ClassImp(solve_data);
ClassImp(input_sample);
int ftttest(int imap) {
	TString rootname=Form("%s%d%s","data/",imap,".root");
	TFile *f1=new TFile(rootname);
	TTree *tree=(TTree*)f1->Get("tree");
	solve_data *thedata=0;
	tree->SetBranchAddress("thedata",&thedata);
	tree->GetEntry(0);
	tree->Show();

	fftw::maxthreads=get_max_threads();
	unsigned int nx=20, ny=20;
	size_t align=sizeof(Complex);

	array2<Complex> f(nx,ny,align);

	fft2d Forward2(-1,f);
	fft2d Backward2(1,f);

	for(unsigned int i=0; i < nx; i++) 
	  for(unsigned int j=0; j < ny; j++) 
		f(i,j)=thedata->data[0][20*i+j];

	Forward2.fft(f);
    

	TCanvas *c1=new TCanvas("c1","",300,300);
    gStyle->SetPalette(53,0);
	gStyle->SetNumberContours(256);

	Double_t *a=new Double_t[400];
	Double_t *x=new Double_t[400];
	Double_t *y=new Double_t[400];
	for(int iy=0;iy<20;iy++){
		for(int jx=0;jx<20;jx++){
			a[iy*20+jx]=f(jx,iy).real();
			if((iy+jx)%2==0) a[iy*20+jx]=1;
			else a[iy*20+jx]=0;
			x[iy*20+jx]=jx+0.5;
			y[iy*20+jx]=iy+0.5;
		}
	}   
	TString name=Form("%s%d%s","data/",imap,"ftt.bmp");

	TGraph2D *graphA=new TGraph2D(400,x,y,a);
	graphA->SetName("A");
	graphA->SetNpx(20);
	graphA->SetNpy(20);
	graphA->SetMaximum(1.0001);
	graphA->SetMinimum(-0.0001);
	graphA->SetTitle(name);
	graphA->Draw("COLZ");
	c1->Update();
    c1->Print(name);

}

int main() {
	ftttest(3243);
	return 0;
}
