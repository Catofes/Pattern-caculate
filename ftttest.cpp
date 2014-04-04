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
#include"TH2D.h"

using namespace std;
using namespace Array;
using namespace fftwpp;

ClassImp(solve_data);
ClassImp(input_sample);
void draw(Double_t *data,Double_t *datain,TString name,TCanvas *c1,TH2D *h,TH2D *hin,Int_t x=20, Int_t y=20){
	c1->Divide(2,1);
	c1->cd(1);
	h->SetStats(kFALSE);
	h->SetMinimum(-0.0001);
	for(Int_t i=0;i<x*y;i++){
		for(Int_t jhight=0;jhight<data[i]*100;jhight++){
			h->Fill(i%x+0.5,i/y+0.5);
		}
	}
	h->Draw("COLZ");
	c1->cd(2);
	hin->SetStats(kFALSE);
	hin->SetMinimum(-0.0001);
	for(Int_t i=0;i<x*y;i++){
		for(Int_t jhight=0;jhight<datain[i]*100;jhight++){
			hin->Fill(i%x+0.5,i/y+0.5);
		}
	}
	hin->Draw("COLZ");

	c1->Print(name);
}


int ftttest(int size,int imap,TCanvas *c1,TH2D *h1,TH2D *h2) {
	TString rootname=Form("%s%d%s","data/",imap,".root");
	TFile *f1=new TFile(rootname);
	TTree *tree=(TTree*)f1->Get("tree");
	solve_data *thedata=0;
	tree->SetBranchAddress("thedata",&thedata);

	int nn=tree->GetEntries();
	for ( int ni=0;ni<nn;ni++){
		tree->GetEntry(ni);
		tree->Show();

		//	fftw::maxthreads=get_max_threads();
		unsigned int nx=size,ny=size;//nx=20, ny=20;
		size_t align=sizeof(Complex);

		array2<Complex> f(nx,ny,align);

		fft2d Forward2(-1,f);

		for(unsigned int i=0; i < nx; i++) 
		  for(unsigned int j=0; j < ny; j++) 
			f(i,j)=thedata->data[0][2*(nx/2*(i/2)+j/2)];
		Forward2.Shift(f,nx,ny,1,0);
		Forward2.fft(f);

		double  finish[nx*ny];
		for(unsigned int i=0; i < nx; i++) 
		  for(unsigned int j=0; j < ny; j++) 
			finish[nx*i+j]=sqrt(f(i,j).real()*f(i,j).real()+f(i,j).imag()*f(i,j).imag());

		double  old[nx*ny];
		for(unsigned int i=0; i < nx; i++) 
		  for(unsigned int j=0; j < ny; j++) 
			old[nx*i+j]=thedata->data[0][2*(nx/2*(i/2)+j/2)];

		TString name=Form("%s%d%s%d%s","data/",imap,"/",ni,".bmp");
		c1->Clear();
		draw(finish,old,name,c1,h1,h2,nx,ny);
	}
}

int main(int argc, char **argv) {	
	int start=atoi(argv[1]);
	int theend=atoi(argv[2]);
	gStyle->SetPalette(53,0);
	gStyle->SetNumberContours(256);

	TCanvas *c1=new TCanvas("c1","",800,400);
	TH2D *h1=new TH2D("h1","",40,0,40,40,0,40);
	TH2D *h2=new TH2D("h2","",40,0,40,40,0,40);
	h2->SetMaximum(100);
	h2->SetMinimum(0);
	for(int i=start;i<theend;i++)
	{   
		cout<<i<<endl;
		ftttest(40,i,c1,h1,h2);
	}  
	return 0;
}
