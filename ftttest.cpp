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
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cstddef>

using namespace std;
using namespace Array;
using namespace fftwpp;

ClassImp(solve_data);
ClassImp(input_sample);
void draw(Double_t *data,TString name, Int_t x=20, Int_t y=20){
	static Int_t colors[256];
	static Bool_t initialized = kFALSE;
	Double_t Red[3]  ={1.00,0.50,0.00};
	Double_t Green[3]={1.00,0.50,0.00};
	Double_t Blue[3] ={1.00,0.50,0.00};
	Double_t Length[3]={0.00,0.50,1.00};
	if(!initialized){
		Int_t FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 256);
		for(int i=0;i<256;i++)colors[i]=FI+i;
		initialized = kTRUE;
	}
	gStyle->SetNumberContours(256);
	gStyle->SetPalette(256,colors);
	TCanvas *c1=new TCanvas("c1","",400,400);
	//	gStyle->SetC
	TH2D *h=new TH2D("h","",x,0,x,y,0,y);
	h->SetStats(kFALSE);
	h->SetMaximum(101);
	h->SetMinimum(0);
	for(Int_t i=1;i<x*y;i++){
		for(Int_t jhight=0;jhight<data[i]*100;jhight++){
			h->Fill(i%x+0.5,i/y+0.5);
		}
	}
	h->Draw("COLZ");
	c1->Print(name);
}

void drawraw(Double_t *data,TString name, Int_t x=20, Int_t y=20){
	TCanvas *c1=new TCanvas("c1","",400,400);
	gStyle->SetPalette(53,0);
	gStyle->SetNumberContours(256);
	TH2D *h=new TH2D("h","",x,0,x,y,0,y);
	h->SetStats(kFALSE);
	h->SetMaximum(101);
	h->SetMinimum(0);
	for(Int_t i=0;i<x*y;i++){
		for(Int_t jhight=0;jhight<data[i]*100;jhight++){
			h->Fill(i%x+0.5,i/y+0.5);
		}
	}   
	h->Draw("COLZ");
	c1->Print(name);
}

int ftttest(int imap) {
	cout<<"0.1"<<endl;
	TString rootname=Form("%s%d%s","data/",imap,".root");
	TFile *f1=new TFile(rootname);
	TTree *tree=(TTree*)f1->Get("tree");
	solve_data *thedata=0;
	tree->SetBranchAddress("thedata",&thedata);
	int nen=tree->GetEntries();
	for(int count=0;count<nen;count++){
		cout<<(count-0.1)/nen<<endl;
		tree->GetEntry(count);
		fftw::maxthreads=get_max_threads();
		unsigned int nx=20, ny=20;
		size_t align=sizeof(Complex);
		array2<Complex> f(nx,ny,align);
		fft2d Forward2(-1,f);
		for(unsigned int i=0; i < nx; i++) 
		  for(unsigned int j=0; j < ny; j++) 
			f(i,j)=thedata->data[0][2*(20*i+j)];
		Forward2.fft(f);
		double * finish=new double[400];
		for(unsigned int i=0; i < nx; i++) 
		  for(unsigned int j=0; j < ny; j++) 
			finish[20*i+j]=sqrt(f(i,j).real()*f(i,j).real()+f(i,j).imag()*f(i,j).imag());
		double * old=new double[400];
		for(int i=0;i<400;i++){
			//cout<<finish[i]<<"	";
			old[i]=thedata->data[0][2*i];
		}
		TString name0=Form("%s%d%s%d%s","data/",imap,"/",count,".bmp");
		TString name1=Form("%s%d%s%d%s","data/",imap,"/",count,".fft.bmp");
		draw(finish,name0);
		drawraw(old,name1);
	}
}
int main(int argc, char **argv) {
	int start=atoi(argv[1]);
	int theend=atoi(argv[2]);
	gStyle->SetPalette(53,0);
	gStyle->SetNumberContours(256);

	for(int i=start;i<theend;i++)
	{   
		cout<<i<<endl;
		ftttest(i);
	}  
	return 0;
}
