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
void draw(Double_t *data,int imap, Int_t x=20, Int_t y=20){
	TCanvas *c1=new TCanvas("c1","",400,400);
	TH2D *h=new TH2D("h","",x,0,x,y,0,y);
	for(Int_t i=0;i<x*y;i++){
		for(Int_t jhight=0;jhight<data[i]*100;jhight++){
			h->Fill(i%x+0.5,i/y+0.5);
		}
	}
	h->Draw("COLZ");
	TString name=Form("%s%d%s","data/",imap,".fft.bmp");
	c1->Print(name);
}

void drawraw(Double_t *data,int imap, Int_t x=20, Int_t y=20){
	TCanvas *c1=new TCanvas("c1","",400,400);
	TH2D *h=new TH2D("h","",x,0,x,y,0,y);
	for(Int_t i=0;i<x*y;i++){
		for(Int_t jhight=0;jhight<data[i]*100;jhight++){
			h->Fill(i%x+0.5,i/y+0.5);
		}
	}   
	h->Draw("COLZ");
	TString name=Form("%s%d%s","data/",imap,".bmp");
	c1->Print(name);
}

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

	for(unsigned int i=0; i < nx; i++) 
	  for(unsigned int j=0; j < ny; j++) 
		f(i,j)=thedata->data[0][2*(20*i+j)];

	Forward2.Shift(f,nx,ny,1);
	Forward2.fft(f);
	cout<<f<<endl;

	double * finish=new double[400];
	for(unsigned int i=0; i < nx; i++) 
	  for(unsigned int j=0; j < ny; j++) 
		finish[20*i+j]=sqrt(f(i,j).real()*f(i,j).real()+f(i,j).imag()*f(i,j).imag());

	double * old=new double[400];
	for(int i=0;i<400;i++){
		//cout<<finish[i]<<"	";
		old[i]=thedata->data[0][2*i];
		cout<<old[i]<<"	";
	}

	draw(finish,imap);
	drawraw(old,imap);

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
