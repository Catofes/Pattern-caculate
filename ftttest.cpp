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
void draw(Double_t *data,Double_t *datain,int imap, Int_t x=20, Int_t y=20){
	TCanvas *c1=new TCanvas("c1","",800,400);
	c1->Divide(2,1);
	c1->cd(1);
	TH2D *h=new TH2D("h","",x,0,x,y,0,y);
	h->SetStats(kFALSE);
	h->SetMinimum(-0.0001);
	for(Int_t i=0;i<x*y;i++){
		for(Int_t jhight=0;jhight<data[i]*100;jhight++){
			h->Fill(i%x+0.5,i/y+0.5);
		}
	}
	h->Draw("COLZ");
	c1->cd(2);
	TH2D *hin=new TH2D("hin","",x,0,x,y,0,y);
	hin->SetStats(kFALSE);
	hin->SetMinimum(-0.0001);
	for(Int_t i=0;i<x*y;i++){
		for(Int_t jhight=0;jhight<datain[i]*100;jhight++){
			hin->Fill(i%x+0.5,i/y+0.5);
		}
	}
	hin->Draw("COLZ");
	
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

int ftttest(int size=20,int imap=0) {
	TString rootname=Form("%s%d%s","data/",imap,".root");
	TFile *f1=new TFile(rootname);
	TTree *tree=(TTree*)f1->Get("tree");
	solve_data *thedata=0;
	tree->SetBranchAddress("thedata",&thedata);
	tree->GetEntry(0);
	tree->Show();

//	fftw::maxthreads=get_max_threads();
	unsigned int nx=size,ny=size;//nx=20, ny=20;
	size_t align=sizeof(Complex);

	Double_t *datain=new Double_t[nx*ny];
	array2<Complex> f(nx,ny,align);

	fft2d Forward2(-1,f);

	for(unsigned int i=0; i < nx; i++) 
	  for(unsigned int j=0; j < ny; j++) 
		//f(i,j)=thedata->data[0][2*(20*i+j)];
		{
			if(sqrt(pow((i-nx/2.),2)+pow((j-ny/2.),2))>nx/10){
				f(i,j)=(i/2+j/2)%2==0?0:1;
				datain[i*nx+j]=f(i,j).real();
			}
			else {
				f(i,j)=(i/2+j/2)%2==0?1:0;
				datain[i*nx+j]=f(i,j).real();
			}
	  }
	Forward2.Shift(f,nx,ny,1,0);
	Forward2.fft(f);
	//cout<<f<<endl;

	double * finish=new double[nx*ny];
	for(unsigned int i=0; i < nx; i++) 
	  for(unsigned int j=0; j < ny; j++) 
		finish[nx*i+j]=sqrt(f(i,j).real()*f(i,j).real()+f(i,j).imag()*f(i,j).imag());

	double * old=new double[nx*ny];
	for(int i=0;i<nx*ny;i++){
		//cout<<finish[i]<<"	";
		old[i]=thedata->data[0][2*i];
	}

	draw(finish,datain,imap,nx,ny);

}

int main(int argc, char **argv) {	
	//int start=atoi(argv[1]);
	//int theend=atoi(argv[2]);
	int size=atoi(argv[1]);
	gStyle->SetPalette(53,0);
	gStyle->SetNumberContours(256);

	//for(int i=start;i<theend;i++)
	{   
//		cout<<i<<endl;
		ftttest(size,0);
	}  
	return 0;
}
