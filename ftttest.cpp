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
#include"TF1.h"
#include"TMath.h"


using namespace std;
using namespace Array;
using namespace fftwpp;

ofstream fout;
ClassImp(solve_data);
ClassImp(input_sample);
void draw(Double_t *data,Double_t *datain,TString name,TCanvas *c1,TH2D *h,TH2D *hin,Int_t x=20, Int_t y=20){
	c1->Divide(2,1);
	c1->cd(1);
	h->Reset();
	h->SetStats(kFALSE);
	h->SetMinimum(-0.0001);
	for(Int_t i=0;i<x*y;i++){
		for(Int_t jhight=0;jhight<data[i]*100;jhight++){
			h->Fill(i%x+0.5,i/y+0.5);
		}
	}
	h->Draw("COLZ");
	c1->cd(2);
	hin->Reset();
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



int fitfather(double * data,TH1D * h,TF1 * f,int n)
{
	int time=0;
	double mean=-1000;
	while(time<10&&(mean<0||mean>n)){
		time++;
		h->Reset();
		for(int i=0;i<n;i++)
		  for(int k=0;k<data[i];k++)
			h->Fill(i);
		f->SetParameters(h->GetEntries(),n/2+1,1);
		h->Fit(f,"LQ","",0,n);
		mean=f->GetParameter(1);
	}
	fout<<f->GetParameter(0)<<"	"<<f->GetParameter(1)<<"	"<<f->GetParameter(2)<<"	";
}
int datafather(double * old,TH1D * h,TF1 * f,int n,int x,int y,int nx,int ny)
{
	double data[n];
	for(int i=0;i<n;i++)
	  data[i]=0;
	for(int i=0;i<n;i++)
	  for(int k=0;k<n;k++)
		data[i]+=old[nx*(x+n/2-i)+y+n/2-k];
	fitfather(data,h,f,n);
	for(int i=0;i<n;i++)
	  data[i]=0;
	for(int i=0;i<n;i++)
	  for(int k=0;k<n;k++)
		data[i]+=old[nx*(x+n/2-k)+y+n/2-i];
	fitfather(data,h,f,n);
}


int ftttest(int size,int imap,TCanvas *c1,TH2D *h1,TH2D *h2,TH1D *h3,TF1 *fun,int nbin) {
	TString rootname=Form("%s%d%s","data/",imap,".root");
	TFile *f1=new TFile(rootname);
	TTree *tree=(TTree*)f1->Get("tree");
	solve_data *thedata=0;
	tree->SetBranchAddress("thedata",&thedata);

	int nn=tree->GetEntries();
	for ( int ni=0;ni<nn;ni++){
		if(ni%100==1)
		  cout<<imap<<"	"<<(ni+0.01)/nn<<endl;
		tree->GetEntry(ni);
		//	tree->Show();

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
		for(unsigned int i=0; i < nx; i++) {
			for(unsigned int j=0; j < ny; j++){ 
				finish[nx*i+j]=sqrt(f(i,j).real()*f(i,j).real()+f(i,j).imag()*f(i,j).imag());
			}
		}

		double  old[nx*ny];
		for(unsigned int i=0; i < nx; i++) 
		  for(unsigned int j=0; j < ny; j++) 
			old[nx*i+j]=thedata->data[0][2*(nx/2*(i/2)+j/2)];

		//TString name=Form("%s%d%s%d%s","data/",imap,"/",ni,".png");
		//c1->Clear();
		//draw(finish,old,name,c1,h1,h2,nx,ny);
		TString dataname=Form("%s%d%s","data/",imap,".data");
		fout<<imap<<"	"<<ni/3<<"	"<<ni%3<<"	";
		fout<<thedata->params.t_a<<"  "<<thedata->params.t_b<<"  "<<thedata->params.alpha<<"  "<<thedata->params.beta<<"  ";
		for(int mp=0;mp<8;mp++)
		  fout<<thedata->params.k[mp]<<"  ";

		for(int mp=0;mp<8;mp++)
		  fout<<thedata->params.n[mp]<<"  ";

		datafather(finish,h3,fun,nbin,nx/2,ny/2,nx,ny);
		datafather(finish,h3,fun,nbin,nx/4,ny/4,nx,ny);
		fout<<endl;
	}
}


int main(int argc, char **argv) {	
	int start=atoi(argv[1]);
	int theend=atoi(argv[2]);
	gStyle->SetPalette(53,0);
	gStyle->SetNumberContours(256);
	fout.open("data/analsys");
	TCanvas *c1=new TCanvas("c1","",800,400);
	TH2D *h1=new TH2D("h1","",40,0,40,40,0,40);
	TH2D *h2=new TH2D("h2","",40,0,40,40,0,40);
	TH1D *h3=new TH1D("h3","",11,0,11);
	TF1 *f=new TF1("f","abs([0])*TMath::Gaus(x,[1],[2])",0,11);
	h2->SetMaximum(100.01);
	h2->SetMinimum(-0.01);
	fout<<"imap"<<"	"<<"参数编号"<<"	"<<"初始值编号"<<"	";
	fout<<"t_a"<<"	"<<"t_b"<<"	"<<"alpha"<<"	"<<"beta"<<"	";
	for(int mp=0;mp<8;mp++)
	  fout<<"k["<<mp<<"]"<<"	";

	for(int mp=0;mp<8;mp++)
	  fout<<"n["<<mp<<"]"<<"	";
	for(int mp=0;mp<4;mp++)
	  fout<<"a"<<"	"<<"mean"<<"	"<<"theta"<<"	";
	fout<<endl;
	for(int i=start;i<theend;i++)
	{   
		cout<<i<<endl;
		ftttest(40,i,c1,h1,h2,h3,f,11);
	}
	fout.close();
	return 0;
}
