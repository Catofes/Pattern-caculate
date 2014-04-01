#include<iostream>
#include<stdlib.h>
#include<fstream>
#include"TCanvas.h"
#include"TGraph2D.h"
#include"TColor.h"
#include"TStyle.h"
#include"TH2.h"
using namespace std;
void Pal()//set the drawing color
{
   static Int_t  colors[256];
   static Bool_t initialized = kFALSE;

   Double_t Red[3]    = { 1.00, 0.50, 0.00};
   Double_t Green[3]  = { 1.00, 0.50, 0.00};
   Double_t Blue[3]   = { 1.00, 0.50, 0.00};
   Double_t Length[3] = { 0.00, 0.50, 1.00};

   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,256);
      for (int i=0; i<256; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetNumberContours(256);
   gStyle->SetPalette(256,colors);
}

void draw(int imap,TCanvas* c1){

	ifstream fin;
	TString name=Form("%s%d%s","data/",imap,".txt");
	cout<<name.Data()<<endl;
	fin.open(name);
	char *header=new char[400];
	for(int ihead=0;ihead<4;ihead++){
		fin.getline(header,400);
	}
	Double_t *a=new Double_t[402];
	Double_t *b=new Double_t[402];
	Double_t *fft=new Double_t[402];
	Double_t *x=new Double_t[402];
	Double_t *y=new Double_t[402];
	for(int iy=0;iy<20;iy++){
		for(int jx=0;jx<20;jx++){
		double tem=0;
			fin>>tem;
			a[iy*20+jx]=tem;
			//cout<<tem<<" ";
			x[iy*20+jx]=jx+0.5;
			y[iy*20+jx]=iy+0.5;
		}
		//cout<<endl;
	}
	//fin.getline(header,1);
	char bb[2];
	fin>>bb[0]>>bb[1];	
	//cout<<"***************"<<endl;
	for(int iy=0;iy<20;iy++){
		for(int jx=0;jx<20;jx++){
		double tem=0;
			fin>>tem;//cout<<tem<<" ";
			b[iy*20+jx]=tem;
		}
		//cout<<endl;
	}
	//fin.getline(header,1);
	char cc[8];
	fin>>cc[0]>>cc[1]>>cc[2]>>cc[3]>>cc[4]>>cc[5]>>cc[6]>>cc[7];
	//cout<<"**************"<<endl;
	for(int iy=0;iy<20;iy++){
		for(int jx=0;jx<20;jx++){
		double tem=0;
			fin>>tem;
			//cout<<tem<<"  ";
			fft[iy*20+jx]=tem;
		}
		//cout<<endl;
		//getchar();
	}
	a[400]=a[0];a[401]=a[399];b[400]=b[0];b[401]=b[399];fft[400]=fft[0];fft[401]=fft[399];x[400]=y[400]=0;x[401]=y[401]=20;
	
	cout<<imap<<endl;
	TGraph2D *graphA=new TGraph2D(402,x,y,a);
	graphA->SetName("A");
	graphA->SetTitle(name);
	//graphA->Set
	graphA->SetNpx(20);
	graphA->SetNpy(20);
	//graphA->GetXaxis()->Set(20,0,20);
	//graphA->GetYaxis()->Set(20,0,20);
	graphA->SetMaximum(1.0001);
	graphA->SetMinimum(-0.0001);
	TGraph2D *graphB=new TGraph2D(402,x,y,b);
	graphB->SetName("B");
	graphB->SetTitle(name);
	graphB->SetNpx(20);
	graphB->SetNpy(20);
	//graphB->GetXaxis()->SetLimits(0,20);
	//graphB->GetYaxis()->SetLimits(0,20);
	graphB->SetMaximum(1.0001);
	graphB->SetMinimum(-0.0001);
	TGraph2D *graphFFT=new TGraph2D(402,x,y,fft);
	graphFFT->SetName("ArowFFT");
	graphFFT->SetTitle(name);
	graphFFT->SetNpx(20);
	graphFFT->SetNpy(20);
	//graphFFT->GetXaxis()->SetLimits(0,20);
	//graphFFT->GetYaxis()->SetLimits(0,20);
	graphFFT->SetMaximum(10.);
	graphFFT->SetMinimum(-10.);
	
	c1->cd(1);
	graphA->Draw("COLZ");
	c1->Update();
	c1->cd(2);
	graphB->Draw("COLZ");
	c1->Update();
	c1->cd(3);
	graphFFT->Draw("COLZ");
	c1->Update();
	
	TString plotName=Form("%s%d%s","data/",imap,".bmp");
	c1->Print(plotName);
}
int main(int argc, char **argv)
{
	int start=atoi(argv[1]);
	int theend=atoi(argv[2]);
	TCanvas *c1=new TCanvas("c1","",600,600);
	gStyle->SetPalette(53,0);
	gStyle->SetNumberContours(256);
	c1->Divide(2,2);
	for(int i=start;i<theend;i++){
		draw(i,c1);
	}
}
