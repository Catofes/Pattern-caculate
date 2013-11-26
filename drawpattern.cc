#include<fstream>
#include"TCanvas.h"
#include"TGraph2D.h"
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
   gStyle->SetPalette(200,colors);
}

void draw(int imap,TCanvas* c1){

	ifstream fin;
	TString name=Form("%s%d%s","../pattern/",imap,".txt");
	cout<<name.Data()<<endl;
	fin.open(name);
	char *header=new char[400];
	for(int ihead=0;ihead<4;ihead++){
		getline(fin,header);
	}
	Double_t *a=new Double_t[400];
	Double_t *b=new Double_t[400];
	Double_t *x=new Double_t[400];
	Double_t *y=new Double_t[400];
	for(int iy=0;iy<20;iy++){
		for(int jx=0;jx<20;jx++){
		double tem=0;
			fin>>tem;
			a[iy*20+jx]=tem;
			x[iy*20+jx]=jx;
			y[iy*20+jx]=iy;
		}
	}
	
	char aa,bb;
	fin>>aa>>bb;
	
	for(int iy=0;iy<20;iy++){
		for(int jx=0;jx<20;jx++){
		double tem=0;
			fin>>tem;
			b[iy*20+jx]=tem;
		}
	}
	cout<<imap<<endl;
	TGraph2D *graphA=new TGraph2D(400,x,y,a);
	graphA->SetName("A");
	graphA->SetTitle(name);
	graphA->SetNpx(20);
	graphA->SetNpy(20);
	graphA->SetMaximum(1.0001);
	graphA->SetMinimum(-0.0001);
	TGraph2D *graphB=new TGraph2D(400,x,y,a);
	graphB->SetName("B");
	graphB->SetTitle(name);
	graphB->SetNpx(20);
	graphB->SetNpy(20);
	graphB->SetMaximum(1.0001);
	graphB->SetMinimum(-0.0001);
	
	c1->cd(1);
	graphA->Draw("COLZ");
	//TExec *ex1 = new TExec("ex1","Pal();");
	//ex1->Draw();
	//graphA->Draw("COLZ SAME");
	c1->Update();
	c1->cd(2);
	graphB->Draw("COLZ");
	//TExec *ex2 = new TExec("ex2","Pal();");
	//ex2->Draw();
	//graphB->Draw("COLZ SAME");
	c1->Update();
	
	TString plotName=Form("%s%d%s","../pattern/",imap,"A.bmp");
	c1->Print(plotName);
	//getchar();
}
void main(){
	TCanvas *c1=new TCanvas("c1","",600,300);
	//Pal();
	gStyle->SetPalette(53,0);
	gStyle->SetNumberContours(256);
	c1->Divide(2,1);
	for(int i=0;i<81;i++){
		draw(i,c1);
	}
}
