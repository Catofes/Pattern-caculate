#include "TH2D.h"
#include "TCanvas.h"
#include "iostream"
using namespace std;
void draw(Double_t *data, Int_t x=20, Int_t y=20){
	TCanvas *c1=new TCanvas("c1","",400,400);
	TH2D *h=new TH2D("h","",x,0,x,y,0,y);
	for(Int_t i=0;i<x*y;i++){
		for(Int_t jhight=0;jhight<data[i]*100;jhight++){
			h->Fill(i%x+0.5,i/y+0.5);
		}
	}
	h->Draw("COLZ");
	c1->Print("test.bmp");
}
int main(){
	int i=3,j=5;
	cout<<i/j<<endl;
	Double_t *test=new Double_t[400];
	for(Int_t ix=0;ix<20;ix++){
		for(Int_t iy=0;iy<20;iy++){
			test[iy*20+ix]=(iy+ix)%2==1?0:1;
		}
	}
	draw(test);
	return 0;
}
