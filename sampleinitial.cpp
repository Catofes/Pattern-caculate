#include"fstream"
#include"iostream"
#include"TRandom.h"
using namespace std;
int main(){
	ofstream fout;
	fout.open("initial3");
	Double_t x;
	for(int iset=0;iset<3;iset++){
		for(int i=0;i<800;i++){
			x=gRandom->Rndm();
			fout<<x<<"	";
		}
		fout<<endl;
	}
fout.close();
}
