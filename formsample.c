#include<fstream>
void main(){
	ifstream fin;
	fin.open("sample.txt");
	double a;
	ofstream fout;
	fout.open("sample");
	for(int i=0;i<2000;i++){
		for(int j=0;j<20;j++){
			fin>>a;
			fout<<a<<"	";
		}
		fout<<endl;
	}
	fout.close();
	fin.close();
}
