#include <iostream>
#include <cmath>
#include <fstream>

#define twidth 0.00001
#define TMAX 0.04
#define W 500

using namespace std;

double LPF(double insig, double outsig)
{
	double dsig = (insig - outsig)*W;	
	outsig += dsig*twidth;
	return outsig;
}

int main(){
	ofstream ofs_p("signal.txt");
	double x, y, z, insig, outsig=0.0;
	
	for(double t = 0.0; t < TMAX; t += twidth)
	{
		x = 2.0 * sin(2 * M_PI * 20 * t);
		y = 0.5 * sin(2 * M_PI * 20000 * t);
		z = 2.0 * cos(2 * M_PI * 4000 * t);
		
		insig = x + y + z;
		outsig = LPF(insig, outsig);
		
		ofs_p << t << " " << insig << " " << outsig << " " << x + z << endl;
	}
	
}
