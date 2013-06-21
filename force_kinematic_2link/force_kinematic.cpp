#include<iostream>
#include<fstream>
#include<cmath>
#include"constant.h"
#include"function.h"
#include"state.h"
#include"force.h"

using namespace std;

int main(){
	std::ofstream ofs("force_res.dat");
	double t;
	state M1,M2;
	force F1,F2;
	
	int i;
	
	for(t=0.0;t<TMAX;t+=twidth){
		Jaco_t(r,the,phi,J);
		if(sin(phi)<0.000001 && sin(phi)>-0.000001){
			std::cout<<"特異点に到達しました。\n";
			break;
		}else{
			Jaco_in(r,the,phi,J_in);
		}
		
		for(i=0;i<N;i++){
			F_dis[i]=v[i]*D;
		}
		
		F_dis[2]+=m*g;					
		
		cross_31(J,F_dis,tau_dis);		//トルク変換(外乱)
		tau_dis[0]=0.0;					//径方向の外乱０
		
		for(i=0;i<N;i++){
			F_ref[i]=F_cmd[i]+F_dist[i];
		}
		cross_31(J,F_ref,tau_ref);		//トルク変換
		tau_ref[0]=0.0;					//径方向の力０
		
		for(i=0;i<N;i++){
			tau_res[i]=tau_ref[i]-tau_dis[i];
		}
		tau_res[0]=0.0;					//径方向の力０
		cross_31(J_in,tau_res,F_res);
		
		for(i=0;i<N;i++){
			a[i]=F_res[i]/m;
		}
		for(i=0;i<N;i++){
			v[i]+=a[i]*twidth;
		}
		for(i=0;i<N;i++){
			x[i]+=v[i]*twidth;
		}
		phi=acos(x[2]/l);
		the=asin(x[1]/(l*sin(phi)));
		p=sqrt(pow(x[0],2.0)+pow(x[1],2.0)+pow(x[2],2.0));
		for(i=0;i<N;i++){
			x[i]/=p;
		}
		
		
		ofs << t << " " << tau_res[0] << " " << tau_res[1] << " " << tau_res[2] <<" "<< x[0] << " " << x[1] << " " << x[2] << " " << p << " " << phi << " " << the <<std::endl;
		
		for(i=0;i<N;i++){
			F_dif[i]=F_ref[i]-F_res[i];
		}
		for(i=0;i<N;i++){
			dF_dist[i]+=F_dist[i]*twidth*W;
		}
		for(i=0;i<N;i++){
			F_dist[i]=F_dif[i]-dF_dist[i];
		}
		
	}
	return 0;
}
