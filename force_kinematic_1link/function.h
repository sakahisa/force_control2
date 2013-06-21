#include<iostream>
#include<cmath>
#include"constant.h"
/*
//---------------------------------------回転移動行列の生成--------------------------------------------------------
void Trans_round(int i,double the,double R[3][3]){
	if(i==1){ //x軸方向
		R[0][0]=1;
		R[0][1]=0;
		R[0][2]=0;
		R[1][0]=0;
		R[1][1]=cos(the);
		R[1][2]=-sin(the);
		R[2][0]=0;
		R[2][1]=sin(the);
		R[2][2]=cos(the);
	}
	else if(i==2){ //y軸方向
		R[0][0]=cos(the);
		R[0][1]=0;
		R[0][2]=sin(the);
		R[1][0]=0;
		R[1][1]=1;
		R[1][2]=0;
		R[2][0]=-sin(the);
		R[2][1]=0;
		R[2][2]=cos(the);
	}
	else if(i==3){ //z軸方向
		R[0][0]=cos(the);
		R[0][1]=-sin(the);
		R[0][2]=0;
		R[1][0]=sin(the);
		R[1][1]=cos(the);
		R[1][2]=0;
		R[2][0]=0;
		R[2][1]=0;
		R[2][2]=1;
	}
}
*/
/*
//---------------------------------------行列の積((3×3)・(3×3))------------------------------------------------------
void cross_33(double A[N][N],double B[N][N], double C[N][N]){
	int i,j,k;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			C[i][j]=0.0;
		}
	}
	
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			for(k=0;k<N;k++){
				C[i][j]+=A[i][k]*B[k][j];
			}
		}
	}
}
*/
//---------------------------------------行列の積((3×3)・(3×1))------------------------------------------------------
void cross_31(double A[N][N],double B[N], double C[N]){
	int i,j;
	
	for(i=0;i<N;i++){
		C[i]=0.0;
	}
	
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			C[i]+=A[i][j]*B[j];
		}
	}
}
/*
//---------------------------------------同次変換行列の生成---------------------------------------------------------
void Trans(double r[N],double R[N][N],double T[N][N]){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			T[i][j]=R[i][j];
		}
	}
	for(i=0;i<3;i++){
		T[i][3]=r[i];
	}
	for(j=0;j<4;j++){
		if(j!=3){
			T[3][j]=0;
		}else{
			T[3][j]=1;
		}
	}
}*/

//---------------------------------ヤコビアンの転置行列の生成-------------------------------------------------------
void Jaco_t(double r,double the,double phi,double J[3][3]){
	J[0][0]=sin(phi)*cos(the);
	J[0][1]=sin(phi)*sin(the);
	J[0][2]=cos(phi);
	J[1][0]=r*cos(phi)*cos(the);
	J[1][1]=r*cos(phi)*sin(the);
	J[1][2]=-r*sin(phi);
	J[2][0]=-r*sin(phi)*sin(the);
	J[2][1]=r*sin(phi)*cos(the);
	J[2][2]=0.0;
}

/*---------------------------------ヤコビアンの転置行列の逆行列の生成-----------------------------------------------
void Jaco_in(double r,double the,double phi,double J[3][3]){
	J[0][0]=cos(the);
	J[0][1]=cos(phi)*cos(the)/(r*sin(phi));
	J[0][2]=-1/r;
	J[1][0]=sin(the);
	J[1][1]=cos(phi)*sin(the)/(r*sin(phi));
	J[1][2]=r*cos(the);
	J[2][0]=-r*sin(phi)*sin(the);
	J[2][1]=r*sin(phi)*cos(the);
	J[2][2]=0;
}*/

void Jaco_in(double r,double the,double phi,double J[3][3]){
	J[0][0]=sin(phi)*cos(the);
	J[1][0]=sin(phi)*sin(the);
	J[2][0]=cos(phi);
	J[0][1]=cos(phi)*cos(the)/r;
	J[1][1]=cos(phi)*sin(the)/r;
	J[2][1]=-sin(phi)/r;
	J[0][2]=-sin(the)/(r*sin(phi));
	J[1][2]=cos(the)/(r*sin(phi));
	J[2][2]=0.0;
}
