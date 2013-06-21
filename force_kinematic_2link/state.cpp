#include"state.h"

state::state(){
	int i;
	for(i=0;i<3;i++){
		dx[i]=0.0;
		ddx[i]=0.0;
		dr[i]=0.0;
	}
}
