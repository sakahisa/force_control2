#include "Eigen/Dense"

using namespace Eigen;

#define twidth 	0.000001	//sample time
#define TMAX 	0.001		//total time

// length of left leg links
#define L1 		0.1
#define L2 		0.5
#define L3 		0.5

// length of right leg links
#define R1		0.1
#define R2		0.5
#define R3		0.5

// mass of links
#define M1 		0.2
#define M2 		1.0
#define M3 		1.0

/*
#define l1 		0.1
#define l2 		0.5
#define l3 		0.5

#define M 1.0
*/
#define W 		1000			//cut off

#define D 		2.5			//coefficient of viscosity

#define TRANS(x) x.block<3,1>(0,3)
#define ROT(x) x.block<3,3>(0,0) 



/*                                  L1           |       R2
 * 							○__________________Base________________○
 * 							|										 |
 * 							|										 |
 * 							| 										 |
 * 							|  L2									 |  R2
 * 							|										 |
 * 							|										 |
 * 							|           							 |
 * 							○										○
 * 							|										 |
 * 							|										 |
 * 							| 										 |
 * 							|  L3									 |  R3
 * 							|										 |
 * 							|										 |
 * 							|           							 |
 *                        --○--								　--○--
 */
