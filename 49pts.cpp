/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 */

#include "minigl.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>
#include <stack>
#include <map>
#include <iostream>

using namespace std;


/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//VECTOR of Vector of Verticies
vector< vector < vector<MGLfloat> >> GLOBAL_VERTICIES;	

//Global temporary "list" of verticies to be put into the vector above this one
vector<vector<MGLfloat>> TEMP_LIST_VERTICIES;

//GLOBAL MATRIX STACK
stack<vector<MGLfloat>> VMStack; 
stack<vector<MGLfloat>> ProjStack;

vector<MGLfloat> VMCurrentMatrix;
vector<MGLfloat> ProjCurrentMatrix;

MGLint MatrixMode;	//if 1 --> proj if 2 -->VM

//Varible for TYPE OF vertex
MGLfloat GLOBAL_VTYPE;

// GLOBAL COLOR for verticies
MGLfloat GLOBAL_COLOR[3] = {1,1,1};

//MAX FAR AND NEAR VALUES OF Z
MGLfloat Zfar = 10000;
MGLfloat Znear = -10000;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Current matrix * VERTEX
vector<MGLfloat> MultMatrixCUR(vector<MGLfloat> i, MGLint MODE );
//dot product of 2 matricies
MGLfloat dot(vector<MGLfloat> X, vector<MGLfloat> Y);
//4x4 matrix * 4x4 matrix saved in current
vector<MGLfloat> MatrixMult4x4(vector<MGLfloat> i, vector<MGLfloat> j);


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
	//STOREING Z VALUES 
	map<pair<int, int>, float> Z_Map;
	
	
	//STEP 1	Find The Origin in Data Value
	MGLpixel NewOrigin = width *( height +1);
	//because if you divide an odd number it will round down so we added one for height
	//dividing by 2
	NewOrigin >>= 1;
	if(NewOrigin % width != 0)
	{
		NewOrigin = NewOrigin + (NewOrigin%width);
	}
	//sutracting we want to round down so its good 
	NewOrigin = NewOrigin - (width >> 1); 
	
	
	//Do everything in this for loop b/c one shape at a time
	for(unsigned H = 0; H < GLOBAL_VERTICIES.size(); H++)
	{
		int halfwidth = width >> 1;
		int halfheight = height >> 1;
	
		int Xmax = -(halfwidth);
		int Xmin = halfwidth;
		int Ymin = halfheight;
		int Ymax =  -(halfheight);
		
		/*for(unsigned l = 0; l < 3; l++)
		{
			cout << "x: " << GLOBAL_VERTICIES[H][l][0] << " y: " << GLOBAL_VERTICIES[H][l][1] <<  endl;
		}*/
	//STEP 2 	Convert XY percent to X,Y on the Origin store it back into global
	// 			convertion ... X% * .5 width = X ... Y% * .5(height) = Y	
		for(unsigned j = 0; j < GLOBAL_VERTICIES[H].size(); j++)
		{
	
			GLOBAL_VERTICIES[H][j][0] = GLOBAL_VERTICIES[H][j][0]*(width >> 1);
			GLOBAL_VERTICIES[H][j][1] = GLOBAL_VERTICIES[H][j][1]*(height >> 1);

	
		//GET MIN AND MAX POINTS HERE TOO
			if(GLOBAL_VERTICIES[H][j][0] > Xmax && (GLOBAL_VERTICIES[H][j][0] <= halfwidth) )
			{
				Xmax = GLOBAL_VERTICIES[H][j][0];
			}
			if(GLOBAL_VERTICIES[H][j][0] < Xmin && (GLOBAL_VERTICIES[H][j][0] >= -halfwidth) )
			{
				Xmin = GLOBAL_VERTICIES[H][j][0];
			}
			if(GLOBAL_VERTICIES[H][j][1] > Ymax && (GLOBAL_VERTICIES[H][j][1] <= halfheight) )
			{
				Ymax = GLOBAL_VERTICIES[H][j][1];
			}
			if(GLOBAL_VERTICIES[H][j][1] < Ymin && (GLOBAL_VERTICIES[H][j][1] >= -halfheight))
			{
				Ymin = GLOBAL_VERTICIES[H][j][1];
			}			
		}
	
	//STEP 3	Do Barametric Cordinates on these values ... and change values into data and store
		//1st get the vectors V0 V1  
		vector<MGLfloat> V0, V1,V2;
		V0.resize(3);
		V1.resize(3);
		V2.resize(3);
		
		//dot product results
		MGLfloat dot00, dot01, dot02, dot11, dot12;  
		
		if(GLOBAL_VERTICIES[H].size() >= 3)
		{	
			V0[0] = ( GLOBAL_VERTICIES[H][2][0] - GLOBAL_VERTICIES[H][0][0] ); 
			V0[1] = ( GLOBAL_VERTICIES[H][2][1] - GLOBAL_VERTICIES[H][0][1] );
			//V0[2] = ( GLOBAL_VERTICIES[H][2][2] - GLOBAL_VERTICIES[H][0][2] );
			
			V1[0] = ( GLOBAL_VERTICIES[H][1][0] - GLOBAL_VERTICIES[H][0][0] ); 
			V1[1] = ( GLOBAL_VERTICIES[H][1][1] - GLOBAL_VERTICIES[H][0][1] );
			//V1[2] = ( GLOBAL_VERTICIES[H][1][2] - GLOBAL_VERTICIES[H][0][2] );
			
		}
		
		dot00 = dot(V0, V0);
		dot01 = dot(V0, V1);
		
		dot11 = dot(V1,V1);

		/*
		//CODE USED FOR CHECKING
		for(unsigned l = 0; l < 3; l++)
		{
			cout << "x: " << GLOBAL_VERTICIES[H][l][0] << " y: " << GLOBAL_VERTICIES[H][l][1] << " z: " << GLOBAL_VERTICIES[H][l][2] << " w: " << GLOBAL_VERTICIES[0][l][3] <<  endl;
		}
		cerr << endl;*/
		//cerr << Xmin << ' ' << Xmax << ' ' << Ymin << ' ' << Ymax << endl;
		/*
		MGLpixel a; 
		a =  GLOBAL_VERTICIES[0][2][0] + NewOrigin;
		a = a + (GLOBAL_VERTICIES[0][2][1] * width);
		//make it white not teal
		data[a] = 0xFFFFFFFF;
		
		a =  GLOBAL_VERTICIES[0][1][0] + NewOrigin;
		a = a + (GLOBAL_VERTICIES[0][1][1] * width);
		//make it white not teal
		data[a] = 0xFFFFFFFF;
		
		a =  GLOBAL_VERTICIES[0][0][0] + NewOrigin;
		a = a + (GLOBAL_VERTICIES[0][0][1] * width);
		//make it white not teal
		data[a] = 0xFFFFFFFF;
		*/
		if(Xmin == halfwidth)
		{
			Xmin = -halfwidth;
		}
		if(Xmax == -halfwidth)
		{
			Xmax = halfwidth;
		}
		
		if(Ymin == halfheight)
		{
			Ymin = -halfheight;
		}
		if(Ymax == -halfheight)
		{
			Ymax = halfheight;
		}
		
		
		
		cerr << "Ymin " << Ymin << " Ymax " << Ymax << endl;
		cerr << "Xmin " << Xmin << " Xmax " << Xmax << endl;
		cerr << "NewOrigin " << NewOrigin << " width " << width << endl;
		
		int temp3 = Xmin + NewOrigin;
		cerr << "P1 = "<< temp3 << endl;
		//temp3 = Ymin * width;
		cerr << "P2 = " << temp3 << endl;
		temp3 = temp3 + (Ymin * width);

		cerr << "SHOULD BE 0 = " << temp3 << endl;
		
		//call function for dot product
		for(int j = Ymin; j < Ymax; j++)
		{
			
			for(int i = Xmin; i < Xmax; i++)
			{
				
				V2[0] = i - GLOBAL_VERTICIES[H][0][0];
				V2[1] = j - GLOBAL_VERTICIES[H][0][1];
				
				
				dot02 = dot(V0,V2);
				dot12 = dot(V1,V2);

				MGLfloat divisor = 1/(dot00 * dot11 - dot01 * dot01);
				MGLfloat alpha = (dot11 * dot02 - dot01 * dot12) * divisor;
				MGLfloat beta = (dot00 * dot12 - dot01 * dot02) * divisor;
				
				if((alpha >= 0) && (beta >= 0) && (alpha+beta <= 1))
				{			
					pair<int,int> CHECKER = {i,j};	
					if(Z_Map.find(CHECKER) == Z_Map.end())
					{
						MGLfloat tempZ = (1-alpha-beta) * GLOBAL_VERTICIES[H][0][2] + beta * GLOBAL_VERTICIES[H][1][2] + alpha * GLOBAL_VERTICIES[H][2][2];
						/*
						cerr << "Znear: " << Znear;
						cerr << " Zfar: " << Zfar;
						cerr << " Zval at pt " << tempZ;
						cerr << "vertex A B C ... z val " << GLOBAL_VERTICIES[H][0][2] << ' ' << GLOBAL_VERTICIES[H][1][2] << ' ' << GLOBAL_VERTICIES[H][2][2] << endl; 
						*/
						if(tempZ >= Znear && tempZ <= Zfar)
						{
							Z_Map[CHECKER] = tempZ;
							int temp2 = (j * width);
							MGLpixel temp; 
							temp = i + NewOrigin;
							temp2 = temp + temp2;
							if(temp2 >= 0)
							{
								temp = temp2;
							}

							//Need to do something with the COLOR !!!
							//we have alpha and beta need gama so if we do 1-a-b = g
							MGLfloat R1,G1,B1,R2,G2,B2,R3,G3,B3;						
							
							R1 = (1-alpha-beta) * GLOBAL_VERTICIES[H][0][4];
							G1 = (1-alpha-beta) * GLOBAL_VERTICIES[H][0][5];
							B1 = (1-alpha-beta) * GLOBAL_VERTICIES[H][0][6];
							
							R2 = beta * GLOBAL_VERTICIES[H][1][4];
							G2 = beta * GLOBAL_VERTICIES[H][1][5];
							B2 = beta * GLOBAL_VERTICIES[H][1][6];
							
							R3 = alpha * GLOBAL_VERTICIES[H][2][4];
							G3 = alpha * GLOBAL_VERTICIES[H][2][5];
							B3 = alpha * GLOBAL_VERTICIES[H][2][6];

							//cerr << G1 << ' ' << G2 << ' ' << B3 << endl;
							
							R1 = (R1+R2+R3) ;
							G1 = (G1+G2+G3) ;
							B1 = (B1+B2+B3) ;
							
							//cerr << R1 << ' ' << G1 << ' ' << B1 << endl;
							
							MGLpixel R = 0xFF * R1;
							MGLpixel G = 0xFF * G1;
							MGLpixel B = 0xFF * B1;
							
							
							R = (R << 24) + ( G << 16 )+ (B << 8) + 0xFF; 
							
							data[temp] = R;
						}
					}
					else
					{
						
						float CurZ;
						CurZ = (1-alpha-beta) * GLOBAL_VERTICIES[H][0][2] + beta * GLOBAL_VERTICIES[H][1][2] + alpha * GLOBAL_VERTICIES[H][2][2];
						
						//cerr <<  "CurZ " << CurZ << " Z_map point "<< CHECKER.first << ',' << CHECKER.second << ' ' << Z_Map[CHECKER] << endl;
						if(CurZ <= Zfar && CurZ >= Znear)
						{	
							if(CurZ < Z_Map[CHECKER])
							{
								Z_Map[CHECKER] = CurZ;
								MGLpixel temp; 
								temp = i + NewOrigin;
								temp = temp + (j * width);
								
								//Need to do something with the COLOR !!!
								//we have alpha and beta need gama so if we do 1-a-b = g
								MGLfloat R1,G1,B1,R2,G2,B2,R3,G3,B3;						
									
								R1 = (1-alpha-beta) * GLOBAL_VERTICIES[H][0][4];
								G1 = (1-alpha-beta) * GLOBAL_VERTICIES[H][0][5];
								B1 = (1-alpha-beta) * GLOBAL_VERTICIES[H][0][6];
								
								R2 = beta * GLOBAL_VERTICIES[H][1][4];
								G2 = beta * GLOBAL_VERTICIES[H][1][5];
								B2 = beta * GLOBAL_VERTICIES[H][1][6];
								
								R3 = alpha * GLOBAL_VERTICIES[H][2][4];
								G3 = alpha * GLOBAL_VERTICIES[H][2][5];
								B3 = alpha * GLOBAL_VERTICIES[H][2][6];

								//cerr << G1 << ' ' << G2 << ' ' << B3 << endl;
								
								R1 = (R1+R2+R3) ;
								G1 = (G1+G2+G3) ;
								B1 = (B1+B2+B3) ;
								
								//cerr << R1 << ' ' << G1 << ' ' << B1 << endl;
								
								MGLpixel R = 0xFF * R1;
								MGLpixel G = 0xFF * G1;
								MGLpixel B = 0xFF * B1;
								
								
								R = (R << 24) + ( G << 16 )+ (B << 8) + 0xFF; 
								
								data[temp] = R;
								//data[temp] = 0xFFFFFFFF;
							}
						}
					}
					
					
				
				}
					

				}
			}
		}
		
		
} 


/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	switch(mode)
	{
		case MGL_TRIANGLES:
			GLOBAL_VTYPE = 1;
			break;
				
		case MGL_QUADS:
			GLOBAL_VTYPE = 2;
			break;
				
		default:
			GLOBAL_VTYPE = 0;
			break;
	}
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	
	/*for(unsigned l = 0; l < TEMP_LIST_VERTICIES.size(); l++)
	{
		cout << "x: " << TEMP_LIST_VERTICIES[l][0] << " y: " << TEMP_LIST_VERTICIES[l][1] << " z: " << TEMP_LIST_VERTICIES[l][2] << " w: " << TEMP_LIST_VERTICIES[l][3] <<  endl;
	}
	cout << endl;*/
	

	if(GLOBAL_VTYPE == 2)
	{
		//cerr << "IN SQUARE \n";
		for(unsigned i = 0; i < TEMP_LIST_VERTICIES.size(); i = i + 4 )
		{
			vector <vector<MGLfloat>> T;
			vector<MGLfloat> TEMP_Vertex;
			TEMP_Vertex = TEMP_LIST_VERTICIES[i];
			T.push_back(TEMP_Vertex);
			TEMP_Vertex = TEMP_LIST_VERTICIES[i+1];
			T.push_back(TEMP_Vertex);
			TEMP_Vertex = TEMP_LIST_VERTICIES[i+2];
			T.push_back(TEMP_Vertex);
			GLOBAL_VERTICIES.push_back(T);
			T.clear();
			//cerr << "in TEMP_LIST_VERTICIES Z " << TEMP_LIST_VERTICIES[i][2];
			//cerr << " in END z:"<< TEMP_Vertex[2] << endl; 

			TEMP_Vertex = TEMP_LIST_VERTICIES[i];
			T.push_back(TEMP_Vertex);
			TEMP_Vertex = TEMP_LIST_VERTICIES[i+2];
			T.push_back(TEMP_Vertex);
			TEMP_Vertex = TEMP_LIST_VERTICIES[i+3];
			T.push_back(TEMP_Vertex);
			GLOBAL_VERTICIES.push_back(T);
		}
		TEMP_LIST_VERTICIES.clear();
		GLOBAL_VTYPE = 0;
	}
	else
	{
		for(unsigned i = 0; i < TEMP_LIST_VERTICIES.size(); i = i + 3 )
		{
			vector <vector<MGLfloat>> T;
			vector<MGLfloat> TEMP_Vertex;
			TEMP_Vertex = TEMP_LIST_VERTICIES[i];
			//cerr << "in TEMP_LIST_VERTICIES Z " << TEMP_LIST_VERTICIES[i][2];
			//cerr << " in END z:"<< TEMP_Vertex[2] << endl; 
			T.push_back(TEMP_Vertex);
			TEMP_Vertex = TEMP_LIST_VERTICIES[i+1];
			T.push_back(TEMP_Vertex);
			TEMP_Vertex = TEMP_LIST_VERTICIES[i+2];
			T.push_back(TEMP_Vertex);
			GLOBAL_VERTICIES.push_back(T);
		}
		TEMP_LIST_VERTICIES.clear();
		GLOBAL_VTYPE = 0;
	}
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	vector<MGLfloat> TEMP_Vertex;
	TEMP_Vertex.resize(8);
	TEMP_Vertex[0] = x;
	TEMP_Vertex[1] = y;
	TEMP_Vertex[2] = 0;
	TEMP_Vertex[3] = 1;
	TEMP_Vertex[4] = GLOBAL_COLOR[0];
	TEMP_Vertex[5] = GLOBAL_COLOR[1];
	TEMP_Vertex[6] = GLOBAL_COLOR[2];
	TEMP_Vertex[7] = GLOBAL_VTYPE;
	
	//cerr << "0riginal: " << TEMP_Vertex[0] << ' ' << TEMP_Vertex[1] << ' ' << TEMP_Vertex[2] << ' ' << TEMP_Vertex[3] << endl;
		
	//MV mult 
	MGLint i = 2;
	TEMP_Vertex = MultMatrixCUR(TEMP_Vertex, i);
	
	//cerr << "0: " << TEMP_Vertex[0] << ' ' << TEMP_Vertex[1] << ' ' << TEMP_Vertex[2] << ' ' << TEMP_Vertex[3] << endl;
	
	//proj mult
	i = 1;
	TEMP_Vertex = MultMatrixCUR(TEMP_Vertex, i);
	
	//cerr << "1: "<< TEMP_Vertex[0] << ' ' << TEMP_Vertex[1] << ' ' << TEMP_Vertex[2] << ' ' << TEMP_Vertex[3] << endl;
	
	//Divide points by the weight
	TEMP_Vertex[0] = TEMP_Vertex[0] / TEMP_Vertex[3];
	TEMP_Vertex[1] = TEMP_Vertex[1] / TEMP_Vertex[3];
	
	TEMP_Vertex[3] = 1;
	
		
	//cerr << "2: "<< TEMP_Vertex[0] << ' ' << TEMP_Vertex[1] << " z:" << TEMP_Vertex[2] << ' ' << TEMP_Vertex[3] << endl;	
	
	if(TEMP_Vertex[0] > 1)
	{
		TEMP_Vertex[0] = 1;
	}
	if(TEMP_Vertex[0] < -1)
	{
		TEMP_Vertex[0] = -1;
	}
	if(TEMP_Vertex[1] > 1)
	{
		TEMP_Vertex[1] = 1;
	}
	if(TEMP_Vertex[1] < -1)
	{
		TEMP_Vertex[1] = -1;
	}
	
	//cerr << "inputed z:" << TEMP_Vertex[2] << ' ';
	
	
	//push into temporary list of verticies
	TEMP_LIST_VERTICIES.push_back(TEMP_Vertex);
	//cerr << "after Z stored in LIST "<< TEMP_LIST_VERTICIES.back()[2] << endl;
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	vector<MGLfloat> TEMP_Vertex;
	TEMP_Vertex.resize(8);
	TEMP_Vertex[0] = x;
	TEMP_Vertex[1] = y;
	TEMP_Vertex[2] = z;
	TEMP_Vertex[3] = 1;
	TEMP_Vertex[4] = GLOBAL_COLOR[0];
	TEMP_Vertex[5] = GLOBAL_COLOR[1];
	TEMP_Vertex[6] = GLOBAL_COLOR[2];
	TEMP_Vertex[7] = GLOBAL_VTYPE;
	
	//MV mult 
	MGLint i = 2;
	TEMP_Vertex = MultMatrixCUR(TEMP_Vertex, i);
	
	//proj mult	
	i = 1;
	TEMP_Vertex = MultMatrixCUR(TEMP_Vertex, i);
	
	//Divide points by the weight
	TEMP_Vertex[0] = TEMP_Vertex[0] / TEMP_Vertex[3];
	TEMP_Vertex[1] = TEMP_Vertex[1] / TEMP_Vertex[3];
	TEMP_Vertex[2] = TEMP_Vertex[2] / TEMP_Vertex[3];
	TEMP_Vertex[3] = 1;
	
	//cerr << "vertex 3 after everything "<< TEMP_Vertex[0] << ' ' << TEMP_Vertex[1] << " " << TEMP_Vertex[2] << ' ' << TEMP_Vertex[3] << endl;	

	/*
	if(TEMP_Vertex[0] > 1)
	{
		TEMP_Vertex[0] = 1;
	}
	if(TEMP_Vertex[0] < -1)
	{
		TEMP_Vertex[0] = -1;
	}
	if(TEMP_Vertex[1] > 1)
	{
		TEMP_Vertex[1] = 1;
	}
	if(TEMP_Vertex[1] < -1)
	{
		TEMP_Vertex[1] = -1;
	}*/

	
	TEMP_LIST_VERTICIES.push_back(TEMP_Vertex);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	if(mode == MGL_MODELVIEW)
	{
		MatrixMode = 2;
	}
	else if(mode == MGL_PROJECTION)
	{
		MatrixMode = 1;
	}
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	if(MatrixMode == 1)
	{
		ProjStack.push(ProjCurrentMatrix);
	}
	else if(MatrixMode == 2)
	{
		VMStack.push(VMCurrentMatrix);
	}
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	if(MatrixMode == 1)
	{
		if(!ProjStack.empty())
		{			
			ProjCurrentMatrix = ProjStack.top();
			ProjStack.pop();
		}

	}
	else if(MatrixMode == 2)
	{
		if(!VMStack.empty())
		{
			VMCurrentMatrix = VMStack.top();
			VMStack.pop();
			
		}
		//VMCurrentMatrix = VMStack.top();
	}
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	
	vector<MGLfloat> CurrentMatrix;
	CurrentMatrix.resize(16);
	CurrentMatrix.at(0) = 1;
	CurrentMatrix.at(5) = 1;
	CurrentMatrix.at(10) = 1;
	CurrentMatrix.at(15) = 1;

	if(MatrixMode == 1)
	{
		ProjCurrentMatrix = CurrentMatrix;
	}
	else if(MatrixMode == 2)
	{
		VMCurrentMatrix = CurrentMatrix;
	}
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
	vector<MGLfloat> CurrentMatrix;
	for(unsigned i = 0; i < 16; i++)
	{
		CurrentMatrix.push_back(matrix[i]);
	}

	if(MatrixMode == 1)
	{
		ProjCurrentMatrix = CurrentMatrix;
	}
	else if(MatrixMode == 2)
	{
		VMCurrentMatrix = CurrentMatrix;
	}
	
/*	
cerr << "mglLoadMatrix \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << CurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
	
cerr << "projMatrix \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << ProjCurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
	*/
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
	vector<MGLfloat> CurrentMatrix;
	for(unsigned i = 0; i < 16; i++)
	{
		CurrentMatrix.push_back(matrix[i]);
	}
	
	/*
cerr << "mglMultMatrix \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << CurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
	
cerr << "projMatrix \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << ProjCurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}*/
	
	if(MatrixMode == 1)
	{
		ProjCurrentMatrix = MatrixMult4x4(ProjCurrentMatrix,CurrentMatrix);
	}
	else if(MatrixMode == 2)
	{
		VMCurrentMatrix = MatrixMult4x4(VMCurrentMatrix,CurrentMatrix);
	}
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	vector<MGLfloat> CurrentMatrix;
	CurrentMatrix.resize(16);
	CurrentMatrix.at(0) = 1;
	CurrentMatrix.at(5) = 1;
	CurrentMatrix.at(10) = 1;
	CurrentMatrix.at(15) = 1;
	CurrentMatrix.at(12) = x;
	CurrentMatrix.at(13) = y;
	CurrentMatrix.at(14) = z;
/*	
cerr << "ORIGINAL MV MATRIX  \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << VMCurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
*/	
	if(MatrixMode == 1)
	{
		ProjCurrentMatrix = MatrixMult4x4(ProjCurrentMatrix,CurrentMatrix);
	}
	else if(MatrixMode == 2)
	{
		VMCurrentMatrix = MatrixMult4x4(VMCurrentMatrix,CurrentMatrix);
	}

/*	
	
cerr << "GL TRANS --> MV \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << CurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
	
cerr << "VM * GLtrans  \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << VMCurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
*/	
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{/*
	if(angle < 0)
	{
		angle = 360 + angle;
	}*/
	//NORMALIZE THE X Y Z VALUES
	// sqrt(x^2 + y^2 + z^2) 
	float temp_x, temp_y, temp_z;
	temp_x = x*x;
	temp_y = y*y;
	temp_z = z*z;
	//store x^2 + y^2 + z^2  in temp_x
	temp_x = temp_x + temp_y + temp_z;

	//store mag of x y z in temp_x
	temp_x = sqrt(temp_x);
	if(temp_x != 1)
	{
		x = x/temp_x;
		y = y/temp_x;
		z = z/temp_x;
	}
	
	
	vector<MGLfloat> CurrentMatrix;
	CurrentMatrix.resize(16);
	MGLfloat PI = 3.14159265; 
	MGLfloat c = cos(PI * angle / 180);
	MGLfloat s = sin(PI * angle / 180);
	
	CurrentMatrix[0] = (x*x) * (1-c) + c;
	CurrentMatrix[1] = y*x*(1-c) + z*s;
	CurrentMatrix[2] = x*z*(1-c) - y*s;
	CurrentMatrix[4] = x*y*(1-c) - z*s;
	CurrentMatrix[5] = y*y*(1-c) + c;
	CurrentMatrix[6] = y*z*(1-c) + x*s;
	CurrentMatrix[8] = x*z*(1-c) + y*s;
	CurrentMatrix[9] = y*z*(1-c) - x*s;
	CurrentMatrix[10] = z*z*(1-c) + c;
	CurrentMatrix[15] = 1;
	
	
	
		
	if(MatrixMode == 1)
	{
		ProjCurrentMatrix = MatrixMult4x4(ProjCurrentMatrix,CurrentMatrix);
	}
	else if(MatrixMode == 2)
	{
		VMCurrentMatrix = MatrixMult4x4(VMCurrentMatrix,CurrentMatrix);
	}
	
cerr << "Rotate Matrix " << endl;
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << CurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
cerr << "VMatrix * glrotate  \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << VMCurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}	

}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	vector<MGLfloat> k = {x,0,0,0,0,y,0,0,0,0,z,0,0,0,0,1};
	if(MatrixMode == 1)
	{
		ProjCurrentMatrix = MatrixMult4x4(ProjCurrentMatrix,k);
	}
	else if(MatrixMode == 2)
	{
		VMCurrentMatrix = MatrixMult4x4(VMCurrentMatrix,k);
	}
/*	
cerr << "Scale Matrix " << endl;
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << k[i + j*4] << ' ';
		}
		cerr << endl;
	}
cerr << "VMatrix * GLScale  \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << VMCurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
*/	
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{

	vector<MGLfloat> temp;
	temp.resize(16);
	MGLfloat A = (right + left) / (right - left); 
	MGLfloat B = (top + bottom) / (top - bottom);
	MGLfloat C = (far + near) / (far - near);
	MGLfloat D =  (2 * far * near) / (far - near);
	
	temp[0] = (2 * near) / (right - left);
	temp[5] = (2 * near) / (top - bottom);
	temp[8] = A;
	temp[9] = B;
	temp[10] = -C;
	temp[11] = -1;
	temp[14] = -D;
	
	if(MatrixMode == 1)
	{
		ProjCurrentMatrix = MatrixMult4x4(ProjCurrentMatrix, temp);
	}
	else if(MatrixMode == 2)
	{
		VMCurrentMatrix = MatrixMult4x4(VMCurrentMatrix, temp);
	}
	/*
	cerr << " Projection Matrix \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << ProjCurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}*/
	
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{

	Zfar = far;
	Znear = near; 

	vector<MGLfloat> originalMatrix;
	originalMatrix.resize(16);
		
	originalMatrix[0] = 2/(right - left);
	originalMatrix[1] = 0;
	originalMatrix[2] = 0;
	originalMatrix[3] = 0;
	
	originalMatrix[4] = 0;
	originalMatrix[5] = 2/(top -bottom);
	originalMatrix[6] = 0;
	originalMatrix[7] = 0;
	
	originalMatrix[8] = 0;
	originalMatrix[9] = 0;
	originalMatrix[10] = -2/(far - near);
	originalMatrix[11] = 0;
	
	originalMatrix[12] = -(right + left) / (right - left);
	originalMatrix[13] = -(top + bottom) / (top - bottom);
	if(-near != far)
	{	
		originalMatrix[14] = (near + far) / (near - far);
	}
	else
	{
		originalMatrix[14] = 0;
	}
	originalMatrix[15] = 1;
/*	
cerr << "Original ProjMatrix "<< endl;
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << ProjCurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
	cerr << endl;
*/	
	if(MatrixMode == 1)
	{
		ProjCurrentMatrix =	MatrixMult4x4(ProjCurrentMatrix, originalMatrix);
	}
	else if(MatrixMode == 2)
	{
		VMCurrentMatrix = MatrixMult4x4(VMCurrentMatrix, originalMatrix);
	}
/*	
cerr << " Ortho \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << originalMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}

cerr << " Projection Matrix * Ortho \n";
	for(unsigned i = 0; i < 4; i++)
	{ 
		for(unsigned j = 0; j < 4; j++)
		{
			cerr << ProjCurrentMatrix[i + j*4] << ' ';
		}
		cerr << endl;
	}
*/	
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	GLOBAL_COLOR[0] = red;
	GLOBAL_COLOR[1] = green;
	GLOBAL_COLOR[2] = blue;
}

//4x4 matrix * 4x4 matrix saved in current
vector<MGLfloat> MatrixMult4x4(vector<MGLfloat> i, vector<MGLfloat> j)
{
	vector<MGLfloat> POP;
	POP.resize(16);
	
	POP[0] = i[0]* j[0] + i[4]* j[1] + i[8] * j[2] + i[12] * j[3] ;
	POP[1] = i[1]* j[0] + i[5]* j[1] + i[9] * j[2] + i[13] * j[3];
	POP[2] = i[2]* j[0] + i[6]* j[1] + i[10] * j[2] + i[14] * j[3];
	POP[3] = i[3]* j[0] + i[7]* j[1] + i[11] * j[2] + i[15] * j[3];
	
	POP[4] = i[0]* j[4] + i[4]* j[5] + i[8] * j[6] + i[12] * j[7];
	POP[5] = i[1]* j[4] + i[5]* j[5] + i[9] * j[6] + i[13] * j[7];
	POP[6] = i[2]* j[4] + i[6]* j[5] + i[10] * j[6] + i[14] * j[7];
	POP[7] = i[3]* j[4] + i[7]* j[5] + i[11] * j[6] + i[15] * j[7];
	
	POP[8] = i[0]* j[8] + i[4]* j[9] + i[8] * j[10] + i[12] * j[11];
	POP[9] = i[1]* j[8] + i[5]* j[9] + i[9] * j[10] + i[13] * j[11];
	POP[10] = i[2]* j[8] + i[6]* j[9] + i[10] * j[10] + i[14] * j[11];
	POP[11] = i[3]* j[8] + i[7]* j[9] + i[11] * j[10] + i[15] * j[11];
	
	POP[12] = i[0]* j[12] + i[4]* j[13] + i[8] * j[14] + i[12] * j[15];
	POP[13] = i[1]* j[12] + i[5]* j[13] + i[9] * j[14] + i[13] * j[15];
	POP[14] = i[2]* j[12] + i[6]* j[13] + i[10] * j[14] + i[14] * j[15];
	POP[15] = i[3]* j[12] + i[7]* j[13] + i[11] * j[14] + i[15] * j[15];
	
	return POP;
	
}


//dot product
MGLfloat dot(vector<MGLfloat> X, vector<MGLfloat> Y)
{
	MGLfloat T = 0;
	for(unsigned int i = 0; i < X.size(); i++)
	{
		T = T + (X[i] * Y[i]);
	}
	return T;
	
}


//Current matrix * VERTEX
vector<MGLfloat> MultMatrixCUR(vector<MGLfloat> i, MGLint MODE)
{

	MGLfloat temp0, temp1, temp2, temp3;
	temp0 = 0;
	temp1 = 0;
	temp2 = 0;
	temp3 = 0;
	vector<MGLfloat> CurrentMatrix;
	if(MODE == 1)
	{
		CurrentMatrix = ProjCurrentMatrix;
	}
	else if(MODE == 2)
	{
		CurrentMatrix = VMCurrentMatrix;
	}
	

	temp0 = CurrentMatrix.at(0) * i[0];
	temp0 = temp0 + CurrentMatrix.at(4)*i[1];
	temp0 = temp0 + CurrentMatrix.at(8)*i[2];
	temp0 = temp0 + CurrentMatrix.at(12)*i[3];
	

	temp1 = CurrentMatrix.at(1) * i[0];
	temp1 = temp1 + CurrentMatrix.at(5)*i[1];
	temp1 = temp1 + CurrentMatrix.at(9)*i[2];
	temp1 = temp1 + CurrentMatrix.at(13)*i[3];

	temp2 = CurrentMatrix.at(2) * i[0];
	temp2 = temp2 + CurrentMatrix.at(6)*i[1];
	temp2 = temp2 + CurrentMatrix.at(10)*i[2];
	temp2 = temp2 + CurrentMatrix.at(14)*i[3];
	

	temp3 = CurrentMatrix.at(3) * i[0];
	temp3 = temp3 + CurrentMatrix.at(7)*i[1];
	temp3 = temp3 + CurrentMatrix.at(11)*i[2];
	temp3 = temp3 + CurrentMatrix.at(15)*i[3];
	
	
	i[0] = temp0;
	i[1] = temp1;
	i[2] = temp2;
	i[3] = temp3;
	return i;
	//cerr << temp0 << ' '<< temp1 << ' '<< temp2 << ' '<< temp3 << endl << endl;
	
}


