//
//  CS6555 Computer Animation, Fall 2017
//  Assignment 5
//  Flying Bird Wings Simulator
//  main.cpp
//  SimpleGLUT
//  Created by Mengqi Xie on 2017/12/13.
//  avaqi@gwu.edu
//
#include "stdafx.h"
#include <ctime> 
#include <cstdlib>
// standard
#include <assert.h>
#include <math.h>

// glut
#include <GL/glut.h>


//================================
// global variables
//================================
// screen size
int g_screenWidth = 0;
int g_screenHeight = 0;

// frame index
int g_frameIndex = 0;

// angle for rotation
int g_angle = 0;

//define PI
#define PI 3.1415926

// rotation matrix
GLfloat M[16] = { 1,0,0,0,
0,1,0,0,
0,0,1,0,
0,0,0,1 };

// translation matrix
GLfloat Tl[16] = { 1,0,0,0,
0,1,0,0,
0,0,1,0,
0,0,0,1 };

// translation matrix
GLfloat Tr[16] = { 1,0,0,0,
0,1,0,0,
0,0,1,0,
0,0,0,1 };

// translation matrix
GLfloat T[16] = { 1,0,0,0,
0,1,0,0,
0,0,1,0,
0,0,0,1 };

// translation matrix
GLfloat T_ring[16] = { 1,0,0,0,
0,1,0,0,
0,0,1,0,
0,0,0,1 };

//Euler angels
float x_Phi = 0;  //Phi
float y_Theta = 0; //Theta
float z_Psi = 0; //Psi

//interpolation
//Catmull-Rom spline
float spline_t = 0;
float mov_x = 0, mov_y = 0;
float pos_x = 0, pos_y = 0;
GLfloat posMatx[16] = { 0 };

float px0 = -50, px1 = -20, px2 = 15, px3 = 42;
float py0 = 40, py1 = -20, py2 = 10, py3 = -32;

float p2x0 = -50, p2x1 = 15, p2x2 = -20, p2x3 = 42;
float p2y0 = 40, p2y1 = 10, p2y2 = -20, p2y3 = -32;
bool CRSpflag = true;
bool CRSpswitch = true;
float CRSstep = 0.002;

//method switch
bool isEuler = true;
bool isCRSpline = true;
int stepFlag = 1;

//closed B-spline
GLfloat bspMatx[1200][2] = { 0,0 };
int step = 0;
int speed = 0;

#if 0
// the points of the curve - these are the same as the bezier curve
// points demonstrated in the bezier curve example.

float Points[4][3] = {
	{ 10,10,0 },
	{ 5,10,2 },
	{ -5,0,0 },
	{ -10,5,-2 }
};

#define NUM_POINTS 4

//    The following sets of 4 indices are the curves that need to
//    be drawn to create a clamped cubic b-spline. In total there
//    are 5 curve segments to draw.
//
//    0 0 0 1
//    0 0 1 2
//      0 1 2 3
//        1 2 3 3
//          2 3 3 3
//
//    Remember this when trying to understand knot vectors!!
//

#else

float Points[9][3] = {
	{ 10,5,0 },
	{ 5,10,0 },
	{ -5,15,0 },
	{ -10,-5,0 },
	{ 4,-4,0 },
	{ 10,5,0 },
	{ 5,10,0 },
	{ -5,15,0 },
	{ -10,-5,0 }
};

#define NUM_POINTS 9
#endif

// the level of detail for the curve
unsigned int LOD = 200;

#define NUM_SEGMENTS (NUM_POINTS-3)
//

//================================
// init
//================================
void init(void) {
	// init something before main loop...

}

//================================
// update
//================================
void update(void) {
	// do something before rendering...

	// rotation angle
	//g_angle = ( g_angle + 5 ) % 360;
}

//================================================
// Rotation based on Euler angle and Quaternion
//================================================
void EulerAngleRotate() {
	M[0] = cos(y_Theta)*cos(z_Psi);
	M[1] = cos(z_Psi)*sin(x_Phi)*sin(y_Theta) - cos(x_Phi)*sin(z_Psi);
	M[2] = cos(x_Phi)*cos(z_Psi)*sin(y_Theta) + sin(x_Phi)*sin(z_Psi);

	M[4] = sin(z_Psi)*cos(y_Theta);
	M[5] = sin(x_Phi)*sin(y_Theta)*sin(z_Psi) + cos(x_Phi)*cos(z_Psi);
	M[6] = cos(x_Phi)*sin(y_Theta)*sin(z_Psi) - sin(x_Phi)*cos(z_Psi);

	M[8] = -sin(y_Theta);
	M[9] = sin(x_Phi)*cos(y_Theta);
	M[10] = cos(x_Phi)*cos(y_Theta);
}

void EulertoQuaternion() {
	//GLfloat quaternion[4] = { 0 };
	float w = cos(x_Phi / 2)*cos(y_Theta / 2)*cos(z_Psi / 2) + sin(x_Phi / 2)*sin(y_Theta / 2)*sin(z_Psi / 2);
	float x = sin(x_Phi / 2)*cos(y_Theta / 2)*cos(z_Psi / 2) - cos(x_Phi / 2)*sin(y_Theta / 2)*sin(z_Psi / 2);
	float y = cos(x_Phi / 2)*sin(y_Theta / 2)*cos(z_Psi / 2) + sin(x_Phi / 2)*cos(y_Theta / 2)*sin(z_Psi / 2);
	float z = cos(x_Phi / 2)*cos(y_Theta / 2)*sin(z_Psi / 2) - sin(x_Phi / 2)*sin(y_Theta / 2)*cos(z_Psi / 2);

	M[0] = w * w + x * x - y * y - z * z;
	M[1] = 2 * x * y - 2 * w * z;
	M[2] = 2 * x * z + 2 * w * y;

	M[4] = 2 * x * y + 2 * w * z;
	M[5] = w * w - x * x + y * y - z * z;
	M[6] = 2 * y * z - 2 * w * x;

	M[8] = 2 * x * z - 2 * w * y;
	M[9] = 2 * y * z + 2 * w * x;
	M[10] = w * w - x * x - y * y + z * z;
}

//=======================================================
// interpolation done by Catmull-Rom spline and B spline
//=======================================================

float CatmullRom(float t, float x0, float x1, float x2, float x3)
{
	//caculate in between frames by control points
	float t2 = t * t;
	float t3 = t2 * t;
	return ((2 * x1) +
		(-x0 + x2) * t +
		(2 * x0 - 5 * x1 + 4 * x2 - x3) * t2 +
		(-x0 + 3 * x1 - 3 * x2 + x3) * t3) * 0.5f;
}

void getCRSpline(float pos_x, float pos_y) {
	//do interpolation
	if (CRSpswitch) {
		mov_x = CatmullRom(spline_t, px0, px1, px2, px3) - pos_x;
		mov_y = CatmullRom(spline_t, py0, py1, py2, py3) - pos_y;
	}
	else {
		mov_x = CatmullRom(spline_t, p2x0, p2x1, p2x2, p2x3) - pos_x;
		mov_y = CatmullRom(spline_t, p2y0, p2y1, p2y2, p2y3) - pos_y;
	}
}


float* GetPoint(int i)
{
	// return 1st point
	if (i<0)
	{
		return    Points[0];
	}

	if (i<NUM_POINTS)
	{
		return Points[i];
	}
	// return last point

	return Points[NUM_POINTS - 1];
}
void bSpline()
{
	int index = 0;
	for (int start_cv = 0, j = 0; j<NUM_SEGMENTS; j++, start_cv++)
	{
		// for each section of curve, draw LOD number of divisions
		for (int i = 0; i != LOD; ++i)
		{
			// use the parametric time value 0 to 1 for this curve
			// segment.
			float t = (float)i / LOD;
			// the t value inverted
			float it = 1.0f - t;

			// calculate blending functions for cubic bspline
			float b0 = it*it*it / 6.0f;
			float b1 = (3 * t*t*t - 6 * t*t + 4) / 6.0f;
			float b2 = (-3 * t*t*t + 3 * t*t + 3 * t + 1) / 6.0f;
			float b3 = t*t*t / 6.0f;

			// calculate the x,y and z of the curve point
			float x = b0 * GetPoint(start_cv + 0)[0] +
				b1 * GetPoint(start_cv + 1)[0] +
				b2 * GetPoint(start_cv + 2)[0] +
				b3 * GetPoint(start_cv + 3)[0];

			float y = b0 * GetPoint(start_cv + 0)[1] +
				b1 * GetPoint(start_cv + 1)[1] +
				b2 * GetPoint(start_cv + 2)[1] +
				b3 * GetPoint(start_cv + 3)[1];

			float z = b0 * GetPoint(start_cv + 0)[2] +
				b1 * GetPoint(start_cv + 1)[2] +
				b2 * GetPoint(start_cv + 2)[2] +
				b3 * GetPoint(start_cv + 3)[2];

			// specify the point
			bspMatx[index][0] = x;
			bspMatx[index][1] = y;
			index++;
		}
	}
}

void getBSpline(float pos_x, float pos_y) {
	mov_x = bspMatx[step][0] - pos_x;
	mov_y = bspMatx[step][1] - pos_y;
}

void applyRotation() {
	if (isEuler) {
		EulerAngleRotate();
		glMultMatrixf(M);
	}
	else {
		EulertoQuaternion();
		glMultMatrixf(M);
	}
}

void applyMove(float pos_x, float pos_y) {
	if (isCRSpline) {
		getCRSpline(pos_x, pos_y);
	}
	else {
		getBSpline(pos_x, pos_y);
	}
	glTranslatef(mov_x, mov_y, 0.0);
}

void drawSphere(GLfloat xx, GLfloat yy, GLfloat zz, GLfloat radius, GLfloat M, GLfloat N, float cr, float cg, float cb)
{
	float step_z = PI / M;
	float step_xy = 2 * PI / N;
	float x[4], y[4], z[4];

	float angle_z = 0.0;
	float angle_xy = 0.0;
	int i = 0, j = 0;

	glBegin(GL_POLYGON);
	for (i = 0; i<M; i++)
	{
		angle_z = i * step_z;

		for (j = 0; j<N; j++)
		{
			angle_xy = j * step_xy;

			x[0] = radius * sin(angle_z) * cos(angle_xy);
			y[0] = radius * sin(angle_z) * sin(angle_xy);
			z[0] = radius * cos(angle_z);
			glColor3f(0.0, 1.0, 0.0);

			x[1] = radius * sin(angle_z + step_z) * cos(angle_xy);
			y[1] = radius * sin(angle_z + step_z) * sin(angle_xy);
			z[1] = radius * cos(angle_z + step_z);
			glColor3f(1.0, 0.0, 0.0);

			x[2] = radius*sin(angle_z + step_z)*cos(angle_xy + step_xy);
			y[2] = radius*sin(angle_z + step_z)*sin(angle_xy + step_xy);
			z[2] = radius*cos(angle_z + step_z);
			glColor3f(0.0, 0.0, 1.0);

			x[3] = radius * sin(angle_z) * cos(angle_xy + step_xy);
			y[3] = radius * sin(angle_z) * sin(angle_xy + step_xy);
			z[3] = radius * cos(angle_z);


			for (int k = 0; k<4; k++)
			{
				glColor3f(cr, cg, cb);
				glNormal3f(xx + x[k], yy + y[k], zz + z[k]);
				glTexCoord2f(0.0, k / 4);
				glVertex3f(xx + x[k], yy + y[k], zz + z[k]);
			}
		}
	}
	glEnd();
}

GLfloat centerVertices[][3] = { { -1, -1, -1 },
{ 1, -1, -1 },{ 1, 1, -1 },{ -1, 1, -1 },
{ -1, -1, 1 },{ 1, -1, 1 },{ 1, 1, 1 },
{ -1, 1, 1 } };

float colR, colG, colB;
void polygon(int a,int b, int c, int d, GLfloat vertices[][3]) {
	glColor3f(colR, colG, colB);
	glBegin(GL_POLYGON);
	//v means the parameter is a  array.

	glNormal3fv(vertices[a]);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(vertices[a]);
	 

	glNormal3fv(vertices[b]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(vertices[b]);


	glNormal3fv(vertices[c]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(vertices[c]);


	glNormal3fv(vertices[d]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(vertices[d]);


	glEnd();


}

//draw cube
void drawCube(GLfloat vertices[][3]) {
	polygon(0, 3, 2, 1, vertices);
	polygon(2, 3, 7, 6, vertices);
	polygon(0, 4, 7, 3, vertices);
	polygon(1, 2, 6, 5, vertices);
	polygon(4, 5, 6, 7, vertices);
	polygon(0, 1, 5, 4, vertices);
}

void drawBody() {
	GLfloat larmVertices[8][3] = { 0 };
	GLfloat rarmVertices[8][3] = { 0 };
	GLfloat llegVertices[8][3] = { 0 };
	GLfloat rlegVertices[8][3] = { 0 };

	for (int i = 0; i < 8; i++) {
		if (centerVertices[i][0] > 0) {
			larmVertices[i][0] = centerVertices[i][0] - 2.5;
			llegVertices[i][0] = centerVertices[i][0] - 1.5;
			rarmVertices[i][0] = centerVertices[i][0] + 1;
			rlegVertices[i][0] = centerVertices[i][0];
		}
		if (centerVertices[i][0] < 0) {
			larmVertices[i][0] = centerVertices[i][0] - 1;
			llegVertices[i][0] = centerVertices[i][0];
			rarmVertices[i][0] = centerVertices[i][0] + 2.5;
			rlegVertices[i][0] = centerVertices[i][0] + 1.5;
		}
		if (centerVertices[i][1] < 0) {
			larmVertices[i][1] = centerVertices[i][1] - 1;
			rarmVertices[i][1] = centerVertices[i][1] - 1;
		}
		if (centerVertices[i][1] > 0) {
			larmVertices[i][1] = centerVertices[i][1];
			rarmVertices[i][1] = centerVertices[i][1];
		}

		llegVertices[i][1] = centerVertices[i][1] - 2.5;
		rlegVertices[i][1] = centerVertices[i][1] - 2.5;

		larmVertices[i][2] = centerVertices[i][2];
		rarmVertices[i][2] = centerVertices[i][2];
		llegVertices[i][2] = centerVertices[i][2];
		rlegVertices[i][2] = centerVertices[i][2];
	}
	EulertoQuaternion();
	glPushMatrix();
	glMultMatrixf(M);
	drawCube(larmVertices);
	drawCube(rlegVertices);
	glPopMatrix();


	x_Phi = -x_Phi;
	EulertoQuaternion();
	glPushMatrix();
	glMultMatrixf(M);
	drawCube(rarmVertices);
	drawCube(llegVertices);
	glPopMatrix();
	
	x_Phi = -x_Phi;
}

bool hitground = false;
float v_up = -0.08;
int hit_step = 1;

float v_hr = -0.1;
float a_vr = 0.001;
float v_vr = 0;

bool detecthit(float radius,float vx) {
	float ver_d = 0;
	ver_d = 10 - radius + vx;
	//cen_d = -16 + Tl[12]Tl[12];
	if (ver_d <= -10) {
		return true;
	}
	else return false;
}
void applyPhysics_R(int mystep) {
	if (mystep == 1)
		Tr[12] += v_hr;
	else
		Tr[12] += -v_hr;
	if (detecthit(1, Tr[13])) {
		v_vr = -v_vr + 0.015;
		if (v_vr >= 0)
			v_vr = 0;
		v_hr = -0.01;
	}
	v_vr += a_vr;
	Tr[13] -= v_vr;
	glMultMatrixf(Tr);
}

float v_hl = 0.1;
float a_vl = 0.001;
float v_vl = 0;
void applyPhysics_L(int mystep) {
	if(mystep == 1)
		Tl[12] += v_hl;
	else if(mystep == 2)
		Tl[12] += -v_hl;
	
	if (detecthit(1, Tl[13])) {
		v_vl = -v_vl + 0.02;
		if (v_vl >= 0)
			v_vl = 0;
		v_hl = 0.01;
	}
	v_vl += a_vl;
	Tl[13] -= v_vl;
	glMultMatrixf(Tl);
}

bool collided = false;
void detectCollision() {
	
	float cen_d = 0;
	cen_d = sqrt(pow((-16+Tl[12] - 16-Tr[12]), 2) + pow((Tr[13] - Tl[13]), 2));
	//cen_d = -16 + Tl[12]Tl[12];
	printf("cend = %f \n",cen_d);
	if (cen_d<=1.8) {
		collided = true;
	}
}


float cen_v_v = -0.02;
float cen_h_v = 0;
float cen_h_a = 0.005;

float rotateQuan = PI / 90;
float rotateQuan_c = PI / 90;
int rostep = 0;
void centerMove() {
	if(cen_h_v >=0.3|| cen_h_v <= -0.3)
		cen_h_a = -cen_h_a;
	cen_h_v += cen_h_a;
	T[12] += cen_h_v;
	T[13] += cen_v_v;
	if (T[13] <= -25 || T[13] >= 4)
		cen_v_v = -cen_v_v;
	glMultMatrixf(T);
}

GLfloat cubeVertices[][3] = { { -0.5, -0.5, -0.5 },
{ 0.5, -0.5, -0.5 },{ 0.5, 0.5, -0.5 },{ -0.5, 0.5, -0.5 },
{ -0.5, -0.5, 0.5 },{ 0.5, -0.5, 0.5 },{ 0.5, 0.5, 0.5 },
{ -0.5,0.5, 0.5 } };


void ringrun(float posx,float myring_ta, float myring_r, float rox, float roy, float roz) {
	for (int i = 0; i<4; i++) {
		myring_ta += 2 * PI / 4;
		T_ring[12] = posx;
		T_ring[13] = myring_r * sin(myring_ta);
		T_ring[14] = myring_r * cos(myring_ta);
		glMultMatrixf(T_ring);
		//drawSphere(0.0, 0.0, 0.0, 0.5f, 100, 100, 0.4, 0.1, 0.4);
		glPushMatrix();
		srand(time(NULL));
		x_Phi = rox + rostep * rotateQuan;
		y_Theta = roy + rostep * rotateQuan;
		z_Psi = roz + rostep * rotateQuan;
		EulerAngleRotate();
		glMultMatrixf(M);
		drawCube(cubeVertices);
		glPopMatrix();
		T_ring[12] = -posx;
		T_ring[13] = -myring_r * sin(myring_ta);
		T_ring[14] = -myring_r * cos(myring_ta);
		glMultMatrixf(T_ring);
	}
}



float ring_ta = PI / 4;
float rota_de = PI / 90;

float ring_r_1 = 6;
float ring_r_2 = 4;
float ring_r_3 = 2.1;
float ring_r_4 = 4;
float ring_r_5 = 6;

float ring_rd_1 = 0.1;
float ring_rd_2 = 0.1;
float ring_rd_3 = 0.1;
float ring_rd_4 = 0.1;
float ring_rd_5 = 0.1;

float ring_x_1 = -4;
float ring_x_2 = -2;
float ring_x_3 = 0;
float ring_x_4 = 2;
float ring_x_5 = 4;

float ring_xd_1 = -0.05;
float ring_xd_2 = -0.05;
float ring_xd_3 = 0;
float ring_xd_4 = 0.05;
float ring_xd_5 = 0.05;
void circleObj() {
	rostep++;
	ring_ta += rota_de;
	if (ring_r_1 >= 7 || ring_r_1 <= 2)
		ring_rd_1 = -ring_rd_1;
	ring_r_1 += ring_rd_1;

	if (ring_r_2 >= 7 || ring_r_2 <= 2)
		ring_rd_2 = -ring_rd_2;
	ring_r_2 += ring_rd_2;

	if (ring_r_3 >= 7 || ring_r_3 <= 2)
		ring_rd_3 = -ring_rd_3;
	ring_r_3 += ring_rd_3;

	if (ring_r_4 >= 7 || ring_r_4 <= 2)
		ring_rd_4 = -ring_rd_4;
	ring_r_4 += ring_rd_4;

	if (ring_r_5 >= 7 || ring_r_5 <= 2)
		ring_rd_5 = -ring_rd_5;
	ring_r_5 += ring_rd_5;

	if (ring_x_1 >-4 || ring_x_1 < -7)
		ring_xd_1 = -ring_xd_1;
	ring_x_1 += ring_xd_1;

	if (ring_x_2 >-2 || ring_x_2 < -5)
		ring_xd_2 = -ring_xd_2;
	ring_x_2 += ring_xd_2;

	if (ring_x_5 <4 || ring_x_5 > 7)
		ring_xd_5 = -ring_xd_5;
	ring_x_5 += ring_xd_5;

	if (ring_x_4 <2 || ring_x_4 > 5)
		ring_xd_4 = -ring_xd_4;
	ring_x_4 += ring_xd_4;


	float myring_ta = ring_ta;
	ringrun(ring_x_1, myring_ta, ring_r_1,PI/2,PI/3,PI/5);
	ringrun(ring_x_2, myring_ta + PI / 5, ring_r_2,PI/5,PI/2,PI/3);
	ringrun(ring_x_3, myring_ta + PI*2 / 3, ring_r_3,PI/10,PI/4,PI/2);
	ringrun(ring_x_4, myring_ta + PI*3 / 5, ring_r_4,PI/7,PI/5,PI/6);
	ringrun(ring_x_5, myring_ta + PI*2 / 7, ring_r_5,PI/1,PI/10,PI/8);
	
}

GLfloat R[16] = { 1,0,0,0,
0,1,0,0,
0,0,1,0,
0,0,0,1 };
void EulerRotateMat(float angX, float angY, float angZ) {
	R[0] = cos(angY)*cos(angZ);
	R[1] = cos(angZ)*sin(angX)*sin(angY) - cos(angX)*sin(angZ);
	R[2] = cos(angX)*cos(angZ)*sin(angY) + sin(angX)*sin(angZ);

	R[4] = sin(angZ)*cos(angY);
	R[5] = sin(angX)*sin(angY)*sin(angZ) + cos(angX)*cos(angZ);
	R[6] = cos(angX)*sin(angY)*sin(angZ) - sin(angX)*cos(angZ);

	R[8] = -sin(angY);
	R[9] = sin(angX)*cos(angY);
	R[10] = cos(angX)*cos(angY);
}

// translation matrix
GLfloat Wingt[16] = { 1,0,0,0,
0,1,0,0,
0,0,1,0,
0,0,0,1 };

float leftWing[8][5] = { 0 };
float moveAng = PI / 90;
void drawLWing() {
	if (leftWing[0][1] == 0)
		leftWing[0][1] = 1;
	for (int i = 0; i < 8; i++) {
		if (leftWing[i][2] == 0) {
			leftWing[i][2] = PI / 4 + i*PI / 40;
		}
		glPushMatrix();
		Wingt[12] = -0.7;
		Wingt[13] = 0.5;
		glMultMatrixf(Wingt);
		if (leftWing[i][1] == 1) {
			if (leftWing[i][0] >= leftWing[i][2]) {
				leftWing[i][1] = 2;
				leftWing[i][4] = 0;
			}
		}
		if (leftWing[i][1] == -1) {
			if (leftWing[i][0] <= -leftWing[i][2] / 2)
			{
				leftWing[i][1] = -2;
				leftWing[i][4] = 0;
			}
		}
		if (leftWing[i][1] == 1) {
			leftWing[i][0] += moveAng;
			leftWing[i][4]++;
		}
		if (leftWing[i][1] == -1) {
			leftWing[i][0] -= moveAng;
			leftWing[i][4]++;
		}
		if (i != 7) {
			if (leftWing[i][4] >= 5)
				leftWing[i + 1][1] = leftWing[i][1];
		}
		EulerRotateMat(0, 0, leftWing[i][0]);
		glMultMatrixf(R);
		Wingt[12] = (i + 1)*(-1.5);
		Wingt[13] = (i + 1)*(-0.5);
		glMultMatrixf(Wingt);
		drawCube(cubeVertices);
		glPopMatrix();
	}
	if (leftWing[4][1] == 2)
		leftWing[0][1] = -1;
	if (leftWing[4][1] == -2)
		leftWing[0][1] = 1;
}

float rightWing[8][5] = { 0 };
float moveRAng = -PI / 90;
void drawRWing() {
	if (rightWing[0][1] == 0)
		rightWing[0][1] = 1;
	for (int i = 0; i < 8; i++) {
		if (rightWing[i][2] == 0) {
			rightWing[i][2] = (PI / 4 + i*PI / 40);
		}
		glPushMatrix();
		Wingt[12] = 0.7;
		Wingt[13] = 0.5;
		glMultMatrixf(Wingt);
		if (rightWing[i][1] == 1) {
			if (rightWing[i][0] <= -rightWing[i][2]) {
				rightWing[i][1] = 2;
				rightWing[i][4] = 0;
			}
		}
		if (rightWing[i][1] == -1) {
			if (rightWing[i][0] >= rightWing[i][2] / 2)
			{
				rightWing[i][1] = -2;
				rightWing[i][4] = 0;
			}
		}
		if (rightWing[i][1] == 1) {
			rightWing[i][0] += moveRAng;
			rightWing[i][4]++;
		}
		if (rightWing[i][1] == -1) {
			rightWing[i][0] -= moveRAng;
			rightWing[i][4]++;
		}
		if (i != 7) {
			if (rightWing[i][4] >= 5)
				rightWing[i + 1][1] = rightWing[i][1];
		}
		EulerRotateMat(0, 0, rightWing[i][0]);
		glMultMatrixf(R);
		Wingt[12] = (i + 1)*(1.5);
		Wingt[13] = (i + 1)*(-0.5);
		glMultMatrixf(Wingt);
		drawCube(cubeVertices);
		glPopMatrix();
	}
	if (rightWing[4][1] == 2)
		rightWing[0][1] = -1;
	if (rightWing[4][1] == -2)
		rightWing[0][1] = 1;
}
//================================
// render
//================================
void render(void) {
	// clear buffer
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClearDepth(1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// render state
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);

	// enable lighting
	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);

	// light source attributes
	GLfloat LightAmbient[] = { 0.4f, 0.2f, 0.4f, 1.0f };
	GLfloat LightDiffuse[] = { 0.3f, 0.6f, 0.3f, 1.0f };
	GLfloat LightSpecular[] = { 0.1f, 0.4f, 0.4f, 1.0f };
	GLfloat LightPosition[] = { -5.0f, 5.0f, 5.0f, 1.0f };

	glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, LightSpecular);
	glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);

	// surface material attributes
	GLfloat material_Ka[] = { 0.71f, 0.16f, 0.11f, 1.0f };
	GLfloat material_Kd[] = { 0.13f, 0.6f, 0.24f, 1.0f };
	GLfloat material_Ks[] = { 0.33f, 0.13f, 0.72f, 1.0f };
	GLfloat material_Ke[] = { 0.1f , 0.0f , 0.1f , 1.0f };
	GLfloat material_Se = 10;

	glMaterialfv(GL_FRONT, GL_AMBIENT, material_Ka);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, material_Kd);
	glMaterialfv(GL_FRONT, GL_SPECULAR, material_Ks);
	glMaterialfv(GL_FRONT, GL_EMISSION, material_Ke);
	glMaterialf(GL_FRONT, GL_SHININESS, material_Se);

	// modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(
		0, 0, 0.5,
		0, 0, 0,
		0.0f, 1.0f, 0.0f);

	//glTranslatef(0.0, 0.0, -5.0);
	
	//perform rotation and interpolation
	
	//applyRotation();
	//if (!collided)
	//detectCollision();

	// render objects
	
	glTranslatef(-14.0, 14.0, -5.0);
	centerMove();
	//glTranslatef(0.0, 1.0, 0.0);
	EulerAngleRotate();
	glMultMatrixf(M);
	//drawSphere(0.0, 0.0, 0.0, 1.0f, 100, 100, 0.2, 0.1, 0.5);
	colR = 0.5;
	colG = 0.3;
	colB = 0.9;
	drawCube(centerVertices);
	colR = 0.2;
	colG = 0.2;
	colB = 0.3;
	drawLWing();
	drawRWing();
	glTranslatef(0.0, 0.0, -2.0);
	colR = 0.4;
	colG = 0.9;
	colB = 0.9;
	drawLWing();
	drawRWing();
	glTranslatef(0.0, 0.0, 4.0);
	colR = 0.9;
	colG = 0.4;
	colB = 0.4;
	drawLWing();
	drawRWing();
	//glLoadIdentity();
	//drawCube(cubeVertices);
	/*
	glTranslatef(-16.0, 10.0, -5.0);
	if(!collided)
		applyPhysics_L(1);
	else
		applyPhysics_L(2);
	drawSphere(0.0, 0.0, 0.0, 1.0f, 100, 100);
	glLoadIdentity();
	glTranslatef(16.0, 10.0, -5.0);
	if (!collided)
		applyPhysics_R(1);
	else
		applyPhysics_R(2);
	drawSphere(0.0, 0.0, 0.0, 0.8f, 100, 100);
	*/
	
	//glGetFloatv(GL_MODELVIEW_MATRIX, posMatx);
	pos_x = posMatx[12];
	pos_y = posMatx[13];
	// disable lighting
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);

	// swap back and front buffers
	glutSwapBuffers();
}

//================================
// keyboard input
//================================
void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 'e':
		isEuler = true;
		break;
	case 'q':
		isEuler = false;
		break;
	case 'b':
		isCRSpline = false;
		break;
	case 'c':
		isCRSpline = true;
		break;
	case 'w':
		x_Phi += PI/180;
		break;
	case 's':
		x_Phi -= PI / 180;

	case 'a':
		y_Theta += PI / 180;
		break;
	case 'd':
		y_Theta -= PI / 180;
		break;
	default:
		break;
	}
}

//================================
// reshape : update viewport and projection matrix when the window is resized
//================================
void reshape(int w, int h) {
	// screen size
	g_screenWidth = w;
	g_screenHeight = h;

	// viewport
	//glViewport( 0, 0, (GLsizei)w, (GLsizei)h );

	// projection matrix
	//glMatrixMode( GL_PROJECTION );


	// prevents division by zero when minimising window
	if (h == 0)
	{
		h = 1;
	}

	// set the drawable region of the window
	glViewport(0, 0, w, h);
	// set up the projection matrix
	glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	// just use a perspective projection
	//gluPerspective(45,(float)w/h,0.1,100);
	if (w <= h)
	{
		glOrtho(-20.0, 20.0, -20.0*(GLfloat)h / (GLfloat)w, 20.0*(GLfloat)h / (GLfloat)w, 0.0, 100.0);
	}
	else
	{
		glOrtho(-20.0, 20.0, -20.0*(GLfloat)h / (GLfloat)w, 20.0*(GLfloat)h / (GLfloat)w, 0.0, 100.0);
	}
	// go back to modelview matrix so we can move the objects about



	glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();

	glLoadIdentity();
	//gluPerspective(45.0, (GLfloat)w/(GLfloat)h, 1.0, 2000.0);
}

int degreestep = 0;
//================================
// timer : triggered every 16ms ( about 60 frames per second )
//================================
void timer(int value) {
	// increase frame index
	//g_frameIndex++;
	/*
	if (spline_t<1)
		spline_t += CRSstep;
	else {
		CRSpswitch = !CRSpswitch;
		spline_t = 0;
	}
	speed++;
	if (speed % 1 == 0) {
		if (step<1198)
			step += 2;
		else
			step = 200;
	}
	
	if (x_Phi>= PI / 4)
		stepFlag = -stepFlag;
	if (x_Phi <= -PI / 4)
		stepFlag = -stepFlag;
	if (stepFlag==1)
		x_Phi += PI / 180;
	else
		x_Phi -= PI / 180;
	y_Theta += PI / 180;
	z_Psi += PI / 180;
	
	
	if(step%2==0)
	timestep++;
	step++;
*/

	//y_Theta += PI / 180;
	//z_Psi += PI / 180;
	//x_Phi += PI / 180;
	update();

	// render
	glutPostRedisplay();

	// reset timer
	// 16 ms per frame ( about 60 frames per second )
	glutTimerFunc(32, timer, 0);
}

//================================
// main
//================================
int main(int argc, char** argv) {
	// create opengL window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(argv[0]);

	//bSpline();
	// init
	init();

	// set callback functions
	glutDisplayFunc(render);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutTimerFunc(16, timer, 0);

	// main loop
	glutMainLoop();

	return 0;
}