#include <Windows.h>
#include <GL/glut.h>
#include <math.h>
#include "ScanConvert.h"
#include "PolygonDrawer.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <climits>
#include <string>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* strtof */
#include <sstream> 

using namespace std;

/******************************************************************
	Notes:
	Image size is 400 by 400 by default.  You may adjust this if
		you want to.
	You can assume the window will NOT be resized.
	Call clearFramebuffer to clear the entire framebuffer.
	Call setFramebuffer to set a pixel.  This should be the only
		routine you use to set the color (other than clearing the
		entire framebuffer).  drawit() will cause the current
		framebuffer to be displayed.
	As is, your scan conversion should probably be called from
		within the display function.  There is a very short sample
		of code there now.
	You may add code to any of the subroutines here,  You probably
		want to leave the drawit, clearFramebuffer, and
		setFramebuffer commands alone, though.
  *****************************************************************/

float framebuffer[ImageH][ImageW][3];
float zbuffer[ImageH][ImageW];

vector<Face> faces;
vector<Pt> pts;
bool paintshape = true;
string filename;
int shading = 1;

void keyboard(unsigned char key, int x, int y)
{
	paintshape = true;

	switch (key)
	{
	case '1':
		shading = 1;
		//clearZBuffer();
		//Flat();
		break;
	case '2':
		shading = 2;
		//clearZBuffer();
		//Gouraud();
		break;
	case '3':
		shading = 3;
		//clearZBuffer();
		//Phong();
		break;
	}
	glutPostRedisplay();
}

// Draws the scene
void drawit(void)
{
   glDrawPixels(ImageW,ImageH,GL_RGB,GL_FLOAT,framebuffer);
   glFlush();
}

// Clears framebuffer to black
void clearFramebuffer()
{
	int i,j;

	for(i=0;i<ImageH;i++) {
		for (j=0;j<ImageW;j++) {
			framebuffer[i][j][0] = 0.0;
			framebuffer[i][j][1] = 0.0;
			framebuffer[i][j][2] = 0.0;
			zbuffer[i][j]= INT_MIN;
		}
	}
}

// Clears zBuffer to min distance
void clearZBuffer()
{
	int i, j;
	for (i = 0; i<ImageH; i++) {
		for (j = 0; j<ImageW; j++) {
			zbuffer[i][j] = INT_MIN;
		}
	}
}

// Sets pixel x,y to the color RGB
// I've made a small change to this function to make the pixels match
// those returned by the glutMouseFunc exactly - Scott Schaefer 
void setFramebuffer(int x, int y, int z, float R, float G, float B)
{
	// changes the origin from the lower-left corner to the upper-left corner
	y = ImageH - 1 - y;
	if (zbuffer[y][x] < z)
	{
		zbuffer[y][x] = z;
		if (R <= 1.0)
			if (R >= 0.0)
				framebuffer[y][x][0] = R;
			else
				framebuffer[y][x][0] = 0.0;
		else
			framebuffer[y][x][0] = 1.0;
		if (G <= 1.0)
			if (G >= 0.0)
				framebuffer[y][x][1] = G;
			else
				framebuffer[y][x][1] = 0.0;
		else
			framebuffer[y][x][1] = 1.0;
		if (B <= 1.0)
			if (B >= 0.0)
				framebuffer[y][x][2] = B;
			else
				framebuffer[y][x][2] = 0.0;
		else
			framebuffer[y][x][2] = 1.0;
	}
}

void Flat(void)
{
	if (paintshape)
	{
		for (int i = 0; i < faces.size(); i++)
		{
			if (faces[i].normal.z < 0) //back face culling
			{
				float r, g, b;

				double Lx = 1 / sqrt(3);
				double Ly = 1 / sqrt(3);
				double Lz = -1 / sqrt(3);

				double LN = (faces[i].normal.x*Lx + faces[i].normal.y*Ly + faces[i].normal.z*Lz);

				double Rx = 2 * LN*faces[i].normal.x - Lx;
				double Ry = 2 * LN*faces[i].normal.y - Ly;
				double Rz = 2 * LN*faces[i].normal.z - Lz;
				double RE = (Rx * 0 + Ry * 0 + Rz*-1);
				int n = 5;

				r = .5*.1;
				g = .5 * 0;
				b = .5 * 0;
				if (LN > 0)
				{
					r += (1) * (.7*LN + .5*pow(RE, n));
					g += (1) * (0 * LN + .5*pow(RE, n));
					b += (1) * (0 * LN + .5*pow(RE, n));
				}

				color c;
				c.r = r; c.g = g; c.b = b;
				vector<Pt> pts;
				pts.push_back(faces[i].point1.screen());
				pts.push_back(faces[i].point2.screen());
				pts.push_back(faces[i].point3.screen());
				drawPolygon(pts, c);
				drawit();
			}
			//cout << "points: " << pts[0].x << " " << pts[1].y << " " << pts[2].x << endl;
		}
	}
}

void setVertexNormals(void)
{
	for (int i = 0; i < pts.size(); i++)
	{
		for (int j = 0; j < pts[i].faces_included.size(); j++)
		{
			pts[i].nx += faces[pts[i].faces_included[j]].normal.x;
			pts[i].ny += faces[pts[i].faces_included[j]].normal.y;
			pts[i].nz += faces[pts[i].faces_included[j]].normal.z;
		}
		pts[i].normalize();
	}
}

void Gouraud(void)
{
	setVertexNormals();
	if (paintshape)
	{
		for (int i = 0; i < faces.size(); i++)
		{
			if (faces[i].normal.z <= 0) //back face culling
			{
				vector<Pt> points;
				points.push_back(pts[faces[i].point1.index].screen());
				points.push_back(pts[faces[i].point2.index].screen());
				points.push_back(pts[faces[i].point3.index].screen());
				drawGouraudPolygon(points);
				drawit();
			}
			//cout << "points: " << pts[0].x << " " << pts[1].y << " " << pts[2].x << endl;
		}
	}
}

void Phong(void)
{
	setVertexNormals();
	if (paintshape)
	{
		for (int i = 0; i < faces.size(); i++)
		{
			if (faces[i].normal.z <= 0) //back face culling
			{
				vector<Pt> points;
				points.push_back(pts[faces[i].point1.index].screen());
				points.push_back(pts[faces[i].point2.index].screen());
				points.push_back(pts[faces[i].point3.index].screen());
				drawPhongPolygon(points);
				drawit();
			}
			//cout << "points: " << pts[0].x << " " << pts[1].y << " " << pts[2].x << endl;
		}
	}
}

void display(void)
{
	if (shading == 1)
	{
		clearZBuffer();
		Flat();
	}

	if (shading == 2)
	{
		clearZBuffer();
		Gouraud();
	}
	
	if (shading == 3)
	{
		clearZBuffer();
		Phong();
	}
	
	paintshape = false;
}

void readInput(void)
{

	//filename = "cow.obj";
	string input;
	ifstream file(filename);				//input file name
	bool skip;
	int count = 0;
	while (!file.eof())						//runs until the end of the file is reached
	{						//three variables used to fill in the data members of the student
		getline(file, input);						//attempts to fill in the first variable
		if (file.fail())					//if the input has an error it lets the user know that there is an error and breaks out of the while loop 
		{
			//cout<<"\n There was an error in the input, data not read in, will exit \n\"ERROR\"\n\n";
			break;
		}
		stringstream ss(input);

		if (input[0] == 'v')
		{
			string temp, point1, point2, point3;
			ss >> temp >>point1 >> point2 >> point3;
			Pt point(atof(point1.c_str()), atof(point2.c_str()), atof(point3.c_str()), pts.size());
			//cout << input<< endl;
			//cout << "test: "<< atof(point1.c_str()) <<" "<< point1 << endl;
			pts.push_back(point);
		}
		else if (input[0] == 'f')
		{
			string temp, point1, point2, point3;
			ss >> temp >> point1 >> point2 >> point3;
			//cout << input << endl;
			//cout << "New Face: " << point1 << " " << point2 << " " << point3 << " " << endl;
			Face f(pts[stoi(point1)-1], pts[stoi(point2)-1], pts[stoi(point3)-1]);
			faces.push_back(f);
			pts[stoi(point1) - 1].add(count);
			pts[stoi(point2) - 1].add(count);
			pts[stoi(point3) - 1].add(count);
			count++;
		}
	}
	file.close();			//closes the file
}

void init(void) {
	clearFramebuffer();
	readInput();
	display();
}




int main(int argc, char* argv[])
{
	filename = argv[1];

	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);
	glutInitWindowSize(ImageW,ImageH);
	glutInitWindowPosition(100,100);
	glutCreateWindow("Inaki Rosa - Homework 4");
	init();
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}
