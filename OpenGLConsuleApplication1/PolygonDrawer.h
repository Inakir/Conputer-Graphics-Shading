#ifndef POLYGON_DRAWER_H
#define POLYGON_DRAWER_H

#include <vector>
#include "ScanConvert.h"

using namespace std;

class Pt
{
public:
	float x, y, z;
	float r, g, b;
	float nx, ny, nz;
	int index;
	vector<int> faces_included;
	Pt ( void )
	{
		x = y = z = r = g = b = nx = ny = nz = 0;
	}

	Pt ( float nX, float nY, float nZ, int i)
	{
		x = nX;
		y = nY;
		z = nZ;
		index = i;
		r = g = b = nx = ny = nz = 0;
	}

	Pt(float nX, float nY, float nZ, float nR, float nG, float nB, int i)
	{
		x = nX;
		y = nY;
		z = nZ;

		r = nR;
		g = nG;
		b = nB;

		index = i;
		nx = ny = nz = 0;
	}

	Pt(float nX, float nY, float nZ, float nR, float nG, float nB, float xn, float yn, float zn, int i)
	{
		x = nX;
		y = nY;
		z = nZ;

		r = nR;
		g = nG;
		b = nB;

		index = i;
		nx = xn;
		ny = yn;
		nz = zn;
	}

	Pt screen() {
		int sx = 200 + x*200;
		int sy = 200 - y*200;
		return Pt(sx, sy, z, r, g, b, nx, ny, nz, index);
	}

	int getX(int mult) { return 200 + x*mult; }
	int getY(int mult) { return 200 - y*mult; }
	int getZ(int mult) { return z*mult;	}

	void add(int count)
	{
		faces_included.push_back(count);
	}

	void normalize()
	{
		float magnitude = sqrt(nx*nx + ny*ny + nz*nz);
		nx = nx / magnitude;
		ny = ny / magnitude;
		nz = nz / magnitude;

		double Lx = 1 / sqrt(3);
		double Ly = 1 / sqrt(3);
		double Lz = -1 / sqrt(3);

		double LN = (nx*Lx + ny*Ly + nz*Lz);

		double Rx = 2 * LN*nx - Lx;
		double Ry = 2 * LN*ny - Ly;
		double Rz = 2 * LN*nz - Lz;
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
	}
};

class normalized_vec
{
public: 
	float x, y, z, angle;
	
	normalized_vec() {}

	normalized_vec(Pt point1, Pt point2, Pt point3)
	{
		Pt v(point2.x - point1.x, point2.y - point1.y, point2.z - point1.z, -1);
		Pt w(point3.x - point1.x, point3.y - point1.y, point3.z - point1.z, -1);

		x = -v.y*w.z + v.z*w.y;
		y = -v.z*w.x + v.x*w.z;
		z = -v.x*w.y + v.y*w.x;
		float magnitude = sqrt(x*x + y*y + z*z);
		x = x / magnitude;
		y = y / magnitude;
		z = z / magnitude;
	}
};

class Face
{
public:
	Pt point1;
	Pt point2;
	Pt point3;

	normalized_vec normal;

	Face(Pt one, Pt two, Pt three)
	{
		point1 = one;
		point2 = two;
		point3 = three;

		normal = normalized_vec(point1, point2, point3);
	}

};

void drawPolygon(vector<Pt> points, color c);
void drawGouraudPolygon(vector<Pt> points);
void drawPhongPolygon( vector<Pt> points);

#endif