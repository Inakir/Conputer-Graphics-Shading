#include "PolygonDrawer.h"
#include "ScanConvert.h"
#include <algorithm>
#include <math.h>
#include "ScanConvert.h"

class Edge
{
public:
	float slopeRecip;
	float maxY;
	float currentX;

	float currentFR, currentFG, currentFB;
	float fIncrR, fIncrG, fIncrB;

	float currentNX, currentNY, currentNZ;
	float fIncrNX, fIncrNY, fIncrNZ;

	float currentZ;
	float IncrZ;

	bool operator < ( const Edge &e )
	{
		if ( currentX == e.currentX )
		{
			return slopeRecip < e.slopeRecip;
		}
		else
		{
			return currentX < e.currentX;
		}
	}
};

vector<vector<Edge> > activeEdgeTable;
vector<Edge> activeEdgeList;

void buildActiveEdgeTable ( vector<Pt> &points )
{
	int i;

	activeEdgeTable.clear ( );

	// add rows equal to height of image to active edge table
	for ( i = 0; i < ImageH; i++ )
	{
		vector<Edge> row;

		activeEdgeTable.push_back ( row );
	}

	for ( i = 0; i < points.size ( ); i++ )
	{
		Edge e;
		int next = ( i + 1 ) % points.size ( );

		// ignore horizontal lines
		if ( points [ i ].y == points [ next ].y )
		{
			continue;
		}
		e.maxY = max ( points [ i ].y, points [ next ].y );
		e.slopeRecip = (points[i].x  - points[next].x ) / (float)(points[i].y - points[next].y);
		e.fIncrR     = (points[i].r  - points[next].r ) / (float)(points[i].y - points[next].y);
		e.fIncrG     = (points[i].g  - points[next].g ) / (float)(points[i].y - points[next].y);
		e.fIncrB     = (points[i].b  - points[next].b ) / (float)(points[i].y - points[next].y);
		e.fIncrNX    = (points[i].nx - points[next].nx) / (float)(points[i].y - points[next].y);
		e.fIncrNY    = (points[i].ny - points[next].ny) / (float)(points[i].y - points[next].y);
		e.fIncrNZ    = (points[i].nz - points[next].nz) / (float)(points[i].y - points[next].y);
		e.fIncrNZ	 = (points[i].nz - points[next].nz) / (float)(points[i].y - points[next].y);
		e.IncrZ      = (points[i].z  - points[next].z ) / (float)(points[i].y - points[next].y);

		if ( points [ i ].y == e.maxY )
		{
			e.currentX = points [next].x;

			e.currentFR = points[next].r;
			e.currentFG = points[next].g;
			e.currentFB = points[next].b;

			e.currentNX = points[next].nx;
			e.currentNY = points[next].ny;
			e.currentNZ = points[next].nz;
			
			e.currentZ = points[next].z;

			activeEdgeTable [ points [ next ].y ].push_back ( e );
		}
		else
		{
			e.currentX = points [i].x;

			e.currentFR = points[i].r;
			e.currentFG = points[i].g;
			e.currentFB = points[i].b;

			e.currentNX = points[i].nx;
			e.currentNY = points[i].ny;
			e.currentNZ = points[i].nz;
			
			e.currentZ = points[i].z;

			activeEdgeTable [ points [ i ].y ].push_back ( e );
		}
	}
}

// assumes all vertices are within window!!!
void drawPolygon ( vector<Pt> points, color c )
{
	int x, y, z, i;
	float dZ = 0;
	activeEdgeList.clear ( );
	buildActiveEdgeTable ( points );

	for ( y = 0; y < ImageH; y++ )
	{
		// add edges into active Edge List
		for ( i = 0; i < activeEdgeTable [ y ].size ( ); i++ )
		{
			activeEdgeList.push_back ( activeEdgeTable [ y ] [ i ] );
		}

		// delete edges from active Edge List
		for ( i = 0; i < activeEdgeList.size ( ); i++ )
		{
			if ( activeEdgeList [ i ].maxY <= y )
			{
				activeEdgeList.erase ( activeEdgeList.begin ( ) + i );
				i--;
			}
		}

		// sort according to x value... a little expensive since not always necessary
		sort ( activeEdgeList.begin ( ), activeEdgeList.end ( ) );

		// draw scan line
		for ( i = 0; i < activeEdgeList.size ( ); i += 2 )
		{
			dZ = (activeEdgeList[i + 1].currentZ - activeEdgeList[i].currentZ) / ((int)ceil(activeEdgeList[i + 1].currentX) - (int)ceil(activeEdgeList[i].currentX));
			z = activeEdgeList[i].currentZ;

			for ( x = (int)ceil ( activeEdgeList [ i ].currentX ); x < activeEdgeList [ i + 1 ].currentX; x++ )
			{
				setFramebuffer ( x, y, z, c.r, c.g, c.b );
				z += dZ;
			}
		}

		// update edges in active edge list
		for ( i = 0; i < activeEdgeList.size ( ); i++ )
		{
			activeEdgeList [ i ].currentX += activeEdgeList [ i ].slopeRecip;
			activeEdgeList[i].currentZ += activeEdgeList[i].IncrZ;

		}
	}
}

// assumes all vertices are within window!!!
void drawGouraudPolygon(vector<Pt> points)
{
	int x, y, z, i;
	float r, g, b;
	float dFr = 0;
	float dFg = 0;
	float dFb = 0;
	float dZ = 0;

	activeEdgeList.clear();
	buildActiveEdgeTable(points);

	for (y = 0; y < ImageH; y++)
	{
		// add edges into active Edge List
		for (i = 0; i < activeEdgeTable[y].size(); i++)
		{
			activeEdgeList.push_back(activeEdgeTable[y][i]);
		}

		// delete edges from active Edge List
		for (i = 0; i < activeEdgeList.size(); i++)
		{
			if (activeEdgeList[i].maxY <= y)
			{
				activeEdgeList.erase(activeEdgeList.begin() + i);
				i--;
			}
		}

		// sort according to x value... a little expensive since not always necessary
		sort(activeEdgeList.begin(), activeEdgeList.end());
		
		// draw scan line
		for (i = 0; i < activeEdgeList.size(); i += 2)
		{
			dFr = (activeEdgeList[i + 1].currentFR - activeEdgeList[i].currentFR) / ((int)ceil(activeEdgeList[i + 1].currentX) - (int)ceil(activeEdgeList[i].currentX));
			dFg = (activeEdgeList[i + 1].currentFG - activeEdgeList[i].currentFG) / ((int)ceil(activeEdgeList[i + 1].currentX) - (int)ceil(activeEdgeList[i].currentX));
			dFb = (activeEdgeList[i + 1].currentFB - activeEdgeList[i].currentFB) / ((int)ceil(activeEdgeList[i + 1].currentX) - (int)ceil(activeEdgeList[i].currentX));
			dZ = (activeEdgeList[i + 1].currentZ - activeEdgeList[i].currentZ) / ((int)ceil(activeEdgeList[i + 1].currentX) - (int)ceil(activeEdgeList[i].currentX));

			r = activeEdgeList[i].currentFR;
			g = activeEdgeList[i].currentFG;
			b = activeEdgeList[i].currentFB;
			z = activeEdgeList[i].currentZ;

			for (x = (int)ceil(activeEdgeList[i].currentX); x < activeEdgeList[i + 1].currentX; x++)
			{
				setFramebuffer(x, y, z, r, g, b);
				r += dFr;
				g += dFg;
				b += dFb;
				z += dZ;
			}
		}

		// update edges in active edge list
		for (i = 0; i < activeEdgeList.size(); i++)
		{
			activeEdgeList[i].currentX += activeEdgeList[i].slopeRecip;
			activeEdgeList[i].currentFR += activeEdgeList[i].fIncrR;
			activeEdgeList[i].currentFG += activeEdgeList[i].fIncrG;
			activeEdgeList[i].currentFB += activeEdgeList[i].fIncrB;
			activeEdgeList[i].currentZ += activeEdgeList[i].IncrZ;

		}
	}
}

void setColor(float &r, float &g, float &b, float nx, float ny, float nz)
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
	if (LN >= 0)
	{
		r+=(1) * (.7*LN + .5*pow(RE, n));
		g+=(1) * (0 * LN + .5*pow(RE, n));
		b+=(1) * (0 * LN + .5*pow(RE, n));
	}
}

// assumes all vertices are within window!!!
void drawPhongPolygon(vector<Pt> points)
{
	int x, y, z, i;
	
	float r, g, b;
	float nx, ny, nz = 0;
	float dNX, dNY, dNZ , dZ= 0;

	activeEdgeList.clear();
	buildActiveEdgeTable(points);

	for (y = 0; y < ImageH; y++)
	{
		// add edges into active Edge List
		for (i = 0; i < activeEdgeTable[y].size(); i++)
		{
			activeEdgeList.push_back(activeEdgeTable[y][i]);
		}

		// delete edges from active Edge List
		for (i = 0; i < activeEdgeList.size(); i++)
		{
			if (activeEdgeList[i].maxY <= y)
			{
				activeEdgeList.erase(activeEdgeList.begin() + i);
				i--;
			}
		}

		// sort according to x value... a little expensive since not always necessary
		sort(activeEdgeList.begin(), activeEdgeList.end());

		// draw scan line
		for (i = 0; i < activeEdgeList.size(); i += 2)
		{
			dNX = (activeEdgeList[i + 1].currentNX - activeEdgeList[i].currentNX) / ((int)ceil(activeEdgeList[i + 1].currentX) - (int)ceil(activeEdgeList[i].currentX));
			dNY = (activeEdgeList[i + 1].currentNY - activeEdgeList[i].currentNY) / ((int)ceil(activeEdgeList[i + 1].currentX) - (int)ceil(activeEdgeList[i].currentX));
			dNZ = (activeEdgeList[i + 1].currentNZ - activeEdgeList[i].currentNZ) / ((int)ceil(activeEdgeList[i + 1].currentX) - (int)ceil(activeEdgeList[i].currentX));
			dZ = (activeEdgeList[i + 1].currentZ - activeEdgeList[i].currentZ) / ((int)ceil(activeEdgeList[i + 1].currentX) - (int)ceil(activeEdgeList[i].currentX));

			nx = activeEdgeList[i].currentNX;
			ny = activeEdgeList[i].currentNY;
			nz = activeEdgeList[i].currentNZ;
			z = activeEdgeList[i].currentZ;

			for (x = (int)ceil(activeEdgeList[i].currentX); x < activeEdgeList[i + 1].currentX; x++)
			{
				setColor(r, g, b, nx, ny, nz);
				setFramebuffer(x, y, z, r, g, b);
				nx += dNX;
				ny += dNY;
				nz += dNZ;
				z += dZ;
			}
		}

		// update edges in active edge list
		for (i = 0; i < activeEdgeList.size(); i++)
		{
			activeEdgeList[i].currentX += activeEdgeList[i].slopeRecip;
			activeEdgeList[i].currentNX += activeEdgeList[i].fIncrNX;
			activeEdgeList[i].currentNY += activeEdgeList[i].fIncrNY;
			activeEdgeList[i].currentNZ += activeEdgeList[i].fIncrNZ;
			activeEdgeList[i].currentZ += activeEdgeList[i].IncrZ;
		}
	}
}