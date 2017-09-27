#ifndef SCAN_CONVERT_H
#define SCAN_CONVERT_H

#define ImageW 400
#define ImageH 400

struct color {
	float r, g, b;		// Color (R,G,B values)
};

// Clears framebuffer to black
void clearFramebuffer();

// Sets pixel x,y to the color RGB
void setFramebuffer(int x, int y, int z, float R, float G, float B);

#endif