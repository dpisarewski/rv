/*
   RV01: Affine Transformation
   
   Autor: .....................
   HAW-University of Applied Sciences - Hamburg,Germany

 */ 

#include "ltiObject.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <list>
#include <stdio.h>

#include "RV01.h"
#include "ltiTimer.h"
#include "ltiBMPFunctor.h"
#include "ltiViewer.h"
#include "ltiSplitImg.h"
#include "gtk.h"
#include "ltiGtkServer.h"
#include "RV01.h"

using std::cout;
using std::endl;


namespace lti {

  void RV01::operator()(int argc,char *argv[]) {

	/*********************************************/
	/**** has always to be started (AMei)     ****/
    /**** if program is without gtk-Widgets   ****/
	/*********************************************/

	gtkServer server;
    server.start();

	/******************************************/
	/**** instantiation of used components ****/
	/******************************************/

	/*---------------------*/
	/* loaders and viewers */
    /*---------------------*/
    loadBMP loader;                         // object for loading .bmp-images

    viewer view("Original");                // object for visualizing images
	viewer viewNeighbor("affinNeighbor");
	viewer viewBilinear("affinBilinear");
	viewer viewTransformNeighbor("fourPointTransformNeighbor");
	viewer viewTransformBilinear("fourPointTransformBilinear");

	/*---------------------*/
	/* images & channels   */
    /*---------------------*/
    image img;                              // normalized (color) image
	channel8 src;  // source picture       // 8-bit-image (source)
    channel8 neighbor;  // destination picture  // 8-bit-image (source) 
	channel8 bilinear;
	channel8 transformNeighbor;
	channel8 transformBilinear;


	/*-----------------------------*/
	/* Image processing components */
	/*-----------------------------*/

    // object to split image into hue, saturation and intensity
	// hue        = Farbton
	// saturation = Farbsättigung
	// intensity  = Intensität (Grauwert)
    splitImageToHSI splitter;



	/******************************************/
	/*    the program                         */
	/******************************************/

    // load the source image
    loader.load("Kalib.bmp",img);
    
    // extract the intensity channel only
    splitter.getIntensity(img,src);

    // determine image size
    const int rowSize    = src.rows();
    const int columnSize = src.columns();


    // set destination size to source size 
    neighbor.resize(rowSize,columnSize,0,false,true);
	bilinear.resize(rowSize,columnSize,0,false,true);
	transformNeighbor.resize(rowSize,columnSize,0,false,true);
	transformBilinear.resize(rowSize,columnSize,0,false,true);


    // copy source pixels to destination pixels
	for(int y=0; y<rowSize; y++) {
		for(int x=0; x<columnSize; x++) {
			neighbor[y][x] = src[y][x];
			bilinear[y][x] = src[y][x];
			transformNeighbor[y][x] = src[y][x];
			transformBilinear[y][x] = src[y][x];
		}
	}

	// transform pictures
	affinTransform(src, neighbor, bilinear);
	fourPointTransform(src, transformNeighbor, transformBilinear);

	// view pictures
    view.show(src);
    viewNeighbor.show(neighbor);
	viewBilinear.show(bilinear);
	viewTransformNeighbor.show(transformNeighbor);
	viewTransformBilinear.show(transformBilinear);

    getchar();

    }

	void RV01::affinTransform(const channel8& src, channel8& neighbor, channel8& bilinear) {
		const int PicSizeY = src.rows();
		const int PicSizeX = src.columns();
		// Transformationsparameter
		double a0 = 160.792;
		double a1 = 0.441667;
		double a2 = -0.13;
		double b0 = 218.742;
		double b1 = 0.00166667;
		double b2 = 0.4925;

		// affine transformations ergebnis
		double xa;
		double ya;

		// alle Pixel durchgehen
		for(int y=0; y<PicSizeY; y++) {
			for(int x=0; x<PicSizeX; x++) {
			// formeln für affine transformation
				xa = a0 + a1*x + a2*y;
				ya = b0 + b1*x + b2*y;
				neighbor[y][x] = nearestNeighbor(xa, ya, src);
				bilinear[y][x] = bilinInterpol(xa, ya, src);
			}
		}
    }

    int RV01::nearestNeighbor(const double dX, const double dY, const channel8& src) {
		return src[(int)dY + 0.5][(int)dX + 0.5];
    }

	int RV01::bilinInterpol(const double dX, const double dY, const channel8& src) {

		int fx1 = src[dY][dX] + (dX - (int)dX) * (src[dY][dX+1] - src[dY][dX]);
		int fx2 = src[dY+1][dX] + (dX - (int)dX) * (src[dY+1][dX+1] - src[dY+1][dX]);
		int fy = fx1 + (dY - (int)dY) * (fx2 - fx1);

		return correctBorder(fy, 256);
	}

	bool RV01::checkBorder(int coordinate, int maxsize) {
		return coordinate >= 0 && coordinate < maxsize-1;
	}

	int RV01::correctBorder(int coordinate, int maxsize) {
		if (coordinate < 0) return 0;
		else if (coordinate > maxsize-1) return maxsize-1;
		else return coordinate;
	}

	void RV01::fourPointTransform(const channel8& src,  channel8& transformNeighbor, channel8& transformBilinear) {
		
		const int PicSizeY = src.rows();
		const int PicSizeX = src.columns();

		// check all pixels
		for(int y=0; y<PicSizeY; y++) {
			for(int x=0; x<PicSizeX; x++) {

				// normalize target coordinates
				double xz = ((double) x) / (PicSizeX-1);
				double yz = ((double) y) / (PicSizeY-1);

				double p1 = (1 - xz) * (1 - yz);
				double p2 = xz * (1 - yz);
				double p3 = xz * yz;
				double p4 = (1 - xz) * yz;
				int xq = p1*196 + p2*582 + p3*666 + p4*99;
				int yq = p1*101 + p2*95 + p3*473 + p4*468;
				if (checkBorder(xq, PicSizeX) && checkBorder(yq, PicSizeY)) {
					//transform[y][x] = src[yq][xq];
					transformNeighbor[y][x] = nearestNeighbor(xq, yq, src);
					transformBilinear[y][x] = bilinInterpol(xq, yq, src);
				}
			}
		}
	}
};
