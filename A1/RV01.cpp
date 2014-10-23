/*
   RV02: Median und Sobel
   
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

#include "RV02.h"
#include "ltiTimer.h"
#include "ltiBMPFunctor.h"
#include "ltiViewer.h"
#include "ltiSplitImg.h"
#include "gtk.h"
#include "ltiGtkServer.h"

using std::cout;
using std::endl;


namespace lti {

  void RV02::operator()(int argc,char *argv[]) {

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


	/*---------------------*/
	/* images & channels   */
    /*---------------------*/
    image img;                              // normalized (color) image
	channel8  src;  // source picture       // 8-bit-image (source)
    channel8  med;  // destination picture  // 8-bit-image (source) 
    channel8  grad;  // destination picture  // 8-bit-image (source) 
    channel8  dir;  // destination picture  // 8-bit-image (source) 


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
    loader.load("shaft.bmp",img);
    
    // extract the intensity channel only
    splitter.getIntensity(img,src);
    splitter.getIntensity(img,med);
    splitter.getIntensity(img,grad);
    splitter.getIntensity(img,dir);

    // determine image size
    const int rowSize    = src.rows();
    const int columnSize = src.columns();


    // set destination size to source size 
    med.resize(rowSize,columnSize,0,false,true);
    grad.resize(rowSize,columnSize,0,false,true);
    dir.resize(rowSize,columnSize,0,false,true);

	// MEDIAN
    Median(src,med,9,9);

    view.show(src);
	viewer viewMed("Median");
    viewMed.show(med);

	// SOBEL
	Sobel(src, grad, dir);

    view.show(src);
	viewer viewGrad("Gradient");
	viewer viewDir("Direction");
    viewGrad.show(grad);
    viewDir.show(dir);
	


	// view pictures

    getchar();

  }



  

  /***************************************************************************/
  /* Function definition: ----- Median-operator----                          */
  /***************************************************************************/
  void RV02::Median(  const	     channel8& sPic, 	// source picture 
								 channel8& dPic, 	// destination picture
	                   const int MaskSizeX,		    // mask size in x-direction
					   const int MaskSizeY		 	// mask size in y-direction
					   )
  {
	const int PicSizeY = sPic.rows();
	const int PicSizeX = sPic.columns();

	const int maskWidth			= (MaskSizeX % 2 == 0 ? MaskSizeX + 1 : MaskSizeX);
	const int maskHeight		= (MaskSizeX % 2 == 0 ? MaskSizeY + 1 : MaskSizeY);
	const int maskSize			= maskWidth * maskHeight;
	const int maskHalfWidth		= maskWidth / 2;
	const int maskHalfHeight	= maskHeight / 2;

	// Gesuchert Wert im akkumulierten Histogramm
	ubyte accuHistSearched			= (maskSize + 1) / 2;
	ubyte histogram[256]	= {0};
	ubyte accHistogram[256] = {0};
	
	int x, y, mx, my;

	for(y = maskHalfHeight; y < PicSizeY - maskHalfHeight; y++){
		for(x = maskHalfWidth; x < PicSizeX - maskHalfWidth; x++){

			// Alte Maskenhinstogrammwerte leeren
			std::fill(histogram, histogram + 256, 0);

			// Initialwerte
			ubyte minGrauwert = 255;
			ubyte maxGrauwert = 0;
			ubyte h = 0;

			// Maskenhistrogramm berechnen
			for(my = y - maskHalfHeight; my <= y + maskHalfHeight; my++){
				for(mx = x - maskHalfWidth; mx <= x + maskHalfWidth; mx++){
					histogram[sPic[my][mx]]++;
					h = sPic[my][mx];
					if (h > maxGrauwert) {
						maxGrauwert = h;
					}
					if (h < minGrauwert){
						minGrauwert = h;
					}
				}
			}



			// Akkumuliertes Histogrammwert berechnen
			// Der erste Index dessen Wert im akkumulierten Histogramm den accuHistValue ueberschreitet
			// ist der gesuchte Medianwert.
			int accuHistvalue = 0;
			for(int i = minGrauwert; i <= maxGrauwert; i++){
				accuHistvalue += histogram[i];
				if(accuHistvalue >= accuHistSearched) {
					dPic[y][x] = i;
					break;
				}
			}
		}
	}
  }
	void RV02::Sobel (const channel8& sPic, channel8& GradientPic, channel8& DirectionPic){
		const int PicSizeY = sPic.rows();
		const int PicSizeX = sPic.columns();

		int x, y;
		int gx, gy;

		for(x = 1; x < PicSizeX - 1; x++){
			for(y = 1; y < PicSizeY - 1; y++){

				// Gx berechnen
				gx =   -sPic[y-1][x-1] + 0 + sPic[y-1][x+1]
					 -2*sPic[y][x-1]   + 0 + 2*sPic[y][x+1]
					   -sPic[y+1][x-1] + 0 + sPic[y+1][x+1];
			    gx /= 4;

				// Gy berechnen
				gy =   - sPic[y-1][x-1] - 2*sPic[y-1][x] - sPic[y-1][x+1]
					   +	 0			+		0		 +		 0
					   + sPic[y+1][x-1] + 2*sPic[y+1][x] + sPic[y+1][x+1];
				gy /= 4;

				// Gradientenbetrag
				GradientPic[y][x]	= sqrt(gx*gx + gy*gy);

				// Winkel: 180 <= angle < 540
				double angle = atan2((double)gy,(double)gx)/PI*180+360.0;
				// Abbildung auf die Winkelbereiche
				int dir = ((int)(angle+22.5)%360)/45;
				DirectionPic[y][x]	= dir;
			}
		}
	}

};
