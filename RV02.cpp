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
	viewer viewTransformed("Filtered");

	/*---------------------*/
	/* images & channels   */
    /*---------------------*/
    image img;                              // normalized (color) image
	channel8  src;  // source picture       // 8-bit-image (source)
    channel8  dst;  // destination picture  // 8-bit-image (source) 


	/*-----------------------------*/
	/* Image processing components */
	/*-----------------------------*/

    // object to split image into hue, saturation and intensity
	// hue        = Farbton
	// saturation = Farbsдttigung
	// intensity  = Intensitдt (Grauwert)
    splitImageToHSI splitter;



	/******************************************/
	/*    the program                         */
	/******************************************/

    // load the source image
    loader.load("shaft.bmp",img);
    
    // extract the intensity channel only
    splitter.getIntensity(img,src);

    // determine image size
    const int rowSize    = src.rows();
    const int columnSize = src.columns();


    // set destination size to source size 
    dst.resize(rowSize,columnSize,0,false,true);


    
    Median(src,dst,9,9);
	//Sobel(src, dst, *(new channel8));


	// view pictures
    view.show(src);
    viewTransformed.show(dst);

	

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

		ubyte medianValue		= (maskSize + 1) / 2;
		ubyte histogram[256]	= {0};
		ubyte accHistogram[256] = {0};
		
		int x, y, mx, my;

		for(y = maskHalfHeight + 1; y < PicSizeY - maskHalfHeight - 1; y++){
			for(x = maskHalfWidth + 1; x < PicSizeX - maskHalfWidth - 1; x++){

				std::fill(histogram, histogram + 256, 0);
				
				for(my = y - maskHalfHeight; my <= y + maskHalfHeight; my++){
					for(mx = x - maskHalfWidth; mx <= x + maskHalfWidth; mx++){
						histogram[sPic[my][mx]]++;
					}
				}

				accHistogram[0] = histogram[0];
				for(int i = 1; i < 256; i++){
					accHistogram[i] = histogram[i] + accHistogram[i-1];
				}

				ubyte value		= 0;
				int medianIndex = 0;
				for(int i = 0; i < 256; i++){
					medianIndex = i;
					if(value > medianValue) break;
					value = accHistogram[i];
				}

				dPic[y][x] = medianIndex;
			}
		}
	}

	void RV02::Sobel (const channel8& sPic, channel8& GradientPic, channel8& DirectionPic){
		const int PicSizeY = sPic.rows();
		const int PicSizeX = sPic.columns();

		const int MaskSizeX			= 3;
		const int MaskSizeY			= 3;
		const int maskWidth			= (MaskSizeX % 2 == 0 ? MaskSizeX + 1 : MaskSizeX);
		const int maskHeight		= (MaskSizeX % 2 == 0 ? MaskSizeY + 1 : MaskSizeY);
		const int maskSize			= maskWidth * maskHeight;
		const int maskHalfWidth		= maskWidth / 2;
		const int maskHalfHeight	= maskHeight / 2;

		ubyte Gx[3][3] = {
			-1,0,1,
			-2,0,2,
			-1,0,1
		};
		ubyte Gy[3][3] = {
			-1,-2,-1,
			0,0,0,
			1,2,1
		};

		int x, y, mx, my;

		for(y = maskHalfHeight + 1; y < PicSizeY - maskHalfHeight - 1; y++){
			for(x = maskHalfWidth + 1; x < PicSizeX - maskHalfWidth - 1; x++){
				GradientPic[y][x]	= sPic[y][x];
				DirectionPic[y][x]	= sPic[y][x];
			}
		}
	}
};
