/*
 * Copyright (C) 1998, 1999, 2000, 2001
 * Lehrstuhl fuer Technische Informatik, RWTH-Aachen, Germany
 * 
 * This file is part of the LTI-Computer Vision Library (LTI-Lib)
 *
 * The LTI-Lib is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License (LGPL)
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * The LTI-Lib is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public 
 * License along with the LTI-Lib; see the file LICENSE.  If
 * not, write to the Free Software Foundation, Inc., 59 Temple Place -
 * Suite 330, Boston, MA 02111-1307, USA.  
 */ 

 
#ifndef _LTI_TESTER_H
#define _LTI_TESTER_H



#include "ltiImage.h"


#define WITHOUT_GTK_WIDGETS 


namespace lti {

  /* use this object to test your code:  just reimplement the operator() code
     in the ltiTester.cpp file
     the main program will call tester()!! 
  */
  class RV03 {
  public:
    void operator()(int argc,char *argv[]);


  private:
    void CenterFrequencyDomain   (const channel   real,  const channel  imag, 
		                                channel &creal,        channel &cimag);
    void UnCenterFrequencyDomain (const channel  creal,  const channel  cimag, 
		                                channel &real,         channel &imag);
	double distance(int x, int y, int centerX, int centerY);
	bool isNear(int x, int bound, int dev);


  };
}

#endif
