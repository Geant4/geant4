//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: SbMath.h,v 1.3 2004/12/07 23:40:59 perl Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
#ifndef HEPVis_SbMath_h
#define HEPVis_SbMath_h

#include <cmath>
#ifndef M_PI
#define M_PI       3.1415926535897931160E0
#define M_PI_2     1.5707963267948965580E0  
#endif

#define SbMinimum(a,b) ((a)<(b)?a:b)
#define SbMaximum(a,b) ((a)>(b)?a:b)

#define FCOS(x)   ((float)cos((double)(x)))
#define FSIN(x)   ((float)sin((double)(x)))
#define FACOS(x)  ((float)acos((double)(x)))
#define FASIN(x)  ((float)asin((double)(x)))
#define FTAN(x)   ((float)tan((double)(x)))
#define FATAN(x)  ((float)atan((double)(x)))
#define FSQRT(x)  ((float)sqrt((double)(x)))
#define FPOW(x,y) ((float)pow((double)(x),(double)(y)))
#define FLOG(x)   ((float)log((double)(x)))
#define FLOG10(x) ((float)log10((double)(x)))
#define FFLOOR(x) ((float)floor((double)(x)))
#define FFABS(x)  ((float)fabs((double)(x)))
#define FCEIL(x)  ((float)ceil((double)(x)))

#endif
