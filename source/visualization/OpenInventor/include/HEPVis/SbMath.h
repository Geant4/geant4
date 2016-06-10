//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: SbMath.h 66373 2012-12-18 09:41:34Z gcosmo $
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
