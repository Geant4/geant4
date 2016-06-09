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
#ifdef G4VIS_BUILD_OI_DRIVER

/*----------------------------HEPVis----------------------------------------*/
/*                                                                          */
/* Node:             SbPainter                                              */
/* Author:           Guy Barrand                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

// this :
#include <HEPVis/SbPainter.h>

//////////////////////////////////////////////////////////////////////////////
SbPainter::SbPainter(
)
:fWindowWidth(0)
,fWindowHeight(0)
,fRed(0.)
,fGreen(0.)
,fBlue(0.)
,fEdges(FALSE)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}
//////////////////////////////////////////////////////////////////////////////
SbPainter::~SbPainter(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}
//////////////////////////////////////////////////////////////////////////////
void SbPainter::setWindowSize(
 int aWindowWidth
,int aWindowHeight
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fWindowWidth = aWindowWidth;
  fWindowHeight = aWindowHeight;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainter::enableEdges(
 SbBool aValue
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fEdges = aValue;
}

#endif
