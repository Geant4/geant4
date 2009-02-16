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
#ifdef G4VIS_BUILD_OI_DRIVER

/*----------------------------HEPVis----------------------------------------*/
/*                                                                          */
/* Node:             SoGL2PSAction                                          */
/* Author:           Guy Barrand                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

// this :
#include <HEPVis/actions/SoGL2PSAction.h>

// Inventor :
#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/errors/SoDebugError.h>

#include "Geant4_gl2ps.h"

#include <stdio.h>

SO_ACTION_SOURCE(SoGL2PSAction)
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::initClass(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_ACTION_INIT_CLASS(SoGL2PSAction,SoGLRenderAction);
}
//////////////////////////////////////////////////////////////////////////////
SoGL2PSAction::SoGL2PSAction(
 const SbViewportRegion& aViewPortRegion
)
:SoGLRenderAction(aViewPortRegion)
,G4OpenGL2PSAction()
,fFileName("out.ps")
,fFile(0)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_ACTION_CONSTRUCTOR(SoGL2PSAction);
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::setViewport(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  // Useful ?? L.Garnier 02/2009
#ifdef __COIN__
#else //SGI
  const SbViewportRegion& vpr = getViewportRegion();
  SoViewportRegionElement::set(getState(),vpr);
 
  const SbVec2s& win = vpr.getWindowSize();
  fViewport[0] = 0;
  fViewport[1] = 0;
  fViewport[2] = win[0];
  fViewport[3] = win[1];
#endif
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::beginTraversal(
 SoNode* aNode
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile) {
#ifdef __COIN__
    const SbViewportRegion& vpr = getViewportRegion();
    SoViewportRegionElement::set(getState(),vpr);
    G4gl2psBegin();
    traverse(aNode);
    gl2psEndPage();        
#else //SGI
    SoGLRenderAction::beginTraversal(aNode);
#endif
  } else {
    SoGLRenderAction::beginTraversal(aNode);
  }
}

#endif
