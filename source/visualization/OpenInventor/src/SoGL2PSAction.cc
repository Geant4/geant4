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
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  setFileName("out.ps");
  SO_ACTION_CONSTRUCTOR(SoGL2PSAction);
}

bool SoGL2PSAction::enableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fFile = ::fopen(fFileName,"w");
  if(!fFile) {
    SoDebugError::post("SoGL2PSAction::enableFileWriting",
                       "Cannot open file %s",fFileName);
    return false;
  }
#ifdef __COIN__
#else //SGI
  const SbViewportRegion& vpr = getViewportRegion();
  SoViewportRegionElement::set(getState(),vpr);
  G4gl2psBegin();
#endif
  return true;
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::disableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef __COIN__
#else //SGI
  gl2psEndPage();        
#endif
  ::fclose(fFile);
  fFile = 0;
}

//////////////////////////////////////////////////////////////////////////////
bool SoGL2PSAction::addBitmap(
 int aWidth
,int aHeight
,float aXorig
,float aYorig
,float aXmove
,float aYmove
)
/////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fFile) return false;
  GLboolean valid;
  glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID,&valid);
  if(!valid) return false;
  float pos[4];
  glGetFloatv(GL_CURRENT_RASTER_POSITION,pos);
  int xoff = -(int)(aXmove + aXorig);
  int yoff = -(int)(aYmove + aYorig);
  int x = (int)(pos[0] + xoff);
  int y = (int)(pos[1] + yoff);
  // Should clip against viewport area :
  GLint vp[4];
  glGetIntegerv(GL_VIEWPORT,vp);
  GLsizei w = aWidth;
  GLsizei h = aHeight;
  if(x+w>(vp[0]+vp[2])) w = vp[0]+vp[2]-x;
  if(y+h>(vp[1]+vp[3])) h = vp[1]+vp[3]-y;
  int s = 3 * w * h;
  if(s<=0) return false;
  float* image = (float*)::malloc(s * sizeof(float));
  if(!image) return false;
  glReadPixels(x,y,w,h,GL_RGB,GL_FLOAT,image);
  GLint status = gl2psDrawPixels(w,h,xoff,yoff,GL_RGB,GL_FLOAT,image);
  ::free(image);
  return (status!=GL2PS_SUCCESS ? false : true);
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
    // Should have already do G4gl2psBegin() before
    SoGLRenderAction::beginTraversal(aNode);
    // Should do gl2psEndPage() after
#endif
  } else {
    SoGLRenderAction::beginTraversal(aNode);
  }
}

#endif
