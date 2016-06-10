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
//
// $Id: SoBox.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
/*----------------------------HEPVis----------------------------------------*/
/*                                                                          */
/* Node:             SoBox                                                  */
/* Description:      Represents the G4Box Geant Geometry entity             */
/* Author:           Joe Boudreau Nov 11 1996                               */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#ifdef G4VIS_BUILD_OI_DRIVER

// this :
#include "HEPVis/nodes/SoBox.h"

#include <assert.h>
#include <cmath>

#include <Inventor/SbBox.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoScale.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/elements/SoTextureCoordinateElement.h>

// This statement is required
SO_NODE_SOURCE(SoBox)

// Constructor
SoBox::SoBox() {
  // This statement is required
  SO_NODE_CONSTRUCTOR(SoBox);

  // Data fields are initialized like this:
  SO_NODE_ADD_FIELD(fDx,                (1.0));
  SO_NODE_ADD_FIELD(fDy,                (1.0));
  SO_NODE_ADD_FIELD(fDz,                (1.0));
  SO_NODE_ADD_FIELD(alternateRep,       (NULL));
  children = new SoChildList(this);
}

// Destructor
SoBox::~SoBox() {
 delete children;
}


// initClass
void SoBox::initClass(){
  // This statement is required.
  SO_NODE_INIT_CLASS(SoBox,SoShape,"Shape");
}


// generatePrimitives
void SoBox::generatePrimitives(SoAction *action) {
  // This variable is used to store each vertex
  SoPrimitiveVertex pv;

  // Access the stat from the action
  SoState *state = action->getState();

  // See if we have to use a texture coordinate function,
  // rather than generating explicit texture coordinates.
  SbBool useTexFunction=
    (SoTextureCoordinateElement::getType(state) == 
     SoTextureCoordinateElement::FUNCTION);

  // If we need to generate texture coordinates with a function,
  // we'll need an SoGLTextureCoordinateElement.  Otherwise, we'll
  // set up the coordinates directly.
  const SoTextureCoordinateElement *tce = NULL;
  SbVec4f texCoord;
  if (useTexFunction) {
    tce = SoTextureCoordinateElement::getInstance(state);
  }
  else {
    texCoord[2] = 0.0;
    texCoord[3] = 1.0;
  }
  SbVec3f point, normal;


  //////////////////////////////////////////
  //----------------------------------------
#define GEN_VERTEX(pv,x,y,z,s,t,nx,ny,nz)  \
  point.setValue(x,y,z);                   \
  normal.setValue(nx,ny,nz);               \
  if (useTexFunction) {                    \
    texCoord=tce->get(point,normal);       \
  }                                        \
  else {                                   \
    texCoord[0]=s;                         \
    texCoord[1]=t;                         \
  }                                        \
  pv.setPoint(point);                      \
  pv.setNormal(normal);                    \
  pv.setTextureCoords(texCoord);           \
  shapeVertex(&pv);
  //----------------------------------------
  //////////////////////////////////////////

  const int NPOINTS=8, NFACES=6, NINDICES = NFACES*5;
  int indices[NINDICES] = {3,2,1,0, SO_END_FACE_INDEX,  //z back.
			   4,5,6,7, SO_END_FACE_INDEX,  //z front.
			   0,1,5,4, SO_END_FACE_INDEX,  //y up.
			   1,2,6,5, SO_END_FACE_INDEX,  //x left.
			   2,3,7,6, SO_END_FACE_INDEX,  //y down.
			   3,0,4,7, SO_END_FACE_INDEX}; //x right.

  
  // points for the eight vertices
  float points[NPOINTS][3];
  points[0][0] =  fDx.getValue(); 
  points[0][1] =  fDy.getValue(); 
  points[0][2] = -fDz.getValue();

  points[1][0] = -fDx.getValue();
  points[1][1] =  fDy.getValue();
  points[1][2] = -fDz.getValue();

  points[2][0] = -fDx.getValue(); 
  points[2][1] = -fDy.getValue(); 
  points[2][2] = -fDz.getValue();

  points[3][0] =  fDx.getValue(); 
  points[3][1] = -fDy.getValue(); 
  points[3][2] = -fDz.getValue();

  points[4][0] =  fDx.getValue(); 
  points[4][1] =  fDy.getValue(); 
  points[4][2] =  fDz.getValue();

  points[5][0] = -fDx.getValue(); 
  points[5][1] =  fDy.getValue(); 
  points[5][2] =  fDz.getValue();

  points[6][0] = -fDx.getValue(); 
  points[6][1] = -fDy.getValue(); 
  points[6][2] =  fDz.getValue();

  points[7][0] =  fDx.getValue(); 
  points[7][1] = -fDy.getValue(); 
  points[7][2] =  fDz.getValue();

  float normals[NFACES][3];
  //z back.
  normals[0][0] =  0  ; normals[0][1] =    0; normals [0][2] =  -1;    
  //z front.
  normals[1][0] =  0  ; normals[1][1] =    0; normals [1][2] =   1;    
  //y up.
  normals[2][0] =  0  ; normals[2][1] =    1; normals [2][2] =   0;    
  //x left.
  normals[3][0] = -1  ; normals[3][1] =    0; normals [3][2] =   0;    
  //y down.
  normals[4][0] =  0  ; normals[4][1] =   -1; normals [4][2] =   0;    
  //x right.
  normals[5][0] =  1  ; normals[5][1] =    0; normals [5][2] =   0;    

  float x,y,z;
  int   index;
  for (int nf=0;nf<NFACES;nf++) {
    beginShape(action,TRIANGLE_FAN);
    index = indices[nf * 5];   
    x = points[index][0];
    y = points[index][1];
    z = points[index][2];
    GEN_VERTEX(pv,x,y,z,0.0,0.0,normals[nf][0],normals[nf][1],normals[nf][2]);   
    index = indices[nf * 5 + 1];   
    x = points[index][0];
    y = points[index][1];
    z = points[index][2];
    GEN_VERTEX(pv,x,y,z,0.0,0.0,normals[nf][0],normals[nf][1],normals[nf][2]);   
    index = indices[nf * 5 + 2];   
    x = points[index][0];
    y = points[index][1];
    z = points[index][2];
    GEN_VERTEX(pv,x,y,z,0.0,0.0,normals[nf][0],normals[nf][1],normals[nf][2]);   
    index = indices[nf * 5 + 3];   
    x = points[index][0];
    y = points[index][1];
    z = points[index][2];
    GEN_VERTEX(pv,x,y,z,0.0,0.0,normals[nf][0],normals[nf][1],normals[nf][2]);   
    endShape();
  }
}

// getChildren
SoChildList *SoBox::getChildren() const {
  return children;
}


// computeBBox
void SoBox::computeBBox(SoAction *, SbBox3f &box, SbVec3f &center ){
  SbVec3f vmin(-fDx.getValue(),-fDy.getValue(),-fDz.getValue()), 
          vmax( fDx.getValue(), fDy.getValue(), fDz.getValue());
  center.setValue(0,0,0);
  box.setBounds(vmin,vmax);
}




// updateChildren
void SoBox::updateChildren() {


  // Redraw the G4Box....

  assert(children->getLength()==1);
  SoSeparator       *sep                = (SoSeparator *)  ( *children)[0];
  SoScale           *scale              = (SoScale *)( sep->getChild(0));
  //SoCube            *cube               = (SoCube  *)( sep->getChild(1));
  scale->scaleFactor.setValue(fDx.getValue(), fDy.getValue(), fDz.getValue());
}

// generateChildren
void SoBox::generateChildren() {

  // A box consists of a set of scale factors and a 
  // cube.

  assert(children->getLength() ==0);
  SoSeparator      *sep              = new SoSeparator(); 
  SoScale          *scale            = new SoScale();
  SoCube           *cube             = new SoCube();

  sep->addChild(scale);
  sep->addChild(cube);
  children->append(sep);
}

// generateAlternateRep
void SoBox::generateAlternateRep() {

  // This routine sets the alternate representation to the child
  // list of this mode.  

  if (children->getLength() == 0) generateChildren();
  updateChildren();
  alternateRep.setValue((SoSeparator *)  ( *children)[0]);
}

// clearAlternateRep
void SoBox::clearAlternateRep() {
  alternateRep.setValue(NULL);
}

#endif
