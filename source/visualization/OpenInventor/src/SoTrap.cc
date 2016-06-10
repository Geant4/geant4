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
// $Id: SoTrap.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
/*-----------------------------HEPVis----------------------------------------*/
/*                                                                           */
/* Node:             SoTrap                                                  */
/* Description:      Represents the G4Trap Geant Geometry entity             */
/* Author:           Joe Boudreau Nov 11 1996                                */
/*                                                                           */
/*                                                                           */
/*---------------------------------------------------------------------------*/

#ifdef G4VIS_BUILD_OI_DRIVER

// this :
#include "HEPVis/nodes/SoTrap.h"

#include <assert.h>
#include <cmath>
#include <Inventor/SbBox.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoNormal.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoNormalBinding.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/elements/SoTextureCoordinateElement.h>

#include "HEPVis/SbMath.h"

// This statement is required
SO_NODE_SOURCE(SoTrap)

// Constructor
SoTrap::SoTrap() {
  // This statement is required
  SO_NODE_CONSTRUCTOR(SoTrap);

  // Data fields are initialized like this:
  SO_NODE_ADD_FIELD(pDz,                 (1.0));
  SO_NODE_ADD_FIELD(pTheta,              (0.0));
  SO_NODE_ADD_FIELD(pPhi,                (0.0));
  SO_NODE_ADD_FIELD(pDy1,                (1.0));
  SO_NODE_ADD_FIELD(pDx1,                (1.0));
  SO_NODE_ADD_FIELD(pDx2,                (1.0));
  SO_NODE_ADD_FIELD(pDy2,                (1.0));
  SO_NODE_ADD_FIELD(pDx3,                (1.0));
  SO_NODE_ADD_FIELD(pDx4,                (1.0));
  SO_NODE_ADD_FIELD(pAlp1,               (0.0));
  SO_NODE_ADD_FIELD(pAlp2,               (0.0));
  SO_NODE_ADD_FIELD(alternateRep,        (NULL));
  children = new SoChildList(this);
}

// Destructor
SoTrap::~SoTrap() {
 delete children;
}


// initClass
void SoTrap::initClass(){
  // This statement is required.
  SO_NODE_INIT_CLASS(SoTrap,SoShape,"Shape");
}


// generatePrimitives
void SoTrap::generatePrimitives(SoAction *action) {
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
  float TthetaCphi = FTAN(pTheta.getValue())*FCOS(pPhi.getValue());
  float TthetaSphi = FTAN(pTheta.getValue())*FSIN(pPhi.getValue());
  float Talp1 = FTAN(pAlp1.getValue()); 
  float Talp2 = FTAN(pAlp2.getValue());

  float points[NPOINTS][3];
  points[0][0] =  pDx2.getValue()+pDy1.getValue()*Talp1;
  points[0][1] =  pDy1.getValue();
  points[0][2] = -pDz.getValue();

  points[1][0] = -pDx2.getValue()+pDy1.getValue()*Talp1;
  points[1][1] =  pDy1.getValue();
  points[1][2] = -pDz.getValue();

  points[2][0] = -pDx1.getValue()-pDy1.getValue()*Talp1;
  points[2][1] = -pDy1.getValue();
  points[2][2] = -pDz.getValue();

  points[3][0] =  pDx1.getValue()-pDy1.getValue()*Talp1;
  points[3][1] = -pDy1.getValue();
  points[3][2] = -pDz.getValue();

  points[4][0] =  pDx4.getValue()+pDy2.getValue()*Talp2;
  points[4][1] =  pDy2.getValue();
  points[4][2] =  pDz.getValue();

  points[5][0] = -pDx4.getValue()+pDy2.getValue()*Talp2;
  points[5][1] =  pDy2.getValue();
  points[5][2] =  pDz.getValue();

  points[6][0] = -pDx3.getValue()-pDy2.getValue()*Talp2;
  points[6][1] = -pDy2.getValue();
  points[6][2] =  pDz.getValue();

  points[7][0] =  pDx3.getValue()-pDy2.getValue()*Talp2;
  points[7][1] = -pDy2.getValue();
  points[7][2] =  pDz.getValue();

  int i;
  for (i=0;i<4;i++) {
    points[i][0] -= pDz.getValue()*TthetaCphi;
    points[i][1] -= pDz.getValue()*TthetaSphi;
  }
  for (i=4;i<8;i++) {
    points[i][0] += pDz.getValue()*TthetaCphi;
    points[i][1] += pDz.getValue()*TthetaSphi;
  }

  SbVec3f normals[NFACES];
  int nf;
  for (nf=0;nf<NFACES;nf++) {
    int j0 = indices[5*nf + 0];
    int j1 = indices[5*nf + 1];
    int j2 = indices[5*nf + 2];
    SbVec3f p0(points[j0][0],points[j0][1],points[j0][2]); 
    SbVec3f p1(points[j1][0],points[j1][1],points[j1][2]); 
    SbVec3f p2(points[j2][0],points[j2][1],points[j2][2]); 
    normals[nf] = (p1-p0).cross(p2-p0);
    normals[nf].normalize();
  }

  float x,y,z;
  int   index;
  for (nf=0;nf<NFACES;nf++) {
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
SoChildList *SoTrap::getChildren() const {
  return children;
}


// computeBBox
void SoTrap::computeBBox(SoAction *, SbBox3f &box, SbVec3f &center ){
  float pDx= pDx1.getValue(),pDy=pDy1.getValue();
  
  if (pDx2.getValue() > pDx) pDx = pDx2.getValue(); 
  if (pDx3.getValue() > pDx) pDx = pDx3.getValue(); 
  if (pDx4.getValue() > pDx) pDx = pDx4.getValue(); 
  if (pDy2.getValue() > pDy) pDy = pDy2.getValue(); 
  float TthetaCphi = FTAN(pTheta.getValue())*FCOS(pPhi.getValue());
  float TthetaSphi = FTAN(pTheta.getValue())*FSIN(pPhi.getValue());
  float Xalp = FFABS(std::tan(pAlp1.getValue())*pDy1.getValue());
  float Xalp2 = FFABS(std::tan(pAlp2.getValue())*pDy2.getValue());
  if (Xalp< Xalp2) Xalp=Xalp2;
  pDx += FFABS(TthetaCphi*pDz.getValue());
  pDx += Xalp;
  pDy += FFABS(TthetaSphi*pDz.getValue());  


  center.setValue(0,0,0);
  box.setBounds(SbVec3f(-pDx,-pDy,-pDz.getValue()),
	  SbVec3f( pDx, pDy, pDz.getValue()));
}




// updateChildren
void SoTrap::updateChildren() {


  // Redraw the G4Trap....

  assert(children->getLength()==1);
  SoSeparator       *sep                = (SoSeparator *)  ( *children)[0];
  SoCoordinate3     *theCoordinates     = (SoCoordinate3 *)      ( sep->getChild(0));
  SoNormal          *theNormals         = (SoNormal *)           ( sep->getChild(1)); 
  SoNormalBinding   *theNormalBinding   = (SoNormalBinding *)    ( sep->getChild(2));
  SoIndexedFaceSet  *theFaceSet         = (SoIndexedFaceSet *)   ( sep->getChild(3));

  const int NPOINTS=8, NFACES=6, NINDICES = NFACES*5;
  float points[NPOINTS][3];
  // Indices for the eight faces
#ifdef INVENTOR2_0
  static long  
#else 
  static int32_t
#endif
  indices[NINDICES] = {3,2,1,0, SO_END_FACE_INDEX, // bottom
                     4,5,6,7, SO_END_FACE_INDEX, // top
                     0,1,5,4, SO_END_FACE_INDEX, 
                     1,2,6,5, SO_END_FACE_INDEX,
                     2,3,7,6, SO_END_FACE_INDEX,
                     3,0,4,7, SO_END_FACE_INDEX};

  
  // points for the eight vertices
  float TthetaCphi = FTAN(pTheta.getValue())*FCOS(pPhi.getValue());
  float TthetaSphi = FTAN(pTheta.getValue())*FSIN(pPhi.getValue());
  float Talp1 = FTAN(pAlp1.getValue()); 
  float Talp2 = FTAN(pAlp2.getValue());

  points[0][0] =  pDx2.getValue()+pDy1.getValue()*Talp1;
  points[0][1] =  pDy1.getValue();
  points[0][2] = -pDz.getValue();

  points[1][0] = -pDx2.getValue()+pDy1.getValue()*Talp1;
  points[1][1] =  pDy1.getValue();
  points[1][2] = -pDz.getValue();

  points[2][0] = -pDx1.getValue()-pDy1.getValue()*Talp1;
  points[2][1] = -pDy1.getValue();
  points[2][2] = -pDz.getValue();

  points[3][0] =  pDx1.getValue()-pDy1.getValue()*Talp1;
  points[3][1] = -pDy1.getValue();
  points[3][2] = -pDz.getValue();

  points[4][0] =  pDx4.getValue()+pDy2.getValue()*Talp2;
  points[4][1] =  pDy2.getValue();
  points[4][2] =  pDz.getValue();

  points[5][0] = -pDx4.getValue()+pDy2.getValue()*Talp2;
  points[5][1] =  pDy2.getValue();
  points[5][2] =  pDz.getValue();

  points[6][0] = -pDx3.getValue()-pDy2.getValue()*Talp2;
  points[6][1] = -pDy2.getValue();
  points[6][2] =  pDz.getValue();

  points[7][0] =  pDx3.getValue()-pDy2.getValue()*Talp2;
  points[7][1] = -pDy2.getValue();
  points[7][2] =  pDz.getValue();

  int i;
  for (i=0;i<4;i++) {
    points[i][0] -= pDz.getValue()*TthetaCphi;
    points[i][1] -= pDz.getValue()*TthetaSphi;
  }
  for (i=4;i<8;i++) {
    points[i][0] += pDz.getValue()*TthetaCphi;
    points[i][1] += pDz.getValue()*TthetaSphi;
  }

  for (int np=0;np<NPOINTS;np++) theCoordinates->point.set1Value(np,points[np][0],points[np][1],points[np][2]);
  theFaceSet->coordIndex.setValues(0,NINDICES,indices);
  theNormals->vector.deleteValues(0);
  theNormals->vector.insertSpace(0,6);
  for (int n=0;n<6;n++) {
    int i0 = 5*n+0,i1=5*n+1,i2=5*n+2;
    int j0 = theFaceSet->coordIndex[i0];
    int j1 = theFaceSet->coordIndex[i1];
    int j2 = theFaceSet->coordIndex[i2];
    SbVec3f p0= theCoordinates->point[j0]; 
    SbVec3f p1= theCoordinates->point[j1]; 
    SbVec3f p2= theCoordinates->point[j2]; 
    SbVec3f normal = (p1-p0).cross(p2-p0);
    normal.normalize();
    theNormals->vector.set1Value(n,normal);
  }
  theNormalBinding->value=SoNormalBinding::PER_FACE;
}

// generateChildren
void SoTrap::generateChildren() {

  // This routines creates one SoSeparator, one SoCoordinate3, and
  // one SoLineSet, and puts it in the child list.  This is done only
  // once, whereas redrawing the position of the coordinates occurs each
  // time an update is necessary, in the updateChildren routine. 

  assert(children->getLength() ==0);
  SoSeparator      *sep              = new SoSeparator(); 
  SoCoordinate3    *theCoordinates   = new SoCoordinate3();
  SoNormal         *theNormals       = new SoNormal(); 
  SoNormalBinding  *theNormalBinding = new SoNormalBinding();
  SoIndexedFaceSet *theFaceSet       = new SoIndexedFaceSet();
  // 
  // This line costs some in render quality! but gives speed.
  // 
  sep->addChild(theCoordinates);
  sep->addChild(theNormals);
  sep->addChild(theNormalBinding);
  sep->addChild(theFaceSet);
  children->append(sep);
}

// generateAlternateRep
void SoTrap::generateAlternateRep() {

  // This routine sets the alternate representation to the child
  // list of this mode.  

  if (children->getLength() == 0) generateChildren();
  updateChildren();
  alternateRep.setValue((SoSeparator *)  ( *children)[0]);
}

// clearAlternateRep
void SoTrap::clearAlternateRep() {
  alternateRep.setValue(NULL);
}

#endif
