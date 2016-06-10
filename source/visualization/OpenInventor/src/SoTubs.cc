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
// $Id: SoTubs.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
/*-----------------------------HEPVis---------------------------------------*/
/*                                                                          */
/* Node:             SoTubs                                                 */
/* Description:      Represents the G4Tubs Geant Geometry entity            */
/* Author:           Joe Boudreau Nov 11 1996                               */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#ifdef G4VIS_BUILD_OI_DRIVER

// this :
#include "HEPVis/nodes/SoTubs.h"

#include <assert.h>
#include <cmath>

#include <Inventor/SbBox.h>
#include <Inventor/actions/SoGLRenderAction.h>
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
SO_NODE_SOURCE(SoTubs)

// Constructor
SoTubs::SoTubs() {


  // This statement is required
  SO_NODE_CONSTRUCTOR(SoTubs);

  // Data fields are initialized like this:
  SO_NODE_ADD_FIELD(pRMin,               (0));
  SO_NODE_ADD_FIELD(pRMax,               (1));
  SO_NODE_ADD_FIELD(pDz,                 (10));
  SO_NODE_ADD_FIELD(pSPhi,               (0));
  SO_NODE_ADD_FIELD(pDPhi,               ((float)(2*M_PI)));
  SO_NODE_ADD_FIELD(alternateRep,        (NULL));
  children = new SoChildList(this);
}

// Destructor
SoTubs::~SoTubs() {
  delete children;
}


// initClass
void SoTubs::initClass(){
  // This statement is required.
  SO_NODE_INIT_CLASS(SoTubs,SoShape,"Shape");
}

// generatePrimitives
void SoTubs::generatePrimitives(SoAction *action) {
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
  const SoTextureCoordinateElement* tce = NULL;
  SbVec4f texCoord;
  if (useTexFunction) {
    tce = SoTextureCoordinateElement::getInstance(state);
  }
  else {
    texCoord[2] = 0.0;
    texCoord[3] = 1.0;
  }
  SbVec3f point, normal;


  ///////////////////////////////////////////////////////
  //-----------------------------------------------------
#define GEN_VERTEX(pv,x,y,z,s,t,nx,ny,nz)               \
  point.setValue((float)(x),(float)(y),(float)(z));     \
  normal.setValue((float)(nx),(float)(ny),(float)(nz)); \
  if (useTexFunction) {                                 \
    texCoord=tce->get(point,normal);                    \
  } else {                                              \
    texCoord[0]=(float)(s);                             \
    texCoord[1]=(float)(t);                             \
  }                                                     \
  pv.setPoint(point);                                   \
  pv.setNormal(normal);                                 \
  pv.setTextureCoords(texCoord);                        \
  shapeVertex(&pv);
  //-----------------------------------------------------
  ///////////////////////////////////////////////////////

  int NPHI = (int)(2+22*std::fabs(pDPhi.getValue()/(2.0*M_PI)));
  double deltaPhi = pDPhi.getValue()/NPHI, phi0 = pSPhi.getValue(),phi1=phi0+pDPhi.getValue();
  double rMax=pRMax.getValue(),rMin=pRMin.getValue();
  double zMax=pDz.getValue(),zMin=-zMax;
  double cosPhi0=std::cos(phi0), sinPhi0=std::sin(phi0);
  double cosPhi1=std::cos(phi1), sinPhi1=std::sin(phi1);
  double cosDeltaPhi=std::cos(deltaPhi),sinDeltaPhi=std::sin(deltaPhi);
  //
  // The outer surface!
  //
  int i;
  double sinPhi,cosPhi;
  beginShape(action,TRIANGLE_STRIP);
  sinPhi=sinPhi0;
  cosPhi=cosPhi0;
  for (i = 0; i<=NPHI; i++) {
    GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi,zMax,0.0,0.0,cosPhi,sinPhi,0);   
    GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi,zMin,1.0,1.0,cosPhi,sinPhi,0);   
    inc(sinPhi, cosPhi, sinDeltaPhi, cosDeltaPhi);    
  }
  endShape();
  //
  // The inner surface!
  //
  if(rMin!=0.F) {
    beginShape(action,TRIANGLE_STRIP);
    sinPhi=sinPhi0;
    cosPhi=cosPhi0;
    for (i = 0; i<=NPHI; i++) {
      GEN_VERTEX(pv,rMin*cosPhi,rMin*sinPhi,zMax,0.0,0.0,-cosPhi,-sinPhi,0);   
      GEN_VERTEX(pv,rMin*cosPhi,rMin*sinPhi,zMin,1.0,1.0,-cosPhi,-sinPhi,0);   
      inc(sinPhi, cosPhi, sinDeltaPhi, cosDeltaPhi);    
    } 
    endShape();
  }
  if (std::fabs(deltaPhi)<2.0*M_PI) { 
    //
    // The end 
    //
    beginShape(action,TRIANGLE_STRIP);
    sinPhi=sinPhi0;
    cosPhi=cosPhi0;
    GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi,zMax,0.0,0.0,sinPhi,-cosPhi,0);   
    GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi,zMin,1.0,1.0,sinPhi,-cosPhi,0);   
    GEN_VERTEX(pv,rMin*cosPhi,rMin*sinPhi,zMax,1.0,0.0,sinPhi,-cosPhi,0);   
    GEN_VERTEX(pv,rMin*cosPhi,rMin*sinPhi,zMin,0.0,1.0,sinPhi,-cosPhi,0);   
    endShape();
    //
    // The other end 
    //
    beginShape(action,TRIANGLE_STRIP);
    sinPhi=sinPhi1;
    cosPhi=cosPhi1;
    GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi, zMax,0.0,0.0,-sinPhi,+cosPhi,0);   
    GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi, zMin,1.0,1.0,-sinPhi,+cosPhi,0);   
    GEN_VERTEX(pv,rMin*cosPhi,rMin*sinPhi, zMax,1.0,0.0,-sinPhi,+cosPhi,0);   
    GEN_VERTEX(pv,rMin*cosPhi,rMin*sinPhi, zMin,0.0,1.0,-sinPhi,+cosPhi,0);   
    endShape();
  }
  //
  // The outer surface at z=+PDZ
  //
  if(rMin==0.F) {
    beginShape(action,TRIANGLE_FAN);
    sinPhi=sinPhi0;
    cosPhi=cosPhi0;
    GEN_VERTEX(pv,0,0,zMax,0.0,0.0,0,0,1);   
    for (i = 0; i<=NPHI; i++) {
      GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi,zMax,1.0,1.0,0,0,1);   
      inc(sinPhi, cosPhi, sinDeltaPhi, cosDeltaPhi);    
    }
    endShape();
    //
    // The outer surface at z=-PDZ
    //
    beginShape(action,TRIANGLE_FAN);
    sinPhi=sinPhi0;
    cosPhi=cosPhi0;
    GEN_VERTEX(pv,0,0,zMin,0.0,0.0,0,0,-1);   
    for (i = 0; i<=NPHI; i++) {
      GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi,zMin,1.0,1.0,0,0,-1);   
      inc(sinPhi, cosPhi, sinDeltaPhi, cosDeltaPhi);    
    }
    endShape();
  } else {
    beginShape(action,TRIANGLE_STRIP);
    sinPhi=sinPhi0;
    cosPhi=cosPhi0;
    for (i = 0; i<=NPHI; i++) {
      GEN_VERTEX(pv,rMin*cosPhi,rMin*sinPhi,zMax,0.0,0.0,0,0,1);   
      GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi,zMax,1.0,1.0,0,0,1);   
      inc(sinPhi, cosPhi, sinDeltaPhi, cosDeltaPhi);    
    }
    endShape();
    //
    // The outer surface at z=-PDZ
    //
    beginShape(action,TRIANGLE_STRIP);
    sinPhi=sinPhi0;
    cosPhi=cosPhi0;
    for (i = 0; i<=NPHI; i++) {
      GEN_VERTEX(pv,rMin*cosPhi,rMin*sinPhi,zMin,0.0,0.0,0,0,-1);   
      GEN_VERTEX(pv,rMax*cosPhi,rMax*sinPhi,zMin,1.0,1.0,0,0,-1);   
      inc(sinPhi, cosPhi, sinDeltaPhi, cosDeltaPhi);    
    }
    endShape();
  }
}

// getChildren
SoChildList *SoTubs::getChildren() const {
  return children;
}


// computeBBox
void SoTubs::computeBBox(SoAction *, SbBox3f &box, SbVec3f &center ){
  SbVec3f vmin(-pRMax.getValue(),-pRMax.getValue(),-pDz.getValue()), 
          vmax( pRMax.getValue(), pRMax.getValue(), pDz.getValue());
  center.setValue(0,0,0);
  box.setBounds(vmin,vmax);
}


// updateChildren
void SoTubs::updateChildren() {

  // Redraw the G4Tubs....

  assert(children->getLength()==1);
  SoSeparator       *sep                = (SoSeparator *)  ( *children)[0];
  SoCoordinate3     *theCoordinates     = (SoCoordinate3 *)      ( sep->getChild(0));
  SoNormal          *theNormals         = (SoNormal *)           ( sep->getChild(1)); 
  SoNormalBinding   *theNormalBinding   = (SoNormalBinding *)    ( sep->getChild(2));
  SoIndexedFaceSet  *theFaceSet         = (SoIndexedFaceSet *)   ( sep->getChild(3));
  
  
  const int NPHI=24, NPOINTS=2*(2*NPHI+2), NFACES=4*NPHI+2, NINDICES = NFACES*5;
  float points[NPOINTS][3],normals[NFACES][3];
#ifdef INVENTOR2_0
  static long     indices[NINDICES];
#else
  static int32_t  indices[NINDICES];
#endif
    
  static int init=0;
  double phi, pp, DeltaPhi;
    
  // Indices need to be generated once! This is here to keep it close to the point
  // generation, since otherwise it will be confusing.
    
  int i;
  if (!init) {
    init = 1;
    // Outer face
    for (i = 0; i< NPHI; i++) {
      // 0 1 3 2;
      indices[5*i+0] =  2*i+0;
      indices[5*i+1] =  2*i+1;
      indices[5*i+2] =  2*i+3;
      indices[5*i+3] =  2*i+2;
      indices[5*i+4] = SO_END_FACE_INDEX;
    }
    // the inner face
    for (i=0;i<NPHI;i++) {
      indices[5*1*NPHI + 5*i+0] = 2*NPHI+2 + 2*i+0;  
      indices[5*1*NPHI + 5*i+1] = 2*NPHI+2 + 2*i+1;
      indices[5*1*NPHI + 5*i+2] = 2*NPHI+2 + 2*i+3;
      indices[5*1*NPHI + 5*i+3] = 2*NPHI+2 + 2*i+2;
      indices[5*1*NPHI + 5*i+4] = SO_END_FACE_INDEX;
    }
    // the top side
    for (i=0;i<NPHI;i++) {
      indices[5*2*NPHI + 5*i+0] = 2*i+0;
      indices[5*2*NPHI + 5*i+1] = 2*i+2;
      indices[5*2*NPHI + 5*i+2] = NPOINTS - (2*i+4);
      indices[5*2*NPHI + 5*i+3] = NPOINTS - (2*i+2);
      indices[5*2*NPHI + 5*i+4] = SO_END_FACE_INDEX;
    }
    // the bottom side
    for (i=0;i<NPHI;i++) {
      indices[5*3*NPHI + 5*i+0] = 2*i+1;
      indices[5*3*NPHI + 5*i+1] = NPOINTS - (2*i+1);
      indices[5*3*NPHI + 5*i+2] = NPOINTS - (2*i+3);
      indices[5*3*NPHI + 5*i+3] = 2*i+3;
      indices[5*3*NPHI + 5*i+4] = SO_END_FACE_INDEX;
    }
    // the odd side
    indices[5*4*NPHI +0] = 2*NPHI;
    indices[5*4*NPHI +1] = 2*NPHI+1;
    indices[5*4*NPHI +2] = 2*NPHI+3;
    indices[5*4*NPHI +3] = 2*NPHI+2;
    indices[5*4*NPHI +4] = SO_END_FACE_INDEX;
    // aother odd side
    indices[5*4*NPHI +5 +0] = 0;
    indices[5*4*NPHI +5 +1] = NPOINTS-2;
    indices[5*4*NPHI +5 +2] = NPOINTS-1;
    indices[5*4*NPHI +5 +3] = 1;
    indices[5*4*NPHI +5 +4] = SO_END_FACE_INDEX;
  }
  // Points need to be generated each time:
  if (pDPhi.getValue()<2*M_PI) {
    // the odd side
    indices[5*4*NPHI +0] = 2*NPHI;
    indices[5*4*NPHI +1] = 2*NPHI+1;
    indices[5*4*NPHI +2] = 2*NPHI+3;
    indices[5*4*NPHI +3] = 2*NPHI+2;
    indices[5*4*NPHI +4] = SO_END_FACE_INDEX;
    // aother odd side
    indices[5*4*NPHI +5 +0] = 0;
    indices[5*4*NPHI +5 +1] = NPOINTS-2;
    indices[5*4*NPHI +5 +2] = NPOINTS-1;
    indices[5*4*NPHI +5 +3] = 1;
    indices[5*4*NPHI +5 +4] = SO_END_FACE_INDEX;
  } 
  else {
    // the odd side
    indices[5*4*NPHI +0] = SO_END_FACE_INDEX;
    indices[5*4*NPHI +1] = SO_END_FACE_INDEX;
    indices[5*4*NPHI +2] = SO_END_FACE_INDEX;
    indices[5*4*NPHI +3] = SO_END_FACE_INDEX;
    indices[5*4*NPHI +4] = SO_END_FACE_INDEX;
    // aother odd side
    indices[5*4*NPHI +5 +0] = SO_END_FACE_INDEX;
    indices[5*4*NPHI +5 +1] = SO_END_FACE_INDEX;
    indices[5*4*NPHI +5 +2] = SO_END_FACE_INDEX;
    indices[5*4*NPHI +5 +3] = SO_END_FACE_INDEX;
    indices[5*4*NPHI +5 +4] = SO_END_FACE_INDEX;
  }
  // The outer surface
  DeltaPhi = pDPhi.getValue()/NPHI, phi = pSPhi.getValue();
  for (i = 0; i<=NPHI; i++) {
    points[2*i+0][0] = pRMax.getValue()*FCOS(phi); 
    points[2*i+0][1]= pRMax.getValue()*FSIN(phi); 
    points[2*i+0][2] = +pDz.getValue();

    points[2*i+1][0] = pRMax.getValue()*FCOS(phi); 
    points[2*i+1][1]= pRMax.getValue()*FSIN(phi); 
    points[2*i+1][2] = -pDz.getValue();

    pp = phi+DeltaPhi/2.0;
    if (i!=NPHI) {
      normals[i][0] = FCOS(pp); 
      normals[i][1] = FSIN(pp); 
      normals[i][2] = 0;
    }
    phi+=DeltaPhi;
  }
  // The inner surface
  phi = pSPhi.getValue() + pDPhi.getValue();
  for (i = 0; i<=NPHI; i++) {
    points[2*NPHI+2+2*i+0][0] = pRMin.getValue()*FCOS(phi); 
    points[2*NPHI+2+2*i+0][1] = pRMin.getValue()*FSIN(phi); 
    points[2*NPHI+2+2*i+0][2] = +pDz.getValue();
    points[2*NPHI+2+2*i+1][0] = pRMin.getValue()*FCOS(phi); 
    points[2*NPHI+2+2*i+1][1] = pRMin.getValue()*FSIN(phi); 
    points[2*NPHI+2+2*i+1][2] = -pDz.getValue();
    pp = phi-DeltaPhi/2.0;
    if (i!=NPHI) {
      normals[NPHI+i][0] = -FCOS(pp); 
      normals[NPHI+i][1] = -FSIN(pp); 
      normals[NPHI+i][2] = 0;
    }
    phi-=DeltaPhi;
  }
  // The top side
  for (i=0;i<NPHI;i++) {
    normals[2*NPHI+i][0]=normals[2*NPHI+i][1]=0; 
    normals[2*NPHI+i][2]=  1.0;
  } 
  // The bottom side
  for (i=0;i<NPHI;i++) {
    normals[3*NPHI+i][0]=normals[3*NPHI+i][1]=0; 
    normals[3*NPHI+i][2]= -1.0;
  } 
  // The odd side
  phi = pSPhi.getValue(); 
  normals[4*NPHI+0][0]=  FSIN(phi); 
  normals[4*NPHI+0][1]= -FCOS(phi); 
  normals[4*NPHI+0][2]=0;
    
    // Another odd side
  phi = pSPhi.getValue()+pDPhi.getValue(); 
  normals[4*NPHI+1][0]= -FSIN(phi); 
  normals[4*NPHI+1][1]= +FCOS(phi); 
  normals[4*NPHI+1][2]=0;
    
  for (int np=0;np<NPOINTS; np++) theCoordinates->point.set1Value(np,points[np][0],points[np][1],points[np][2]);
  for (int ni=0;ni<NINDICES;ni++) theFaceSet->coordIndex.set1Value(ni,indices[ni]);
  for (int nf=0;nf<NFACES;nf++) theNormals->vector.set1Value(nf,normals[nf][0],normals[nf][1],normals[nf][2]);
  theNormalBinding->value=SoNormalBinding::PER_FACE;
}

// generateChildren
void SoTubs::generateChildren() {

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
void SoTubs::generateAlternateRep() {

  // This routine sets the alternate representation to the child
  // list of this mode.  

  if (children->getLength() == 0) generateChildren();
  updateChildren();
  alternateRep.setValue((SoSeparator *)  ( *children)[0]);
}

// clearAlternateRep
void SoTubs::clearAlternateRep() {
  alternateRep.setValue(NULL);
}

#endif
