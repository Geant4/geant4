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
/* Node:             SoPolyhedron                                           */
/* Description:      SoNode to represent HepPolyhedron                      */
/* Author:           Guy Barrand                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

// this :
#include "Geant4_SoPolyhedron.h"

#include <Inventor/SbBox.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/elements/SoTextureCoordinateElement.h>

//#include <HEPVis/SbMath.h>
#define SbMinimum(a,b) ((a)<(b)?a:b)
#define SbMaximum(a,b) ((a)>(b)?a:b)

#include "HepPolyhedron.h"

//typedef SbVec3f HVPoint3D;
//typedef SbVec3f HVNormal3D;

typedef HepPoint3D HVPoint3D;
typedef HepNormal3D HVNormal3D;

SO_NODE_SOURCE(Geant4_SoPolyhedron) 
//////////////////////////////////////////////////////////////////////////////
void Geant4_SoPolyhedron::initClass(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_INIT_CLASS(Geant4_SoPolyhedron,SoShape,"Shape");
}
//////////////////////////////////////////////////////////////////////////////
Geant4_SoPolyhedron::Geant4_SoPolyhedron(
)
:fPolyhedron(0)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_CONSTRUCTOR(Geant4_SoPolyhedron);
  SO_NODE_ADD_FIELD(solid,(TRUE));
  SO_NODE_ADD_FIELD(reducedWireFrame,(TRUE));
}
//////////////////////////////////////////////////////////////////////////////
Geant4_SoPolyhedron::Geant4_SoPolyhedron(
 const HepPolyhedron& aPolyhedron
)
:fPolyhedron(0)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_CONSTRUCTOR(Geant4_SoPolyhedron);
  SO_NODE_ADD_FIELD(solid,(TRUE));
  fPolyhedron = new HepPolyhedron(aPolyhedron);
}
//////////////////////////////////////////////////////////////////////////////
Geant4_SoPolyhedron::Geant4_SoPolyhedron(
 HepPolyhedron* aPolyhedron
)
:fPolyhedron(aPolyhedron)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_CONSTRUCTOR(Geant4_SoPolyhedron);
  SO_NODE_ADD_FIELD(solid,(TRUE));
}
//////////////////////////////////////////////////////////////////////////////
Geant4_SoPolyhedron::~Geant4_SoPolyhedron(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  delete fPolyhedron;
}
//////////////////////////////////////////////////////////////////////////////
void Geant4_SoPolyhedron::generatePrimitives(
 SoAction* aAction
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fPolyhedron) return;
  if(fPolyhedron->GetNoFacets()<=0) return; // Abnormal polyhedron.

  SoState *state = aAction->getState();
  SbBool useTexFunction =
    (SoTextureCoordinateElement::getType(state) == 
     SoTextureCoordinateElement::FUNCTION);
  const SoTextureCoordinateElement *tce = NULL;
  SbVec4f texCoord;
  if (useTexFunction) {
    tce = SoTextureCoordinateElement::getInstance(state);
  } else {
    texCoord[2] = 0.0;
    texCoord[3] = 1.0;
  }

  if(solid.getValue()==TRUE) {
    SoPrimitiveVertex pv;
    SbVec3f point, normal;
    //////////////////////////////////////////
    //----------------------------------------
#define GEN_VERTEX(pv,x,y,z,s,t,nx,ny,nz)  \
  point.setValue(x,y,z);                   \
  normal.setValue(nx,ny,nz);               \
  if (useTexFunction) {                    \
    texCoord=tce->get(point,normal);       \
  } else {                                 \
    texCoord[0]=s;                         \
    texCoord[1]=t;                         \
  }                                        \
  pv.setPoint(point);                      \
  pv.setNormal(normal);                    \
  pv.setTextureCoords(texCoord);           \
  shapeVertex(&pv);
  //----------------------------------------
  //////////////////////////////////////////

    // Assume all facets are convex quadrilaterals :
    bool notLastFace;
    do {
      HVNormal3D unitNormal;
      notLastFace = fPolyhedron->GetNextUnitNormal(unitNormal);

      beginShape(aAction,POLYGON);
      bool notLastEdge;
      int edgeFlag = 1;
      do {
        HVPoint3D vertex;
        notLastEdge = fPolyhedron->GetNextVertex(vertex,edgeFlag);
        GEN_VERTEX(pv,
                   vertex[0],
                   vertex[1],
                   vertex[2],
                   0.0,0.0,
                   unitNormal[0],
                   unitNormal[1],
                   unitNormal[2]);
      } while (notLastEdge);
      endShape();
    } while (notLastFace);
  } else {
    SoPrimitiveVertex pvb,pve;
    pve.setTextureCoords(texCoord);
    pvb.setTextureCoords(texCoord);

#ifdef __COIN__ // To bypass a bug in invokeLineSegment when picking.
    beginShape(aAction,POLYGON);
    endShape();
#endif

    SbVec3f point;
    bool notLastFace;
    do {
      HVNormal3D unitNormal;
      notLastFace = fPolyhedron->GetNextUnitNormal(unitNormal);

      SbVec3f normal;
      normal.setValue(unitNormal[0],unitNormal[1],unitNormal[2]);

      // Treat edges :
      int edgeFlag = 1;
      int prevEdgeFlag = edgeFlag;
      bool notLastEdge;
      SbBool firstEdge = TRUE;
      do {
        HVPoint3D vertex;
        notLastEdge = fPolyhedron->GetNextVertex(vertex,edgeFlag);
        if(reducedWireFrame.getValue()==FALSE) edgeFlag = 1;        
        if(firstEdge) {
          if(edgeFlag) {
            pvb.setNormal(normal);
            point.setValue(vertex[0],vertex[1],vertex[2]);
            pvb.setPoint(point);
          } else {
          }
          firstEdge = FALSE;
          prevEdgeFlag = edgeFlag;
        } else {
          if(edgeFlag!=prevEdgeFlag) { 
            if(edgeFlag) { // Pass to a visible edge :
              pvb.setNormal(normal);
              point.setValue(vertex[0],vertex[1],vertex[2]);
              pvb.setPoint(point);
            } else { // Pass to an invisible edge :
              pve.setNormal(normal);
              point.setValue(vertex[0],vertex[1],vertex[2]);
              pve.setPoint(point);
              invokeLineSegmentCallbacks(aAction,&pvb,&pve);
            }
            prevEdgeFlag = edgeFlag;
          } else {
            if(edgeFlag) {
              pve.setNormal(normal);
              point.setValue(vertex[0],vertex[1],vertex[2]);
              pve.setPoint(point);
              invokeLineSegmentCallbacks(aAction,&pvb,&pve);
              pvb = pve;
            } else {
            }
          }
        }
      } while (notLastEdge);
    } while (notLastFace);
  }

}
//////////////////////////////////////////////////////////////////////////////
void Geant4_SoPolyhedron::computeBBox(
 SoAction*
,SbBox3f& aBox
,SbVec3f& aCenter
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fPolyhedron) return;
  if(fPolyhedron->GetNoFacets()<=0) { // Abnormal polyhedron.
    SbVec3f vmin(-1,-1,-1);
    SbVec3f vmax( 1, 1, 1);
    aBox.setBounds(vmin,vmax);
    aCenter.setValue(0,0,0);
  } else {
    SbBool first = TRUE;
    float xmn = 0,ymn = 0,zmn = 0;
    float xmx = 0,ymx = 0,zmx = 0;
    float xct = 0,yct = 0,zct = 0;
    SbVec3f point;
    int count = 0;
    // Assume all facets are convex quadrilaterals :
    bool notLastFace;
    do {
      HVNormal3D unitNormal;
      notLastFace = fPolyhedron->GetNextUnitNormal(unitNormal);
      bool notLastEdge;
      do {
        HVPoint3D vertex;
        int edgeFlag = 1;
        notLastEdge = fPolyhedron->GetNextVertex(vertex,edgeFlag);
        point.setValue(vertex[0],vertex[1],vertex[2]);
        if(first==TRUE) {
          xct = xmx = xmn = point[0];
          yct = ymx = ymn = point[1];
          zct = zmx = zmn = point[2];
          count++;
          first = FALSE;
        } else {
          xmn = SbMinimum(xmn,point[0]);
          ymn = SbMinimum(ymn,point[1]);
          zmn = SbMinimum(zmn,point[2]);
          //
          xmx = SbMaximum(xmx,point[0]);
          ymx = SbMaximum(ymx,point[1]);
          zmx = SbMaximum(zmx,point[2]);
          //
          xct += point[0];
          yct += point[1];
          zct += point[2];
          count++;
        }
        //
      } while (notLastEdge);
    } while (notLastFace);
    SbVec3f vmin(xmn,ymn,zmn);
    SbVec3f vmax(xmx,ymx,zmx);
    aBox.setBounds(vmin,vmax);
    if(count==0)
      aCenter.setValue(0,0,0);
    else
      aCenter.setValue(xct/count,yct/count,zct/count);
  }
}


#endif
