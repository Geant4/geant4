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
#include <Inventor/nodes/SoSeparator.h>

//#include <HEPVis/SbMath.h>
#define SbMinimum(a,b) ((a)<(b)?a:b)
#define SbMaximum(a,b) ((a)>(b)?a:b)

#include <HEPVis/actions/SoAlternateRepAction.h>

#include "G4Polyhedron.hh"

//typedef SbVec3f HVPoint3D;
//typedef SbVec3f HVNormal3D;

typedef HepGeom::Point3D<double> HVPoint3D;
typedef HepGeom::Normal3D<double> HVNormal3D;

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
  SO_NODE_ADD_FIELD(alternateRep,(NULL));
}
//////////////////////////////////////////////////////////////////////////////
Geant4_SoPolyhedron::Geant4_SoPolyhedron(
 const G4Polyhedron& aPolyhedron
)
:fPolyhedron(0)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_CONSTRUCTOR(Geant4_SoPolyhedron);
  SO_NODE_ADD_FIELD(solid,(TRUE));
  SO_NODE_ADD_FIELD(reducedWireFrame,(TRUE));
  SO_NODE_ADD_FIELD(alternateRep,(NULL));

  fPolyhedron = new G4Polyhedron(aPolyhedron);
}
//////////////////////////////////////////////////////////////////////////////
Geant4_SoPolyhedron::Geant4_SoPolyhedron(
 G4Polyhedron* aPolyhedron
)
:fPolyhedron(aPolyhedron)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_CONSTRUCTOR(Geant4_SoPolyhedron);
  SO_NODE_ADD_FIELD(solid,(TRUE));
  SO_NODE_ADD_FIELD(reducedWireFrame,(TRUE));
  SO_NODE_ADD_FIELD(alternateRep,(NULL));
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
  SbVec4f texCoord(0.,0.,0.,0.);
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
          if(edgeFlag > 0) {
            pvb.setNormal(normal);
            point.setValue(vertex[0],vertex[1],vertex[2]);
            pvb.setPoint(point);
          } else {
          }
          firstEdge = FALSE;
          prevEdgeFlag = edgeFlag;
        } else {
          if(edgeFlag!=prevEdgeFlag) { 
            if(edgeFlag > 0) { // Pass to a visible edge :
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
            if(edgeFlag > 0) {
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

#include <Inventor/nodes/SoNormalBinding.h>
#include <Inventor/nodes/SoNormal.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
//////////////////////////////////////////////////////////////////////////////
void Geant4_SoPolyhedron::generateAlternateRep(
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fPolyhedron) return;
  if(fPolyhedron->GetNoFacets()<=0) return; // Abnormal polyhedron.
  if(fPolyhedron->GetNoVertices()<=0) return; // Abnormal polyhedron.

  if(solid.getValue()==TRUE) {

    SoSeparator* separator = new SoSeparator;

    SoNormalBinding* normalBinding = new SoNormalBinding;
    normalBinding->value = SoNormalBinding::PER_FACE;
    separator->addChild(normalBinding);

    SoCoordinate3* coordinate3 = new SoCoordinate3;
    separator->addChild(coordinate3);
    SoNormal* normal = new SoNormal;
    separator->addChild(normal);
    SoIndexedFaceSet* indexedFaceSet = new SoIndexedFaceSet;
    separator->addChild(indexedFaceSet);

    int nvert = fPolyhedron->GetNoVertices();
    int nface = fPolyhedron->GetNoFacets();

    SbVec3f* normals = new SbVec3f[nface];
    //FIXME : have the exact booking.
    SbVec3f* points = new SbVec3f[nvert];
    int32_t* coords = new int32_t[nvert+1];

    int inormal = 0;
    int icoord = 0;
    int iindex = 0;

    // Assume all facets are convex quadrilaterals :
    bool notLastFace;
    do {
      HVNormal3D unitNormal;
      notLastFace = fPolyhedron->GetNextUnitNormal(unitNormal);

      // begin face POLYGON
      int ipoint = 0;

      bool notLastEdge;
      int edgeFlag = 1;
      do {
        HVPoint3D vertex;
        notLastEdge = fPolyhedron->GetNextVertex(vertex,edgeFlag);
        points[ipoint].setValue(vertex[0],vertex[1],vertex[2]);
        coords[ipoint] = icoord + ipoint;
        ipoint++;
      } while (notLastEdge);

      // end face.
      coords[ipoint] = SO_END_FACE_INDEX;
      coordinate3->point.setValues(icoord,ipoint,points);
      icoord += ipoint;

      normals[inormal].setValue(unitNormal[0],unitNormal[1],unitNormal[2]);
      inormal++;  

      indexedFaceSet->coordIndex.setValues(iindex,(ipoint+1),coords);
      iindex += ipoint+1;

    } while (notLastFace);

    normal->vector.setValues(0,inormal,normals);

    delete [] normals;
    delete [] coords;
    delete [] points;

    alternateRep.setValue(separator);

  } else {

    SoSeparator* separator = new SoSeparator;

    int nvert = fPolyhedron->GetNoVertices();

    //FIXME : have the exact booking.
    int nedge = nvert * 3;
    int npoint = nedge*2;
    SbVec3f* points = new SbVec3f[npoint];
    int ncoord = nedge*3;
    int32_t* coords = new int32_t[ncoord];

    SbVec3f pvb(0.,0.,0.), pve(0.,0.,0.);

    SbBool empty = TRUE;
    int ipoint = 0;
    int icoord = 0;

    bool notLastFace;
    do {
      HVNormal3D unitNormal;
      notLastFace = fPolyhedron->GetNextUnitNormal(unitNormal);

      //SbVec3f normal;
      //if( (fProjection==SbProjectionRZ) || (fProjection==SbProjectionZR) ) {
        //normal.setValue(0,0,1);
      //} else {
        //normal.setValue(unitNormal[0],unitNormal[1],unitNormal[2]);
      //}

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
          if(edgeFlag > 0) {
            pvb.setValue(vertex[0],vertex[1],vertex[2]);
          } else {
          }
          firstEdge = FALSE;
          prevEdgeFlag = edgeFlag;
        } else {
          if(edgeFlag!=prevEdgeFlag) { 
            if(edgeFlag > 0) { // Pass to a visible edge :
              pvb.setValue(vertex[0],vertex[1],vertex[2]);
            } else { // Pass to an invisible edge :
              pve.setValue(vertex[0],vertex[1],vertex[2]);

              if((ipoint+1)>=npoint) {
                int new_npoint = 2 * npoint;
                SbVec3f* new_points = new SbVec3f[new_npoint];
                for(int i=0;i<npoint;i++) new_points[i] = points[i];
                delete [] points;
                npoint = new_npoint;
                points = new_points;
              }

              if((icoord+2)>=ncoord) {
                int new_ncoord = 2 * ncoord;
                int32_t* new_coords = new int32_t[new_ncoord];
                for(int i=0;i<ncoord;i++) new_coords[i] = coords[i];
                delete [] coords;
                ncoord = new_ncoord;
                coords = new_coords;
              }

              points[ipoint+0] = pvb;
              points[ipoint+1] = pve;
              coords[icoord+0] = ipoint + 0;
              coords[icoord+1] = ipoint + 1;
              coords[icoord+2] = SO_END_LINE_INDEX;
              ipoint += 2;
              icoord += 3;
              empty = FALSE;
            }
            prevEdgeFlag = edgeFlag;
          } else {
            if(edgeFlag > 0) {
              pve.setValue(vertex[0],vertex[1],vertex[2]);

              if((ipoint+1)>=npoint) {
                int new_npoint = 2 * npoint;
                SbVec3f* new_points = new SbVec3f[new_npoint];
                for(int i=0;i<npoint;i++) new_points[i] = points[i];
                delete [] points;
                npoint = new_npoint;
                points = new_points;
              }

              if((icoord+2)>=ncoord) {
                int new_ncoord = 2 * ncoord;
                int32_t* new_coords = new int32_t[new_ncoord];
                for(int i=0;i<ncoord;i++) new_coords[i] = coords[i];
                delete [] coords;
                ncoord = new_ncoord;
                coords = new_coords;
              }

              points[ipoint+0] = pvb;
              points[ipoint+1] = pve;
              coords[icoord+0] = ipoint + 0;
              coords[icoord+1] = ipoint + 1;
              coords[icoord+2] = SO_END_LINE_INDEX;
              ipoint += 2;
              icoord += 3;
              empty = FALSE;

              pvb = pve;
            } else {
            }
          }
        }
      } while (notLastEdge);
    } while (notLastFace);

    SoCoordinate3* coordinate3 = new SoCoordinate3;
    coordinate3->point.setValues(0,ipoint,points);
    separator->addChild(coordinate3);

    SoIndexedLineSet* indexedLineSet = new SoIndexedLineSet;
    indexedLineSet->coordIndex.setValues(0,icoord,coords);
    separator->addChild(indexedLineSet);

    delete [] coords;
    delete [] points;

    if(empty==TRUE) {
      separator->unref();
    } else {
      alternateRep.setValue(separator);
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
void Geant4_SoPolyhedron::clearAlternateRep(
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  alternateRep.setValue(NULL);
}
//////////////////////////////////////////////////////////////////////////////
void Geant4_SoPolyhedron::doAction(
 SoAction* aAction
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_ALTERNATEREP_DO_ACTION(aAction)
  SoShape::doAction(aAction);
}

#endif
