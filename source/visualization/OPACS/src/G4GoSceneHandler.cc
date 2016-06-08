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
// $Id: G4GoSceneHandler.cc,v 1.13 2001/11/12 18:22:09 johna Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// 
// Guy Barrand 04 November 1996
// Wo stored scene - creates Wo display lists.

//#define DEBUG

#ifdef G4VIS_BUILD_OPACS_DRIVER

//Co
#include <CPrinter.h>
#include <CString.h>
//Go
#include <Go.h>
//G4
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Polyhedron.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Text.hh"
#include "G4NURBS.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
//This
#include "G4GoSceneHandler.hh"

static G4Transform3D transformation;

G4int     G4GoSceneHandler::fSceneIdCount = 0;
G4int     G4GoSceneHandler::fSceneCount   = 0;
ONode     G4GoSceneHandler::fGoNode       = NULL;
OColormap G4GoSceneHandler::fOColormap    = NULL;
/***************************************************************************/
G4GoSceneHandler::G4GoSceneHandler (
 G4VGraphicsSystem& system,
 const G4String& name
)
/***************************************************************************/
// Creator.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
:G4VSceneHandler    (system, fSceneIdCount++, name)
,fSystem(system)
,fRootGoNode(NULL)
,fStaticRootGoNode(NULL)
,fTransientRootGoNode(NULL)
,nodeName(NULL)
/*.........................................................................*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::G4GoSceneHandler" << G4endl;
#endif
  fSceneCount++;
}
/***************************************************************************/
G4GoSceneHandler::~G4GoSceneHandler (
)
/***************************************************************************/
// Destructor.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  //fRootGoNode could be deleted by an XoCamera popup "Erase".
  if(ONodeIsValid(fRootGoNode)==1) 
    ONodeDelete (fRootGoNode); 
  fRootGoNode = NULL;
  fStaticRootGoNode = NULL;
  fTransientRootGoNode = NULL;
  fSceneCount--;
}
/***************************************************************************/
void G4GoSceneHandler::AddPrimitive (
 const G4Polyline& line
) 
/***************************************************************************/
// Method for handling G4Polyline objects.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Polyline&) " << line.size () << G4endl;
#endif
  SetColour        (GetColour(line));
  int              nPoints = line.size ();
  OPointList       points;
  points           = OPointListCreate(nPoints);
  for (int iPoint = 0; iPoint < nPoints; iPoint++) {
    OPointListSetIthEntry   (points,iPoint,
			     line[iPoint].x(), 
			     line[iPoint].y(),
			     line[iPoint].z());
  }
  GoAddLinesToNode (fGoNode,points);
  OPointListDelete (points);
}
/***************************************************************************/
void G4GoSceneHandler::AddPrimitive (
 const G4Polymarker& line
) 
/***************************************************************************/
// Method for handling G4Polymarker objects.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Polymarker&) " << line.size () << G4endl;
#endif
  SetColour  (GetColour(line));
  int        nPoints = line.size ();
  OPointList points;
  points     = OPointListCreate(nPoints);
  for (int iPoint = 0; iPoint < nPoints; iPoint++) {
    OPointListSetIthEntry   (points,iPoint,
			     line[iPoint].x(), 
			     line[iPoint].y(),
			     line[iPoint].z());
  }
  GoAddMarkersToNode (fGoNode,points);
  OPointListDelete   (points);
}
/***************************************************************************/
void G4GoSceneHandler::AddPrimitive (
 const G4Circle& circle
) 
/***************************************************************************/
// Method for handling G4Circle objects.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  SetColour (GetColour(circle));
  G4Point3D center = circle.GetPosition();
  // G4double  radius = circle.GetWorldSize();  ####### Unused. #######
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Circle&) : center : " 
    << " x : " << center.x() 
    << " y : " << center.y() 
    << " z : " << center.z() 
    << " ,radius : " << radius << G4endl;
#endif
  // Direction of circle ?
  //GoAddCircleToNode (fGoNode,radius,center.x(),center.y(),center.z(),0.,0.,1.);
  OContextSetMarkerStyle (OContextGetStaticInstance(),OMarkerStyleCircle);
  OContextSetMarkerSize (OContextGetStaticInstance(),5);
  GoAddMarkerToNode (fGoNode,center.x(),center.y(),center.z());
}
/***************************************************************************/
void G4GoSceneHandler::AddPrimitive (
 const G4Polyhedron& polyhedron
) 
/***************************************************************************/
// Method for handling G4Polyhedron objects for drawing solids.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Polyhedron&) " << G4endl;
#endif

  if (polyhedron.GetNoFacets() == 0) return;

  SetColour     (GetColour(polyhedron));

  OModeling     modeling = OModelingWireFrame;
  switch (GetDrawingStyle (polyhedron)) {
  case (G4ViewParameters::hlhsr):
  case (G4ViewParameters::hsr):
    //not ready yet modeling = OModelingSolid;
    break;
  case (G4ViewParameters::hlr):
  case (G4ViewParameters::wireframe):
  default:
    modeling = OModelingWireFrame;
    break;
  }	

  double xs[5];
  double ys[5];
  double zs[5];
  G4Normal3D SurfaceNormal;
  G4Point3D vertex;

  G4bool notLastFace;
  do {
    G4int edgeFlag = 1;
    G4int lastEdgeFlag = edgeFlag;
    int   iPoint  = 0;
    notLastFace = polyhedron.GetNextNormal (SurfaceNormal);
    G4bool notLastEdge;
    do {
      notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
      if (modeling==OModelingWireFrame) {
	if (edgeFlag != lastEdgeFlag) {
	  if (edgeFlag) {
	    // Pass to visible edges, restart to accumulate :
	    iPoint = 0;
	    xs[iPoint]  = vertex.x();
	    ys[iPoint]  = vertex.y();
	    zs[iPoint]  = vertex.z();
	    iPoint++;
	  } else {
	    // Pass to invisible edges, flush the visible ones :
	    xs[iPoint]  = vertex.x();
	    ys[iPoint]  = vertex.y();
	    zs[iPoint]  = vertex.z();
	    iPoint++;
	    if(iPoint>=2) {
	      ONodeAddLines(fGoNode,OContextGetStaticInstance(),
			    iPoint,xs,ys,zs);
	    }
	  }
	  lastEdgeFlag = edgeFlag;
	} else {
	  if (edgeFlag) {
	    xs[iPoint]  = vertex.x();
	    ys[iPoint]  = vertex.y();
	    zs[iPoint]  = vertex.z();
	    iPoint++;
	  } else {
	  }
	}
      } else {
	xs[iPoint]  = vertex.x();
	ys[iPoint]  = vertex.y();
	zs[iPoint]  = vertex.z();
	iPoint++;
      }
    } while (notLastEdge);
    if (modeling==OModelingWireFrame) {
      if(iPoint!=0) { // Some line started, flush it :
	/*
	xs[iPoint]  = xs[0];
	ys[iPoint]  = ys[0];
	zs[iPoint]  = zs[0];
	iPoint++;
	*/
	if(iPoint>=2) {
	  ONodeAddLines(fGoNode,OContextGetStaticInstance(),iPoint,xs,ys,zs);
	}
      }
    } else if (modeling==OModelingSolid) {
      if( (iPoint!=3) && (iPoint!=4) ) {
	CWarnF ("G4GoSceneHandler::AddPrimitive G4Polyhedron : not a quad or a triangle : %d\n",iPoint);
      } else {
	ONodeAddPolygon (fGoNode,OContextGetStaticInstance(),iPoint,xs,ys,zs);
      }
    }
  } while (notLastFace);
}
/***************************************************************************/
void G4GoSceneHandler::AddPrimitive (
 const G4Text& text
) 
/***************************************************************************/
// Method for handling G4Text objects.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  SetColour (GetTextColour(text));
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Text&) " << G4endl; 
#endif
  CWarnF ("G4GoSceneHandler::AddPrimitive G4Text : not yet implemented.\n");
}
/***************************************************************************/
void G4GoSceneHandler::AddPrimitive (
 const G4Square& Square
) 
/***************************************************************************/
// Method for handling G4Square objects.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  SetColour (GetColour(Square));
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Square&) " << G4endl; 
#endif
  CWarnF ("G4GoSceneHandler::AddPrimitive G4Square : not yet implemented.\n");
}
/***************************************************************************/
void G4GoSceneHandler::AddPrimitive (
 const G4NURBS& nurb
) 
/***************************************************************************/
//Method for handling G4NURBS objects for drawing solids.
//Knots and Ctrl Pnts MUST be arrays of GLfloats.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  SetColour (GetColour(nurb));
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddPrimitive(G4NURBS&) " << G4endl;
#endif
  CWarnF ("G4GoSceneHandler::AddPrimitive G4NURBS : not yet implemented.\n");
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Box& box
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Box&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (box);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Cons& cons
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Cons&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (cons);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Tubs& tubs
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Tubs&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (tubs);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Trd& trd
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Trd&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (trd);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Trap& trap
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Trap&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (trap);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Sphere& sphere
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Sphere&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (sphere);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Para& para
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Para&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (para);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Torus& torus
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Torus&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (torus);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Polycone& polycone
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Polycone&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (polycone);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Polyhedra& polyhedra
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4Polyhedra&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (polyhedra);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4VSolid& solid
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddThis(G4VSolid&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (solid);
}
/***************************************************************************/
void G4GoSceneHandler::PreAddThis (
			    const G4Transform3D& objectTransformation,
			    const G4VisAttributes& visAttribs
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4VSceneHandler::PreAddThis (objectTransformation, visAttribs);

  CStringDelete (nodeName);
  nodeName      = NULL;

  const G4VPhysicalVolume* currentPhysicalVolume = fpCurrentPV;

  if(currentPhysicalVolume==NULL) return;

  nodeName = CStringCreateF (15+64,"PhysicalVolume/%lu",currentPhysicalVolume );

#ifdef DEBUG
  G4cout << "G4GoSceneHandler::PreAddThis : " << nodeName << G4endl;
#endif
}
/***************************************************************************/
void G4GoSceneHandler::PostAddThis (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::PostAddThis" << G4endl;
#endif
  CStringDelete (nodeName);
  nodeName      = NULL;
  G4VSceneHandler::PostAddThis();
}
/***************************************************************************/
void G4GoSceneHandler::BeginPrimitives (
				 const G4Transform3D& objectTransformation
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::BeginPrimitives" << G4endl;
#endif

  G4VSceneHandler::BeginPrimitives (objectTransformation);
  transformation  = objectTransformation;

  fGoNode         = nodeName!=NULL ? ONodeCreate (nodeName) : ONodeMake ();
  if(ONodeIsValid(fRootGoNode)==0) {
#ifdef DEBUG
    G4cout << " build Root, Static, Transient nodes" << G4endl;
#endif
    fRootGoNode = ONodeMake ();
    fStaticRootGoNode = ONodeMake ();
    fTransientRootGoNode = ONodeMake ();
    ONodeAddChild   (fRootGoNode,fStaticRootGoNode);
    ONodeAddChild   (fRootGoNode,fTransientRootGoNode);
  } 
  if (fReadyForTransients) {
#ifdef DEBUG
    G4cout << " add node to transient scene graph" << G4endl;
#endif
    ONodeAddChild   (fTransientRootGoNode,fGoNode);
  } else {
#ifdef DEBUG
    G4cout << " add node to static scene graph" << G4endl;
#endif
    ONodeAddChild   (fStaticRootGoNode,fGoNode);
  }
}
/***************************************************************************/
void G4GoSceneHandler::EndPrimitives (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::EndPrimitives" << G4endl;
#endif
  HepRotation     rot = transformation.getRotation();
  Hep3Vector      tra = transformation.getTranslation();
  OMatrix         omatrix = OMatrixCreate(OMatrixFollowing,
		  rot.xx(),rot.xy(),rot.xz(),tra.x(),
		  rot.yx(),rot.yy(),rot.yz(),tra.y(),
		  rot.zx(),rot.zy(),rot.zz(),tra.z(),
			0.,      0.,      0.,     1.); 
  ONodeSetMatrix  (fGoNode,omatrix);
  OMatrixDelete   (omatrix);
  fGoNode         = NULL;
  G4VSceneHandler::EndPrimitives ();
}
/***************************************************************************/
void G4GoSceneHandler::ClearStore (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::ClearStore" << G4endl;
#endif
  if(ONodeIsValid(fTransientRootGoNode)==0) fTransientRootGoNode = NULL;
  ONodeDestroyChildren (fTransientRootGoNode);
  if(ONodeIsValid(fStaticRootGoNode)==0) fStaticRootGoNode = NULL;
  ONodeDestroyChildren (fStaticRootGoNode);
}
/***************************************************************************/
void G4GoSceneHandler::ClearTransientStore (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::ClearTransientStore" << G4endl;
#endif
  if(ONodeIsValid(fTransientRootGoNode)==0) fTransientRootGoNode = NULL;
  ONodeDestroyChildren (fTransientRootGoNode);
}
/***************************************************************************/
void G4GoSceneHandler::RequestPrimitives (
 const G4VSolid& solid
) 
/***************************************************************************/
// Stop-gap solution for display List re-use.
// A proper implementation would use geometry hierarchy.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::RequestPrimitives" << G4endl;
#endif
  G4VSceneHandler::RequestPrimitives (solid);
}
/***************************************************************************/
ONode G4GoSceneHandler::GetRootNode (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::RequestPrimitives" << G4endl;
#endif
  if(ONodeIsValid(fRootGoNode)==0) {
    fRootGoNode = NULL;
    fStaticRootGoNode = NULL;
    fTransientRootGoNode = NULL;
  }
  return fRootGoNode;
}
/***************************************************************************/
void G4GoSceneHandler::SetColour (
 const G4Colour& a_colour
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(fOColormap==NULL) fOColormap = OColormapGetIdentifier("ocolormap_X");
  int index = OColormapGetRGB_Index(fOColormap,a_colour.GetRed (), a_colour.GetGreen (), a_colour.GetBlue ());
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::SetColour : red : " << a_colour.GetRed ()   <<
                              " green : " << a_colour.GetGreen () <<  
                               " blue : " << a_colour.GetBlue ()  <<
                              " index : " << index                << G4endl;
#endif
  OContextSetColorIndex (OContextGetStaticInstance(),index);
}

#endif
