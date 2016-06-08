// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GoSceneHandler.cc,v 1.3 1999/05/31 17:20:56 barrand Exp $
// GEANT4 tag $Name: geant4-00-01 $
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
  G4cout << "G4GoSceneHandler::G4GoSceneHandler" << endl;
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
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Polyline&) " << line.entries () << endl;
#endif
  SetColour        (GetColour(line));
  int              nPoints = line.entries ();
  OPointList       points;
  points           = OPointListCreate(nPoints);
  for (int iPoint = 0; iPoint < nPoints; iPoint++) {
    OPointListSetIthEntry   (points,iPoint,
			     line(iPoint).x(), 
			     line(iPoint).y(),
			     line(iPoint).z());
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
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Polymarker&) " << line.entries () << endl;
#endif
  SetColour  (GetColour(line));
  int        nPoints = line.entries ();
  OPointList points;
  points     = OPointListCreate(nPoints);
  for (int iPoint = 0; iPoint < nPoints; iPoint++) {
    OPointListSetIthEntry   (points,iPoint,
			     line(iPoint).x(), 
			     line(iPoint).y(),
			     line(iPoint).z());
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
  G4double  radius = circle.GetWorldSize();
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Circle&) : center : " 
    << " x : " << center.x() 
    << " y : " << center.y() 
    << " z : " << center.z() 
    << " ,radius : " << radius << endl;
#endif
  // Direction of circle ?
  //GoAddCircleToNode (fGoNode,radius,center.x(),center.y(),center.z(),0.,0.,1.);
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
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Polyhedron&) " << endl;
#endif
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

  G4bool        notLastFace;
  G4Normal3D    SurfaceNormal;
  double        xs[5];
  double        ys[5];
  double        zs[5];

  do 
    {
      G4bool        notLastEdge;
      G4Point3D     vertex;
      G4int         edgeFlag = 1;
      int           iPoint;
      iPoint        = 0;
      notLastFace   = polyhedron.GetNextNormal (SurfaceNormal);
      do 
	{
	  notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
	  xs[iPoint]  = vertex.x();
	  ys[iPoint]  = vertex.y();
	  zs[iPoint]  = vertex.z();
	  iPoint++;
	} while (notLastEdge);
      if( (iPoint!=3) && (iPoint!=4) )  
	{
	  CWarnF ("G4GoSceneHandler::AddPrimitive G4Polyhedron : not a quad or a triangle : %d\n",iPoint);
	}
      else
	{
	  if (modeling==OModelingWireFrame) 
	    {
	      xs[iPoint]  = xs[0];
	      ys[iPoint]  = ys[0];
	      zs[iPoint]  = zs[0];
	      iPoint++;
	      ONodeAddLines   (fGoNode,OContextGetStaticInstance(),iPoint,xs,ys,zs);
	    }
	  else if (modeling==OModelingSolid) 
	    {
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
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Text&) " << endl; 
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
  G4cout << "G4GoSceneHandler::AddPrimitive(G4Square&) " << endl; 
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
  G4cout << "G4GoSceneHandler::AddPrimitive(G4NURBS&) " << endl;
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
  G4VSceneHandler::AddThis (box);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Cons& cons
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4VSceneHandler::AddThis (cons);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Tubs& tubs
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4VSceneHandler::AddThis (tubs);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Trd& trd
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4VSceneHandler::AddThis (trd);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Trap& trap
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4VSceneHandler::AddThis (trap);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Sphere& sphere
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4VSceneHandler::AddThis (sphere);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Para& para
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4VSceneHandler::AddThis (para);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4Torus& torus
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4VSceneHandler::AddThis (torus);
}
/***************************************************************************/
void G4GoSceneHandler::AddThis (
 const G4VSolid& solid
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
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
  G4cout << "G4GoSceneHandler::PreAddThis : " << nodeName << endl;
#endif
}
/***************************************************************************/
void G4GoSceneHandler::PostAddThis (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4GoSceneHandler::PostAddThis" << endl;
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
  G4cout << "G4GoSceneHandler::BeginPrimitives" << endl;
#endif
  G4VSceneHandler::BeginPrimitives (objectTransformation);
  transformation  = objectTransformation;
  fGoNode         = nodeName!=NULL ? ONodeCreate (nodeName) : ONodeMake ();
  if(ONodeIsValid(fRootGoNode)==0) {
    fRootGoNode = ONodeMake ();
    fStaticRootGoNode = ONodeMake ();
    fTransientRootGoNode = ONodeMake ();
    ONodeAddChild   (fRootGoNode,fStaticRootGoNode);
    ONodeAddChild   (fRootGoNode,fTransientRootGoNode);
  }
  if (fReadyForTransients) {
    ONodeAddChild   (fTransientRootGoNode,fGoNode);
  } else {
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
  G4cout << "G4GoSceneHandler::EndPrimitives" << endl;
#endif
  G4double        g4angle;
  G4ThreeVector   g4axis;
  G4ThreeVector   g4trans;
  transformation.getRotation().getAngleAxis(g4angle,g4axis);
  g4trans         = transformation.getTranslation();
  OMatrix         orot = OMatrixCreate (OMatrixRotationAxis,g4angle,g4axis.x(),g4axis.y(),g4axis.z());
  OMatrix         otra = OMatrixCreate (OMatrixTranslation,g4trans.x(),g4trans.y(),g4trans.z());
  OMatrix         omatrix = OMatrixMultiply (otra,orot);
  ONodeSetMatrix  (fGoNode,omatrix);
  OMatrixDelete   (omatrix);
  OMatrixDelete   (orot);
  OMatrixDelete   (otra);
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
  G4cout << "G4GoSceneHandler::ClearStore" << endl;
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
  G4cout << "G4GoSceneHandler::ClearTransientStore" << endl;
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
  G4cout << "G4GoSceneHandler::RequestPrimitives" << endl;
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
  G4cout << "G4GoSceneHandler::RequestPrimitives" << endl;
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
                              " index : " << index                << endl;
#endif
  OContextSetColorIndex (OContextGetStaticInstance(),index);
}

#endif
