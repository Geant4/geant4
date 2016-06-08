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
// $Id: G4oScene.cc,v 1.7.4.1 2001/06/28 19:06:34 gunter Exp $
// GEANT4 tag $Name:  $
//
/* +---------------------- Copyright notice -------------------------------+ */
/* | Copyright (C) 1995, Guy Barrand, LAL Orsay, (barrand@lal.in2p3.fr)    | */
/* |   Permission to use, copy, modify, and distribute this software       | */
/* |   and its documentation for any purpose and without fee is hereby     | */
/* |   granted, provided that the above copyright notice appear in all     | */
/* |   copies and that both that copyright notice and this permission      | */
/* |   notice appear in supporting documentation.  This software is        | */
/* |   provided "as is" without express or implied warranty.               | */
/* +---------------------- Copyright notice -------------------------------+ */
/*
#define DEBUG
*/

#ifdef DEBUG
#include <stdio.h>
#endif

/*Co*/
#include <CPrinter.h>
#include <CString.h>
#include <OMatrix.h>

/*Go*/
#include <Go.h>
#include <OCamera.h>

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Cons.hh>
#include <G4Trd.hh>
#include <G4Trap.hh>
#include <G4Sphere.hh>
#include <G4Para.hh>
#include <G4Torus.hh>
#include <G4Polycone.hh>
#include <G4Polyhedra.hh>
#include <G4Polyhedron.hh>
#include <G4Polyline.hh>
#include <G4Transform3D.hh>
#include <G4PhysicalVolumeModel.hh>
#include <G4Normal3D.hh>

#include <G4oScene.hh>

ONode G4oScene::node     = NULL;
char* G4oScene::nodeName = NULL;
/***************************************************************************/
G4oScene::G4oScene (
)
:fpAccTransf(NULL)
,fCurrentDepth(0)
,fpCurrentPV(NULL)
,fpCurrentLV(NULL)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
}
/***************************************************************************/
G4oScene::~G4oScene (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
  CStringDelete (nodeName);
  nodeName      = NULL;
}
/***************************************************************************/
void G4oScene::PreAddThis (
 const G4Transform3D&  objectTransformation
,const G4VisAttributes& /*visAttribs*/
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  fpAccTransf = &objectTransformation;
}
/***************************************************************************/
void G4oScene::PostAddThis (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4oScene::EstablishSpecials (
 G4PhysicalVolumeModel& a_model
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  a_model.DefinePointersToWorkingSpace (&fCurrentDepth,
					&fpCurrentPV,
					&fpCurrentLV);
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Box& box
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : box.\n");
#endif
  AddPolyhedron (box.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Cons& cons
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : cons.\n");
#endif
  AddPolyhedron (cons.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Tubs& tubs
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : tubs.\n");
#endif
  AddPolyhedron (tubs.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Trd& trd
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : trd.\n");
#endif
  AddPolyhedron (trd.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Trap& trap
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : trap.\n");
#endif
  AddPolyhedron (trap.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Sphere& sphere
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : sphere.\n");
#endif
  AddPolyhedron (sphere.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Para& para
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : para.\n");
#endif
  AddPolyhedron (para.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Torus& torus
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : torus.\n");
#endif
  AddPolyhedron (torus.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Polycone& polycone
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : polycone.\n");
#endif
  AddPolyhedron (polycone.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4Polyhedra& polyhedra
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : polyhedra.\n");
#endif
  AddPolyhedron (polyhedra.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::AddThis (
 const G4VSolid& solid
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : solid.\n");
#endif
  AddPolyhedron (solid.CreatePolyhedron());
}
/***************************************************************************/
void G4oScene::BeginPrimitives(
 const G4Transform3D& objectTransformation
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  //G4VSceneHandler::BeginPrimitives (objectTransformation);
}
/***************************************************************************/
void G4oScene::EndPrimitives(
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  //G4VSceneHandler::EndPrimitives ();
}
/***************************************************************************/
void G4oScene::AddPrimitive(
 const G4Polyline&
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  printf ("debug : AddPrimitive : polyline.\n");
#endif
}
/***************************************************************************/
void G4oScene::AddPrimitive(
 const G4Text&
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  printf ("debug : AddPrimitive : text.\n");
#endif
}
/***************************************************************************/
void G4oScene::AddPrimitive(
 const G4Circle&
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  printf ("debug : AddPrimitive : circle.\n");
#endif
}
/***************************************************************************/
void G4oScene::AddPrimitive(
 const G4Square&
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  printf ("debug : AddPrimitive : square.\n");
#endif
}
/***************************************************************************/
void G4oScene::AddPrimitive(
 const G4Polymarker&
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  printf ("debug : AddPrimitive : polymarker.\n");
#endif
}
/***************************************************************************/
void G4oScene::AddPrimitive(
 const G4Polyhedron&
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  printf ("debug : AddPrimitive : polyhedron.\n");
#endif
}
/***************************************************************************/
void G4oScene::AddPrimitive(
 const G4NURBS&
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  printf ("debug : AddPrimitive : NURBS.\n");
#endif
}
/***************************************************************************/
void G4oScene::SetNodeName (
 char* a_name
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
  CStringDelete (nodeName);
  nodeName      = CStringDuplicate (a_name);
}
/***************************************************************************/
ONode G4oScene::GetNode (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
  return node;
}
/***************************************************************************/
void G4oScene::AddPolyhedron (
 G4Polyhedron* a_polyhedron
)
/***************************************************************************/
/*
  Loop to get vertices from G4OpenGLScene::AddPrimitive.
*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4bool        notLastFace;
  G4Normal3D    SurfaceUnitNormal;
  double        xs[5];
  double        ys[5];
  double        zs[5];
/*.........................................................................*/
#ifdef DEBUG
  printf ("debug : G4oScene::AddPolyhedron : ONode : %s\n",nodeName);
#endif
  node = ONodeCreate (nodeName);
  do 
    {
      G4int         notLastEdge;
      G4Point3D     vertex;
      G4int         edgeFlag (true);
      int           iPoint;
      iPoint        = 0;
      notLastFace   = a_polyhedron->GetNextNormal (SurfaceUnitNormal);
      do 
	{
	  notLastEdge = a_polyhedron->GetNextVertex (vertex, edgeFlag);
	  xs[iPoint]  = vertex.x();
	  ys[iPoint]  = vertex.y();
	  zs[iPoint]  = vertex.z();
	  iPoint++;
	} while (notLastEdge);
      if( (iPoint!=3) && (iPoint!=4) )  
	{
	  CWarnF ("G4oScene::AddPolyhedron : not a quad or a triangle : %d\n",iPoint);
	}
      else
	{
	  if (OContextGetModeling(OContextGetStaticInstance())==OModelingWireFrame) 
	    {
	      xs[iPoint]  = xs[0];
	      ys[iPoint]  = ys[0];
	      zs[iPoint]  = zs[0];
	      iPoint++;
	      ONodeAddLines   (node,OContextGetStaticInstance(),iPoint,xs,ys,zs);
	    }
	  else if (OContextGetModeling(OContextGetStaticInstance())==OModelingSolid) 
	    {
	      ONodeAddPolygon (node,OContextGetStaticInstance(),iPoint,xs,ys,zs);
	    }
#ifdef DEBUG
  printf ("debug : G4oScene::AddPrimitive : AddLines : %d\n",iPoint);
#endif
	}
    } while (notLastFace);

  if(fpAccTransf==NULL) {
    CWarnF ("G4oScene::AddPolyhedron : no transformation given !\n");
  } else {
    G4double        g4angle;
    G4ThreeVector   g4axis;
    G4ThreeVector   g4trans;
    fpAccTransf->getRotation().getAngleAxis(g4angle,g4axis);
    g4trans         = fpAccTransf->getTranslation();
    OMatrix         orot = OMatrixCreate (OMatrixRotationAxis,g4angle,g4axis.x(),g4axis.y(),g4axis.z());
    OMatrix         otra = OMatrixCreate (OMatrixTranslation,g4trans.x(),g4trans.y(),g4trans.z());
    OMatrix         omatrix = OMatrixMultiply (otra,orot);
    ONodeSetMatrix  (node,omatrix);
    OMatrixDelete   (omatrix);
    OMatrixDelete   (orot);
    OMatrixDelete   (otra);
  }
}
