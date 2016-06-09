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
// $Id: G4OpenInventorSceneHandler.cc,v 1.39 2005/11/15 08:39:03 gbarrand Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor stored scene - creates OpenInventor display lists.
#ifdef G4VIS_BUILD_OI_DRIVER

// this :
#include "G4OpenInventorSceneHandler.hh"

#include <Inventor/SoPath.h>
#include <Inventor/SoNodeKitPath.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoCoordinate4.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoNormal.h>
#include <Inventor/nodes/SoNormalBinding.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/nodes/SoNurbsSurface.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoResetTransform.h>
#include <Inventor/nodes/SoMatrixTransform.h>

#define USE_SOPOLYHEDRON

#ifndef USE_SOPOLYHEDRON
#include "HEPVis/nodes/SoBox.h"
#include "HEPVis/nodes/SoTubs.h"
#include "HEPVis/nodes/SoCons.h"
#include "HEPVis/nodes/SoTrd.h"
#include "HEPVis/nodes/SoTrap.h"
#endif
#include "HEPVis/nodes/SoMarkerSet.h"
typedef HEPVis_SoMarkerSet SoMarkerSet;
#include "HEPVis/nodekits/SoDetectorTreeKit.h"
#include "HEPVis/misc/SoStyleCache.h"

#include "Geant4_SoPolyhedron.h"

#include "G4Scene.hh"
#include "G4NURBS.hh"
#include "G4OpenInventor.hh"
#include "G4OpenInventorTransform3D.hh"
#include "G4ThreeVector.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Transform3D.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polymarker.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ModelingParameters.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

G4Point3D translation;

G4int G4OpenInventorSceneHandler::fSceneIdCount = 0;

G4OpenInventorSceneHandler::G4OpenInventorSceneHandler (G4OpenInventor& system,
                                          const G4String& name)
:G4VSceneHandler (system, fSceneIdCount++, name)
,fRoot(0)
,fDetectorRoot(0)
,fTransientRoot(0)
,fCurrentSeparator(0)
,fModelingSolid(false)
,fReducedWireFrame(true)
,fStyleCache(0)
,fPreviewAndFull(false)
{
  fStyleCache = new SoStyleCache;
  fStyleCache->ref();

  fRoot = new SoSeparator;
  fRoot->ref();
  fRoot->setName("Root");
  
  fDetectorRoot = new SoSeparator;
  fDetectorRoot->setName("StaticRoot");
  fRoot->addChild(fDetectorRoot);
  
  fTransientRoot = new SoSeparator;
  fTransientRoot->setName("TransientRoot");
  fRoot->addChild(fTransientRoot);
  
  fCurrentSeparator = fTransientRoot;
}

G4OpenInventorSceneHandler::~G4OpenInventorSceneHandler ()
{
  fRoot->unref();
  fStyleCache->unref();
}

//
// Method for handling G4Polyline objects (from tracking).
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Polyline& line) {
  G4int nPoints = line.size();
  SbVec3f* pCoords = new SbVec3f[nPoints];

  for (G4int iPoint = 0; iPoint < nPoints ; iPoint++) {
    pCoords[iPoint].setValue((float)line[iPoint].x(),
                             (float)line[iPoint].y(),
                             (float)line[iPoint].z());
  }

  //
  // Color
  //
  const G4Colour& c = GetColour (line);
  SoMaterial* material = 
    fStyleCache->getMaterial((float)c.GetRed(),
                             (float)c.GetGreen(),
                             (float)c.GetBlue(),
                             (float)(1-c.GetAlpha()));
  fCurrentSeparator->addChild(material);
  
  //
  // Point Set
  // 
  SoCoordinate3 *polyCoords = new SoCoordinate3;
  polyCoords->point.setValues(0,nPoints,pCoords);
  fCurrentSeparator->addChild(polyCoords);
  
  //
  // Wireframe
  // 
  SoDrawStyle* drawStyle = fStyleCache->getLineStyle();
  fCurrentSeparator->addChild(drawStyle);

  SoLineSet *pLine = new SoLineSet;
#ifdef INVENTOR2_0
  pLine->numVertices.setValues(0,1,(const long *)&nPoints);
#else 
  pLine->numVertices.setValues(0,1,&nPoints);
#endif
  fCurrentSeparator->addChild(pLine);

  delete [] pCoords;
}

void G4OpenInventorSceneHandler::AddPrimitive (const G4Polymarker& polymarker) {
  G4int pointn = polymarker.size();
  if(pointn<=0) return;

  SbVec3f* points = new SbVec3f[pointn];
  for (G4int iPoint = 0; iPoint < pointn ; iPoint++) {
    points[iPoint].setValue((float)polymarker[iPoint].x(),
                            (float)polymarker[iPoint].y(),
                            (float)polymarker[iPoint].z());
  }

  const G4Colour& c = GetColour (polymarker);
  SoMaterial* material = 
    fStyleCache->getMaterial((float)c.GetRed(),
                             (float)c.GetGreen(),
                             (float)c.GetBlue(),
                             (float)(1-c.GetAlpha()));
  fCurrentSeparator->addChild(material);
  
  SoCoordinate3* coordinate3 = new SoCoordinate3;
  coordinate3->point.setValues(0,pointn,points);
  fCurrentSeparator->addChild(coordinate3);
  
  const G4VMarker& marker = polymarker;
  G4bool userSpecified = (marker.GetWorldSize() || marker.GetScreenSize());
  const G4VMarker& def = fpViewer->GetViewParameters().GetDefaultMarker();
  //const G4Vector3D& viewpointDirection =
  //  fpViewer->GetViewParameters().GetViewpointDirection();
  //const G4Vector3D& up = fpViewer->GetViewParameters().GetUpVector();
  G4double scale = fpViewer->GetViewParameters().GetGlobalMarkerScale();
  G4double size = scale *
    userSpecified ? marker.GetWorldSize() : def.GetWorldSize();
  G4double screenSize = 0;  // Calculate marker in screen size...
  if (size) {  // Size specified in world coordinates.
    screenSize = 10.;  // This needs to be much more sophisticated!
  } else { // Size specified in screen coordinates (pixels).
    screenSize = scale *
      userSpecified ? marker.GetScreenSize() : def.GetScreenSize();
  }

  SoMarkerSet* markerSet = new SoMarkerSet;
  markerSet->numPoints = pointn;

  G4VMarker::FillStyle style = marker.GetFillStyle();
  switch (polymarker.GetMarkerType()) {
  default:
    // Are available 5_5, 7_7 and 9_9
  case G4Polymarker::dots:
    if (screenSize <= 5.) {
      markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_5_5;
    } else if (screenSize <= 7.) {
      markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_7_7;
    } else {
      markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_9_9;
    }
    break;
  case G4Polymarker::circles:
    if (screenSize <= 5.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_5_5;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_5_5;
      }
    } else if (screenSize <= 7.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_7_7;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_7_7;
      }
    } else {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_9_9;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_9_9;
      }
    }
    break;
  case G4Polymarker::squares:
    if (screenSize <= 5.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_5_5;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_5_5;
      }
    } else if (screenSize <= 7.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_7_7;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_7_7;
      }
    } else {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_9_9;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_9_9;
      }
    }
  }
  fCurrentSeparator->addChild(markerSet);

  delete [] points;
}

// ********* NOTE ********* NOTE ********* NOTE ********* NOTE *********
//
//  This method (Text) has not been tested, as it is 
//  innaccessible from the menu in the current configuration
//
// ********* NOTE ********* NOTE ********* NOTE ********* NOTE *********
//
// Method for handling G4Text objects
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Text& text) {
  //
  // Color
  //
  const G4Colour& c = GetTextColour (text);
  SoMaterial* material = 
    fStyleCache->getMaterial((float)c.GetRed(),
                             (float)c.GetGreen(),
                             (float)c.GetBlue(),
                             (float)(1-c.GetAlpha()));
  fCurrentSeparator->addChild(material);

  //
  // Font
  // 

  SoFont *g4Font = new SoFont();
  g4Font->size = 20.0;
  fCurrentSeparator->addChild(g4Font);

  //
  // Text
  // 
  SoText2 *g4String = new SoText2();
  g4String->string.setValue("This is Text");
  g4String->spacing = 2.0;
  g4String->justification = SoText2::LEFT;
  fCurrentSeparator->addChild(g4String);
}

//
// Method for handling G4Circle objects
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Circle& circle) {
  AddCircleSquare(G4OICircle, circle);
}

//
// Method for handling G4Square objects - defaults to wireframe
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Square& square) {
  AddCircleSquare(G4OISquare, square);
}

void G4OpenInventorSceneHandler::AddCircleSquare
(G4OIMarker markerType, const G4VMarker& marker) {
  //
  // Color
  //
  const G4Colour& c = GetColour (marker);
  SoMaterial* material = 
    fStyleCache->getMaterial((float)c.GetRed(),
                             (float)c.GetGreen(),
                             (float)c.GetBlue(),
                             (float)(1-c.GetAlpha()));
  fCurrentSeparator->addChild(material);

  // A few useful quantities...
  G4Point3D centre = marker.GetPosition();
  G4bool userSpecified = (marker.GetWorldSize() || marker.GetScreenSize());
  const G4VMarker& def = fpViewer->GetViewParameters().GetDefaultMarker();
  //const G4Vector3D& viewpointDirection =
  //  fpViewer->GetViewParameters().GetViewpointDirection();
  //const G4Vector3D& up = fpViewer->GetViewParameters().GetUpVector();
  G4double scale = fpViewer->GetViewParameters().GetGlobalMarkerScale();
  G4double size = scale *
    userSpecified ? marker.GetWorldSize() : def.GetWorldSize();
  G4double screenSize = 0;  // Calculate marker in screen size...
  if (size) {  // Size specified in world coordinates.
    screenSize = 10.;  // This needs to be much more sophisticated!
  } else { // Size specified in screen coordinates (pixels).
    screenSize = scale *
      userSpecified ? marker.GetScreenSize() : def.GetScreenSize();
  }

  // Borrowed from AddPrimitive(G4Polymarker) - inefficient?
  SbVec3f* points = new SbVec3f[1];
  points[0].setValue((float)centre.x(),
		     (float)centre.y(),
		     (float)centre.z());
  SoCoordinate3* coordinate3 = new SoCoordinate3;
  coordinate3->point.setValues(0,1,points);
  fCurrentSeparator->addChild(coordinate3);
  SoMarkerSet* markerSet = new SoMarkerSet;
  markerSet->numPoints = 1;

  G4VMarker::FillStyle style = marker.GetFillStyle();
  switch (markerType) {
  case G4OICircle:
    if (screenSize <= 5.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_5_5;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_5_5;
      }
    } else if (screenSize <= 7.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_7_7;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_7_7;
      }
    } else {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_FILLED_9_9;
      } else {
	markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_9_9;
      }
    }
    break;
  case G4OISquare:
    if (screenSize <= 5.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_5_5;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_5_5;
      }
    } else if (screenSize <= 7.) {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_7_7;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_7_7;
      }
    } else {
      if (style == G4VMarker::filled) {
	markerSet->markerIndex = SoMarkerSet::SQUARE_FILLED_9_9;
      } else {
	markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_9_9;
      }
    }
  break;
  }
  fCurrentSeparator->addChild(markerSet);

  delete [] points;
}

//
// Method for handling G4Polyhedron objects for drawing solids.
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {
  if (polyhedron.GetNoFacets() == 0) return;
  Geant4_SoPolyhedron* soPolyhedron = new Geant4_SoPolyhedron(polyhedron);
  SbName sbName(fpCurrentLV?fpCurrentLV->GetName().c_str():"");
  soPolyhedron->setName(sbName);
  soPolyhedron->solid.setValue(fModelingSolid);
  soPolyhedron->reducedWireFrame.setValue(fReducedWireFrame?TRUE:FALSE);
  fCurrentSeparator->addChild(soPolyhedron);  
}

//
// Method for handling G4NURBS objects for drawing solids.
// Knots and Ctrl Pnts MUST be arrays of GLfloats.
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4NURBS& nurb) {

  G4float *u_knot_array, *u_knot_array_ptr;
  u_knot_array = u_knot_array_ptr = new G4float [nurb.GetnbrKnots(G4NURBS::U)];
  G4NURBS::KnotsIterator u_iterator (nurb, G4NURBS::U);
  while (u_iterator.pick (u_knot_array_ptr++));

  G4float *v_knot_array, *v_knot_array_ptr;
  v_knot_array = v_knot_array_ptr = new G4float [nurb.GetnbrKnots(G4NURBS::V)];
  G4NURBS::KnotsIterator v_iterator (nurb, G4NURBS::V);
  while (v_iterator.pick (v_knot_array_ptr++));

  G4float *ctrl_pnt_array, *ctrl_pnt_array_ptr;
  ctrl_pnt_array = ctrl_pnt_array_ptr =
    new G4float [nurb.GettotalnbrCtrlPts () * G4NURBS::NofC*sizeof(G4float)];
  G4NURBS::CtrlPtsCoordsIterator c_p_iterator (nurb);
  while (c_p_iterator.pick (ctrl_pnt_array_ptr++));
  
  SoSeparator *surfSep = new SoSeparator();

  //
  // Color
  //  
  const G4Colour& c = GetColour (nurb);
  SoMaterial* material = 
    fStyleCache->getMaterial((float)c.GetRed(),
                             (float)c.GetGreen(),
                             (float)c.GetBlue(),
                             (float)(1-c.GetAlpha()));
  surfSep->addChild(material);

  //
  // Set up NURBS
  //
  SoComplexity *complexity = new SoComplexity;
  SoCoordinate4 *ctlPts = new SoCoordinate4;
  SoNurbsSurface *oi_nurb = new SoNurbsSurface;
  
  complexity->value = (float)0.6;
  G4int    nPoints = nurb.GettotalnbrCtrlPts ();
  SbVec4f* points  = new SbVec4f[nPoints];
  for (G4int iPoint = 0; iPoint < nPoints ; iPoint++) {
    points[iPoint].setValue(
                            ctrl_pnt_array[iPoint*4 + 0],
                            ctrl_pnt_array[iPoint*4 + 1],
                            ctrl_pnt_array[iPoint*4 + 2],
                            ctrl_pnt_array[iPoint*4 + 3]);
  }
  ctlPts->point.setValues (0,nPoints,points);
  oi_nurb->numUControlPoints = nurb.GetnbrCtrlPts(G4NURBS::U);
  oi_nurb->numVControlPoints = nurb.GetnbrCtrlPts(G4NURBS::V);
  oi_nurb->uKnotVector.setValues(0,nurb.GetnbrKnots(G4NURBS::U),u_knot_array);
  oi_nurb->vKnotVector.setValues(0,nurb.GetnbrKnots(G4NURBS::V),v_knot_array);

  surfSep->addChild(complexity);
  surfSep->addChild(ctlPts);
  surfSep->addChild(oi_nurb);
  
  fCurrentSeparator->addChild(surfSep);

  //
  // Clean-up
  //
  delete [] u_knot_array;
  delete [] v_knot_array;
  delete [] ctrl_pnt_array;
  delete [] points;
}
//
// Generates prerequisites for primitives
//  
void G4OpenInventorSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation) {

  G4VSceneHandler::BeginPrimitives (objectTransformation);

  // Store away for future use, e.g.,
  // AddPrimitive (const G4Polyhedron& polyhedron).
  //  
  fpObjectTransformation = &objectTransformation;

  // The coming primitive is either:
  // a placed detector-type element whose destination
  // in the Scene Database has been predetermined for it.
  // In that case this routinde does absolutely nothing.

  // Or: an unplaced, transient, marker-type of object 
  // which needs to be properly placed, and whose 
  // destination (for now) is the root of the scene 
  // database.  For these types of objects, execute the
  // following code:
  //  
  if (fReadyForTransients) {

    // set the destination to "fTransientRoot"
    fCurrentSeparator = fTransientRoot;

    // place the transient object:
    fCurrentSeparator->addChild(fStyleCache->getResetTransform());

    SoMatrixTransform* matrixTransform = new SoMatrixTransform;
    G4OpenInventorTransform3D oiTran(objectTransformation);
    SbMatrix* sbMatrix = oiTran.GetSbMatrix();

    const G4Vector3D scale = fpViewer->GetViewParameters().GetScaleFactor();
    SbMatrix sbScale;
    sbScale.setScale
      (SbVec3f((float)scale.x(),(float)scale.y(),(float)scale.z()));
    sbMatrix->multRight(sbScale);

    matrixTransform->matrix.setValue(*sbMatrix);
    delete sbMatrix;
    fCurrentSeparator->addChild(matrixTransform);
  }
}

void G4OpenInventorSceneHandler::EndPrimitives () {
  G4VSceneHandler::EndPrimitives ();
}

void G4OpenInventorSceneHandler::EndModeling () {
}

void G4OpenInventorSceneHandler::ClearStore () {
  G4VSceneHandler::ClearStore();

  fDetectorRoot->removeAllChildren();
  fSeparatorMap.clear();

  fTransientRoot->removeAllChildren();
}

void G4OpenInventorSceneHandler::ClearTransientStore () {
  G4VSceneHandler::ClearTransientStore ();

  fTransientRoot->removeAllChildren();
}

void G4OpenInventorSceneHandler::PreAddSolid
(const G4Transform3D& objectTransformation,
 const G4VisAttributes& visAttribs) {

  // This assumes that it is a "Pre" to the "Add" of a
  // G4VPhysicalVolume.  This is true at present, but beware!  We
  // might have to think about the cast in the following code.

  G4VSceneHandler::PreAddSolid (objectTransformation, visAttribs);
  // Stores arguments away for future use, e.g., AddPrimitives.

  // This routines prepares to add solids to the scene database.  
  // We are expecting two
  // kinds of solids: leaf parts and non-leaf parts.  For non-leaf parts,
  // we create a detector tree kit.  This has two separators.  
  // The solid itself
  // goes in the preview separator, the full separator is 
  // forseen for daughters.
  // This separator is not only created--it is also put in a dictionary for 
  // future use by the leaf part.

  // For leaf parts, we first locate the mother volume and find its separator 
  // through the dictionary.

  // The private member fCurrentSeparator is set to the proper 
  // location on in the
  // scene database so that when the solid is actually 
  // added (in addthis), it is
  //  put in the right place.

  //G4PhysicalVolumeModel* pModel = fpModel->GetG4PhysicalVolumeModel();
  //G4bool thisToBeDrawn = pModel->IsThisCulled(fpCurrentLV,fpCurrentMaterial);
  //G4bool visible = !pModel->IsThisCulled(fpCurrentLV,0);

  // First find the color attributes.
  //
  const G4VisAttributes* pVisAttribs =
    fpViewer->GetApplicableVisAttributes (&visAttribs);
  const G4Colour& g4Col =  pVisAttribs->GetColour ();
  const double red = g4Col.GetRed ();
  const double green = g4Col.GetGreen ();
  const double blue = g4Col.GetBlue ();
  double transparency = 1 - g4Col.GetAlpha();

  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle(pVisAttribs);
  switch (drawing_style) {
  case (G4ViewParameters::wireframe):    
    fModelingSolid = false;
    break;
  case (G4ViewParameters::hlr):
  case (G4ViewParameters::hsr):
  case (G4ViewParameters::hlhsr):
    fModelingSolid = true;
    break;
  }	

  G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVisAttribs);
  fReducedWireFrame = !isAuxEdgeVisible;

  //printf("debug : PreAddSolid : %g %g %g : %d\n",
    //red,green,blue,pVisAttribs->IsVisible());
               
  if(!fpCurrentLV || !fpCurrentPV) return; //GB 

  //printf("debug : OIV : LV : %lx %s : %g %g %g\n",
  // fpCurrentLV,
  // fpCurrentLV->GetName().c_str(),
  // red,green,blue);

  if (fpCurrentLV->GetNoDaughters()!=0 ||
      fpCurrentPV->IsReplicated()) {
    // This block of code is executed for non-leaf parts:

    // Make the detector tree kit:
    SoDetectorTreeKit* detectorTreeKit = new SoDetectorTreeKit();  

    SoSeparator* previewSeparator   =  
      (SoSeparator*) detectorTreeKit->getPart("previewSeparator",TRUE);
    previewSeparator->renderCaching = SoSeparator::OFF;

    SoSeparator* fullSeparator =  
      (SoSeparator*) detectorTreeKit->getPart("fullSeparator",   TRUE);
    fullSeparator->renderCaching = SoSeparator::OFF;

    if(fPreviewAndFull) detectorTreeKit->setPreviewAndFull();
    else detectorTreeKit->setPreview(TRUE);

    SoMaterial* material = 
      fStyleCache->getMaterial((float)red,
                               (float)green,
                               (float)blue,
                               (float)transparency);
    detectorTreeKit->setPart("appearance.material",material);

    SoLightModel* lightModel = 
      fModelingSolid ? fStyleCache->getLightModelPhong() : 
                       fStyleCache->getLightModelBaseColor();
    detectorTreeKit->setPart("appearance.lightModel",lightModel);

    // Add the full separator to the dictionary; it is indexed by the 
    // address of the logical volume!
    // NOTE: the map is no longer built iteratively from the whole hierarchy
    //       of volumes since it is no longer possible to retrieve the mother
    //       physical volume. The algorithm requires review !   - GC
    //
    fSeparatorMap[fpCurrentPV->GetLogicalVolume()] = fullSeparator;

    // Find out where to add this volume.  This means locating its mother.  
    // If no mother can be found, it goes under root.
    G4LogicalVolume* MotherVolume = fpCurrentPV->GetMotherLogical();
    if (MotherVolume) {
      if (fSeparatorMap.find(MotherVolume) != fSeparatorMap.end()) {
        //printf("debug : PreAddSolid : mother found in map\n");
        fSeparatorMap[MotherVolume]->addChild(detectorTreeKit);
      } else {
        //printf("debug : PreAddSolid : mother not found in map !!!\n");
        fDetectorRoot->addChild(detectorTreeKit);
      }
    } else {
      //printf("debug : PreAddSolid : has no mother\n");
      fDetectorRoot->addChild(detectorTreeKit);
    }

    fCurrentSeparator = previewSeparator;
  } else {
    // This block of code is executed for leaf parts.

    // Locate the mother volume and find it's corresponding full separator
    G4LogicalVolume* MotherVolume = fpCurrentPV->GetMotherLogical();
    if (MotherVolume) {
      if (fSeparatorMap.find(MotherVolume) != fSeparatorMap.end()) {
        fCurrentSeparator = fSeparatorMap[MotherVolume];
      } else {
        fCurrentSeparator = fDetectorRoot;
      }
    } else {
      fCurrentSeparator = fDetectorRoot;
    }

    SoMaterial* material = 
      fStyleCache->getMaterial((float)red,
                               (float)green,
                               (float)blue,
                               (float)transparency);
    fCurrentSeparator->addChild(material);    

    SoLightModel* lightModel = 
      fModelingSolid ? fStyleCache->getLightModelPhong() : 
                       fStyleCache->getLightModelBaseColor();
    fCurrentSeparator->addChild(lightModel);
  }

  // Set up the geometrical transformation for the coming 
  fCurrentSeparator->addChild(fStyleCache->getResetTransform());

  SoMatrixTransform* matrixTransform = new SoMatrixTransform;
  G4OpenInventorTransform3D oiTran(objectTransformation);
  SbMatrix* sbMatrix = oiTran.GetSbMatrix();

  const G4Vector3D scale = fpViewer->GetViewParameters().GetScaleFactor();
  SbMatrix sbScale;
  sbScale.setScale
    (SbVec3f((float)scale.x(),(float)scale.y(),(float)scale.z()));
  sbMatrix->multRight(sbScale);

  matrixTransform->matrix.setValue(*sbMatrix);
  delete sbMatrix;
  fCurrentSeparator->addChild(matrixTransform);
}


//void G4OpenInventorSceneHandler::AddSolid(const G4VTrajectory& traj) {
//  G4VSceneHandler::AddSolid(traj);  // For now.
//}

//void G4OpenInventorSceneHandler::AddSolid(const G4VHit& hit) {
//  G4VSceneHandler::AddSolid(hit);  // For now.
//}

#endif
