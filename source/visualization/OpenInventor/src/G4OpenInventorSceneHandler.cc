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
// $Id: G4OpenInventorSceneHandler.cc,v 1.20 2004-11-11 16:14:51 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoResetTransform.h>

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

typedef SoDetectorTreeKit SoDetectorTreeKit;
G4Point3D translation;

G4int G4OpenInventorSceneHandler::fSceneIdCount = 0;

G4int G4OpenInventorSceneHandler::fSceneCount = 0;

G4OpenInventorSceneHandler::G4OpenInventorSceneHandler (G4OpenInventor& system,
                                          const G4String& name)
:G4VSceneHandler (system, fSceneIdCount++, name)
,fRoot(0)
,fDetectorRoot(0)
,fTransientRoot(0)
,fCurrentSeparator(0)
,fModelingSolid(false)
{
  fSceneCount++;
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
  fSceneCount--;
}

//
// Method for handling G4Polyline objects (from tracking).
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Polyline& line) {
  G4int nPoints = line.size();
  SbVec3f* pCoords = new SbVec3f[nPoints];

  for (G4int iPoint = 0; iPoint < nPoints ; iPoint++) {
    pCoords[iPoint].setValue(line[iPoint].x(),
                             line[iPoint].y(),
                             line[iPoint].z());
  }

  //
  // Color
  //
  SoMaterial *pMaterial = new SoMaterial;
  const G4Colour& c = GetColour (line);
  pMaterial->diffuseColor.setValue(c.GetRed (), c.GetGreen (), c.GetBlue ());
  fCurrentSeparator->addChild(pMaterial);
  
  //
  // Point Set
  // 
  SoCoordinate3 *polyCoords = new SoCoordinate3;
  polyCoords->point.setValues(0,nPoints,pCoords);
  fCurrentSeparator->addChild(polyCoords);
  
  //
  // Wireframe
  // 
  SoDrawStyle *drawStyle = new SoDrawStyle;
  drawStyle->style = SoDrawStyle::LINES;
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
  //G4VSceneHandler::AddPrimitive (polymarker);

  G4int pointn = polymarker.size();
  if(pointn<=0) return;

  if(polymarker.GetMarkerType()==G4Polymarker::line) {
    // Why a polymarker of type line is not a G4Polyline ???
    G4Polyline polyline (polymarker);
    for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
      polyline.push_back (polymarker[iPoint]);
    }
    AddPrimitive(polyline);
    return;
  }

  SbVec3f* points = new SbVec3f[pointn];
  for (G4int iPoint = 0; iPoint < pointn ; iPoint++) {
    points[iPoint].setValue(polymarker[iPoint].x(),
                            polymarker[iPoint].y(),
                            polymarker[iPoint].z());
  }

  SoMaterial* material = new SoMaterial;
  const G4Colour& c = GetColour (polymarker);
  material->diffuseColor.setValue(c.GetRed (), c.GetGreen (), c.GetBlue ());
  fCurrentSeparator->addChild(material);
  
  SoCoordinate3* coordinate3 = new SoCoordinate3;
  coordinate3->point.setValues(0,pointn,points);
  fCurrentSeparator->addChild(coordinate3);
  
  SoMarkerSet* markerSet = new SoMarkerSet;
  markerSet->numPoints = pointn;
  switch (polymarker.GetMarkerType()) {
  default:
  // Are available 5_5, 7_7 and 9_9
  case G4Polymarker::dots:{
    markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_5_5;
  }break;
  case G4Polymarker::circles:{
    markerSet->markerIndex = SoMarkerSet::CIRCLE_LINE_7_7;
  }break;
  case G4Polymarker::squares:{
    markerSet->markerIndex = SoMarkerSet::SQUARE_LINE_7_7;
  }break;
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
  SoMaterial *pMaterial = new SoMaterial;
  const G4Colour& c = GetTextColour (text);
  pMaterial->diffuseColor.setValue(c.GetRed (), c.GetGreen (), c.GetBlue ());
  //pMaterial->diffuseColor.setValue(1.0, 0.0, 1.0); // Magenta so we can see
  fCurrentSeparator->addChild(pMaterial);

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
  //
  // Color
  //
  SoMaterial *pMaterial = new SoMaterial;
  const G4Colour& c = GetColour (circle);
  pMaterial->diffuseColor.setValue(c.GetRed (), c.GetGreen (), c.GetBlue ());
  fCurrentSeparator->addChild(pMaterial);

  //
  // Dimensions
  //
  // (userSpecified and defaultMarker are unused - what was intended here?)
  //G4double userSpecified = circle.GetWorldSize() || circle.GetScreenSize();
  //const G4VMarker& defaultMarker =
  //   fpViewer -> GetViewParameters().GetDefaultMarker();

  G4double size = GetMarkerSize ( circle );

  //
  // Position
  //
  
  const G4Point3D cCtr = circle.GetPosition();
   
  SoTranslation *cTrans = new SoTranslation;
  cTrans->translation.setValue(cCtr.x(),cCtr.y(),cCtr.z());
  fCurrentSeparator->addChild(cTrans);

  //
  // Sphere
  // 

  SoSphere *g4Sphere = new SoSphere();
  g4Sphere->radius = size;
  fCurrentSeparator->addChild(g4Sphere);
}

//
// Method for handling G4Square objects - defaults to wireframe
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Square& Square) {
  //
  // Color
  //
  SoMaterial *pMaterial = new SoMaterial;
  const G4Colour& c = GetColour (Square);
  pMaterial->diffuseColor.setValue(c.GetRed (), c.GetGreen (), c.GetBlue ());
  fCurrentSeparator->addChild(pMaterial);

  //
  // Size
  //
  G4double sSize = GetMarkerSize(Square);
  sSize = sSize * 2.0;

  //
  // Position
  //
  const G4Point3D sOrig = Square.GetPosition();
   
  SoTranslation *sTrans = new SoTranslation;
  sTrans->translation.setValue(sOrig.x(),sOrig.y(),sOrig.z());
  fCurrentSeparator->addChild(sTrans);
  
  SoCube *g4Square = new SoCube();
  g4Square->width =  sSize;
  g4Square->height = sSize;
  g4Square->depth =  sSize;
  fCurrentSeparator->addChild(g4Square);
}

//
// Method for handling G4Polyhedron objects for drawing solids.
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {
  if (polyhedron.GetNoFacets() == 0) return;

  //
  // Define Material
  //
  SoMaterial *pMaterial = new SoMaterial;
  const G4Colour& c = GetColour (polyhedron);
  pMaterial->diffuseColor.setValue(c.GetRed (), c.GetGreen (), c.GetBlue ());
  fCurrentSeparator->addChild(pMaterial);

  G4bool modelingSolid = false;
  SoLightModel* lightModel = new SoLightModel;
  lightModel->model.setValue(modelingSolid ?
			     SoLightModel::PHONG:
			     SoLightModel::BASE_COLOR);
  fCurrentSeparator->addChild(lightModel);

  Geant4_SoPolyhedron* soPolyhedron = new Geant4_SoPolyhedron(polyhedron);
  soPolyhedron->solid.setValue(modelingSolid);
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
  SoMaterial *pMaterial = new SoMaterial;
  const G4Colour& c = GetColour (nurb);
  pMaterial->diffuseColor.setValue(c.GetRed (), c.GetGreen (), c.GetBlue ());
  surfSep->addChild(pMaterial);

  //
  // Set up NURBS
  //
  SoComplexity *complexity = new SoComplexity;
  SoCoordinate4 *ctlPts = new SoCoordinate4;
  SoNurbsSurface *oi_nurb = new SoNurbsSurface;
  
  complexity->value = 0.6;
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
    //  
    // set the destination to "fTransientRoot"
    //  
    fCurrentSeparator = fTransientRoot;
    //  
    // place the transient object:
    //  
    G4OpenInventorTransform3D oiTran (objectTransformation);
    SoSFMatrix *oiMat = oiTran.GetOIMatrix();
    SoMatrixTransform *xform = new SoMatrixTransform;
    xform->matrix.setValue(oiMat->getValue());
    //  
    // add a transform.
    //  
    fCurrentSeparator->addChild(new SoResetTransform);
    fCurrentSeparator->addChild(xform);
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
  fTransientRoot->removeAllChildren();
}

void G4OpenInventorSceneHandler::ClearTransientStore () {
  fTransientRoot->removeAllChildren();
}

void G4OpenInventorSceneHandler::RequestPrimitives (const G4VSolid& solid) {
    // Stop-gap solution for display List re-use.
    // A proper implementation would use geometry hierarchy.
    //
    G4VSceneHandler::RequestPrimitives (solid);
}


void G4OpenInventorSceneHandler::AddThis (const G4Box & box) {
#ifdef USE_SOPOLYHEDRON
  Geant4_SoPolyhedron* soPolyhedron = 
    new Geant4_SoPolyhedron(box.CreatePolyhedron());
  soPolyhedron->solid.setValue(fModelingSolid?TRUE:FALSE);
  fCurrentSeparator->addChild(soPolyhedron);  
#else
  SoCube *g4Box = new SoCube();
  g4Box->width =2*box.GetXHalfLength();
  g4Box->height=2*box.GetYHalfLength();
  g4Box->depth =2*box.GetZHalfLength(); 
  fCurrentSeparator->addChild(g4Box);
#endif
}
void G4OpenInventorSceneHandler::AddThis (const G4Tubs & tubs) {
#ifdef USE_SOPOLYHEDRON
  Geant4_SoPolyhedron* soPolyhedron = 
    new Geant4_SoPolyhedron(tubs.CreatePolyhedron());
  soPolyhedron->solid.setValue(fModelingSolid?TRUE:FALSE);
  fCurrentSeparator->addChild(soPolyhedron);  
#else
  SoTubs *g4Tubs = new SoTubs();
  g4Tubs->pRMin = tubs.GetRMin();
  g4Tubs->pRMax = tubs.GetRMax();
  g4Tubs->pDz = tubs.GetDz();
  g4Tubs->pSPhi = tubs.GetSPhi();
  g4Tubs->pDPhi = tubs.GetDPhi();
  fCurrentSeparator->addChild(g4Tubs);
#endif
}
void G4OpenInventorSceneHandler::AddThis (const G4Cons &cons) {
#ifdef USE_SOPOLYHEDRON
  Geant4_SoPolyhedron* soPolyhedron = 
    new Geant4_SoPolyhedron(cons.CreatePolyhedron());
  soPolyhedron->solid.setValue(fModelingSolid?TRUE:FALSE);
  fCurrentSeparator->addChild(soPolyhedron);  
#else
  SoCons *g4Cons = new SoCons();
  g4Cons->fRmin1 = cons.GetRmin1();
  g4Cons->fRmin2 = cons.GetRmin2();
  g4Cons->fRmax1 = cons.GetRmax1();
  g4Cons->fRmax2 = cons.GetRmax2();
  g4Cons->fDz    = cons.GetDz();
  g4Cons->fSPhi  = cons.GetSPhi();
  g4Cons->fDPhi  = cons.GetDPhi();
  fCurrentSeparator->addChild(g4Cons);
#endif
}

void G4OpenInventorSceneHandler::AddThis (const G4Trap &trap) {
#ifdef USE_SOPOLYHEDRON
  Geant4_SoPolyhedron* soPolyhedron = 
    new Geant4_SoPolyhedron(trap.CreatePolyhedron());
  soPolyhedron->solid.setValue(fModelingSolid?TRUE:FALSE);
  fCurrentSeparator->addChild(soPolyhedron);  
#else
  G4ThreeVector SymAxis=trap.GetSymAxis();
    SoTrap *g4Trap = new SoTrap();
  g4Trap->pDz  = trap.GetZHalfLength();
  g4Trap->pPhi = atan2(SymAxis(kYAxis),SymAxis(kXAxis));
  g4Trap->pTheta = acos(SymAxis(kZAxis));
  g4Trap->pDy1 = trap.GetYHalfLength1();
  g4Trap->pDx1 = trap.GetXHalfLength1();
  g4Trap->pDx2 = trap.GetXHalfLength2();
  g4Trap->pDy2 = trap.GetYHalfLength2();
  g4Trap->pDx3 = trap.GetXHalfLength3();
  g4Trap->pDx4 = trap.GetXHalfLength4();
  fCurrentSeparator->addChild(g4Trap);
#endif
}

void G4OpenInventorSceneHandler::AddThis (const G4Trd &trd) {
#ifdef USE_SOPOLYHEDRON
  Geant4_SoPolyhedron* soPolyhedron = 
    new Geant4_SoPolyhedron(trd.CreatePolyhedron());
  soPolyhedron->solid.setValue(fModelingSolid?TRUE:FALSE);
  fCurrentSeparator->addChild(soPolyhedron);  
#else
  SoTrd *g4Trd = new SoTrd();
  g4Trd->fDx1 = trd.GetXHalfLength1();
  g4Trd->fDx2 = trd.GetXHalfLength2();
  g4Trd->fDy1 = trd.GetYHalfLength1();
  g4Trd->fDy2 = trd.GetYHalfLength2();
  g4Trd->fDz  = trd.GetZHalfLength();
  fCurrentSeparator->addChild(g4Trd);
#endif
}

void G4OpenInventorSceneHandler::PreAddThis
(const G4Transform3D& objectTransformation,
 const G4VisAttributes& visAttribs) {

  // This assumes that it is a "Pre" to the "Add" of a
  // G4VPhysicalVolume.  This is true at present, but beware!  We
  // might have to think about the cast in the following code.

  G4VSceneHandler::PreAddThis (objectTransformation, visAttribs);
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

  // First find the color attributes.
  //
  const G4VisAttributes* pVisAttribs =
    fpViewer -> GetApplicableVisAttributes (&visAttribs);
  const G4Colour& g4Col =  pVisAttribs -> GetColour ();
  const double red=g4Col.GetRed ();
  const double green=g4Col.GetGreen ();
  const double blue=g4Col.GetBlue ();
  if(pVisAttribs->IsForceDrawingStyle()
    &&(pVisAttribs->GetForcedDrawingStyle()==G4VisAttributes::solid)) {
    fModelingSolid = true;
  } else {
    fModelingSolid = false;
  }
  //printf("debug : PreAddThis : %g %g %g : %d\n",
    //red,green,blue,pVisAttribs->IsVisible());
               
  if(!fpCurrentLV || !fpCurrentPV) return; //GB 

  //printf("debug : OIV : LV : %lx %s : %g %g %g\n",
  // fpCurrentLV,
  // fpCurrentLV->GetName().c_str(),
  // red,green,blue);
  //
  // This block of code is executed for non-leaf parts:
  //
  if (fpCurrentLV->GetNoDaughters()!=0 ||
      fpCurrentPV->IsReplicated()) {
    //printf("debug : PreAddThis : container\n");

    //
    // Make the detector tree kit:
    //
    SoDetectorTreeKit* g4DetectorTreeKit = new SoDetectorTreeKit();  

    SoSeparator* previewSeparator   =  
      (SoSeparator*) g4DetectorTreeKit->getPart("previewSeparator",TRUE);
    previewSeparator->renderCaching = SoSeparator::OFF;

    SoSeparator* fullSeparator =  
      (SoSeparator*) g4DetectorTreeKit->getPart("fullSeparator",   TRUE);
    fullSeparator->renderCaching = SoSeparator::OFF;

    SoMaterial* matColor = (SoMaterial*) 
      g4DetectorTreeKit->getPart("appearance.material", TRUE);
    matColor->diffuseColor.setValue(red,green,blue);

    SoLightModel* lightModel = (SoLightModel*)
      g4DetectorTreeKit->getPart("appearance.lightModel", TRUE);
    lightModel->model.setValue(fModelingSolid ?
			       SoLightModel::PHONG:
			       SoLightModel::BASE_COLOR);
   
    //
    // Add the full separator to the dictionary; it is indexed by the 
    // address of the logical volume!
    // NOTE: the map is no longer built iteratively from the whole hierarchy
    //       of volumes since it is no longer possible to retrieve the mother
    //       physical volume. The algorithm requires review !   - GC
    //
    fSeparatorMap[fpCurrentPV->GetLogicalVolume()] = fullSeparator;

    //
    // Find out where to add this volume.  This means locating its mother.  
    // If no mother can be found, it goes under root.
    //
    G4LogicalVolume* MotherVolume = fpCurrentPV->GetMotherLogical();
    if (MotherVolume) {
      if (fSeparatorMap.find(MotherVolume) != fSeparatorMap.end()) {
        //printf("debug : PreAddThis : mother found in map\n");
        fSeparatorMap[MotherVolume]->addChild(g4DetectorTreeKit);
      } else {
        //printf("debug : PreAddThis : mother not found in map !!!\n");
        fDetectorRoot->addChild(g4DetectorTreeKit);
      }
    } else {
      //printf("debug : PreAddThis : has no mother\n");
      fDetectorRoot->addChild(g4DetectorTreeKit);
    }
    fCurrentSeparator = previewSeparator;
  } else {
    //
    // This block of code is executed for leaf parts.
    //
    //
    // Locate the mother volume and find it's corresponding full separator
    //
    fCurrentSeparator = 0;
    G4LogicalVolume* MotherVolume = fpCurrentPV->GetMotherLogical();
    if (MotherVolume) {
      if (fSeparatorMap.find(MotherVolume) != fSeparatorMap.end()) {
        fCurrentSeparator = fSeparatorMap[MotherVolume];
      } else {
        G4cerr << "ERROR - G4OpenInventorSceneHandler::PreAddThis()"
               << "        OIScene leaf protocol error. Mother volume "
               << MotherVolume->GetName() <<  " missing." << G4endl;
        G4cerr << "        Daughter volume was: "
               << fpCurrentPV->GetName()
               << G4endl;
      }
    }
    //
    // If the mother volume has no full separator, then the solid and its 
    // attributes will go under "fDetectorRoot"
    //
    if (!fCurrentSeparator) {
      SoMaterial* newColor = new SoMaterial();
      newColor->diffuseColor.setValue(red,green,blue);
      fDetectorRoot->addChild(newColor);
      SoLightModel* lightModel = new SoLightModel;
      lightModel->model.setValue(fModelingSolid ?
			         SoLightModel::PHONG:
			         SoLightModel::BASE_COLOR);
      fDetectorRoot->addChild(lightModel);
      fCurrentSeparator = fDetectorRoot;
    } else {
      SoMaterial* newColor = new SoMaterial();
      newColor->diffuseColor.setValue(red,green,blue);
      fCurrentSeparator->addChild(newColor);    
      SoLightModel* lightModel = new SoLightModel;
      lightModel->model.setValue(fModelingSolid ?
			         SoLightModel::PHONG:
			         SoLightModel::BASE_COLOR);
      fCurrentSeparator->addChild(lightModel);
    }
  }
  //
  // Set up the geometrical transformation for the coming 
  //
  G4OpenInventorTransform3D oiTran (objectTransformation);
  SoSFMatrix* oiMat = oiTran.GetOIMatrix();
  SoMatrixTransform* xform = new SoMatrixTransform;
  xform->matrix.setValue(oiMat->getValue());
  fCurrentSeparator->addChild(new SoResetTransform);
  fCurrentSeparator->addChild(xform);
}


G4double  G4OpenInventorSceneHandler::GetMarkerSize ( const G4VMarker& mark ) 
{
        //----- return value ( marker radius in 3d units) 
        G4double size       = 1.0 ; // initialization

        //----- parameters to calculate 3d size from 2d size
        const double HALF_SCREEN_SIZE_2D = 300.0 ; // pixels
        double zoom_factor  = fpViewer->GetViewParameters().GetZoomFactor() ;
        if ( zoom_factor <=  0.0 ) { zoom_factor = 1.0 ; }
        double extent_radius_3d = GetScene()->GetExtent().GetExtentRadius() ;
        if ( extent_radius_3d <= 0.0 ) { extent_radius_3d = 1.0 ; } 

        //----- get marker radius in 3D units
        size = mark.GetWorldSize();
        if        ( size  ) {

                // get mark radius in 3D units
                size = 0.5 * mark.GetWorldSize()  ; 

        } else if ( (size = mark.GetScreenSize()) ) {

                // local
                double mark_radius_2d   = 0.5 * mark.GetScreenSize() ;

                // get mark radius in 3D units
                size \
                 = extent_radius_3d * ( mark_radius_2d / HALF_SCREEN_SIZE_2D );
                size *= zoom_factor ;

        } else {
                // local
                double mark_radius_2d \
                 = fpViewer->GetViewParameters().GetDefaultMarker().GetScreenSize();

                // get mark radius in 3D units
                size \
                 = extent_radius_3d * ( mark_radius_2d / HALF_SCREEN_SIZE_2D );
                size *= zoom_factor ;
        }

                //----- global rescaling
        size *= fpViewer->GetViewParameters().GetGlobalMarkerScale(); 

                //----- return size
        return size ;

} // G4OpenInventorSceneHandler::GetMarkerSize ()

void G4OpenInventorSceneHandler::AddThis(const G4VTrajectory& traj) {
  G4VSceneHandler::AddThis(traj);  // For now.
}

void G4OpenInventorSceneHandler::AddThis(const G4VHit& hit) {
  G4VSceneHandler::AddThis(hit);  // For now.
}

#endif
