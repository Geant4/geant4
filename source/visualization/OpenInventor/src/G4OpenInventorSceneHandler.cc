// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorSceneHandler.cc,v 1.7 2000-04-12 13:09:15 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor stored scene - creates OpenInventor display lists.
#ifdef G4VIS_BUILD_OI_DRIVER

#include <Inventor/SoPath.h>
#include <Inventor/SoNodeKitPath.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoCoordinate4.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/fields/SoMFString.h>
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

#include <HEPVis/nodes/SoG4Box.h>
#include <HEPVis/nodes/SoG4Tubs.h>
#include <HEPVis/nodes/SoG4Cons.h>
#include <HEPVis/nodes/SoG4Trd.h>
#include <HEPVis/nodes/SoG4Trap.h>
#include <HEPVis/nodekits/SoDetectorTreeKit.h>

#include "G4Scene.hh"
#include "G4NURBS.hh"
#include "G4OpenInventor.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4OpenInventorViewer.hh"
#include "G4OpenInventorTransform3D.hh"
#include "G4ThreeVector.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Transform3D.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
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

typedef SoDetectorTreeKit SoG4DetectorTreeKit;
G4Point3D translation;

G4OpenInventorSceneHandler::G4OpenInventorSceneHandler (G4OpenInventor& system,
					  const G4String& name)
:G4VSceneHandler (system, fSceneIdCount++, name)
,root(NULL)
,staticRoot(NULL)
,transientRoot(NULL)
{
  //
  root = new SoSeparator;
  root->ref();
  root->setName("Root");
  //
  staticRoot = new SoSeparator;
  staticRoot->setName("StaticRoot");
  root->addChild(staticRoot);
  //
  transientRoot = new SoSeparator;
  transientRoot->setName("TransientRoot");
  root->addChild(transientRoot);
  //
  fSceneCount++;
}

G4OpenInventorSceneHandler::~G4OpenInventorSceneHandler ()
{
  root->unref();
  fSceneCount--;
}

//
// Method for handling G4Polyline objects (from tracking or wireframe).
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Polyline& line) {
  G4int nPoints = line.entries();
  SbVec3f* pCoords = new SbVec3f[nPoints];

  SoCoordinate3 *polyCoords = new SoCoordinate3;
  SoDrawStyle *drawStyle = new SoDrawStyle;
  SoLineSet *pLine = new SoLineSet;

  for (G4int iPoint = 0; iPoint < nPoints ; iPoint++) {
    pCoords[iPoint].setValue(line(iPoint).x(),
			     line(iPoint).y(),
			     line(iPoint).z());
  }

  //
  // Color
  //
  SoMaterial *pMaterial = new SoMaterial;
  const G4Colour& c = GetColour (line);
  pMaterial->diffuseColor.setValue(c.GetRed (), c.GetGreen (), c.GetBlue ());
  currentSeparator->addChild(pMaterial);
  
  //
  // Point Set
  // 
  polyCoords->point.setValues(0,nPoints,pCoords);
  currentSeparator->addChild(polyCoords);
  
  //
  // Wireframe
  // 
  drawStyle->style = SoDrawStyle::LINES;
  currentSeparator->addChild(drawStyle);
#ifdef INVENTOR2_0
  pLine->numVertices.setValues(0,1,(const long *)&nPoints);
#else 
  pLine->numVertices.setValues(0,1,&nPoints);
#endif
  currentSeparator->addChild(pLine);

  delete [] pCoords;
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
  currentSeparator->addChild(pMaterial);

  //
  // Font
  // 

  SoFont *g4Font = new SoFont();
  g4Font->size = 20.0;
  currentSeparator->addChild(g4Font);

  //
  // Text
  // 

  SoMFString *tString = new SoMFString;
  tString->setValue("This is Text");

  SoText2 *g4String = new SoText2();
  g4String->spacing = 2.0;
  g4String->justification = SoText2::LEFT;
  currentSeparator->addChild(g4String);
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
  currentSeparator->addChild(pMaterial);

  //
  // Dimensions
  //
  G4double userSpecified = circle.GetWorldSize() || circle.GetScreenSize();
  const G4VMarker& defaultMarker =
     fpViewer -> GetViewParameters().GetDefaultMarker();

  G4double size = GetMarkerSize ( circle );

  //
  // Position
  //
  
  const G4Point3D cCtr = circle.GetPosition();
   
  SoTranslation *cTrans = new SoTranslation;
  cTrans->translation.setValue(cCtr.x(),cCtr.y(),cCtr.z());
  currentSeparator->addChild(cTrans);

  //
  // Sphere
  // 

  SoSphere *g4Sphere = new SoSphere();
  g4Sphere->radius = size;
  currentSeparator->addChild(g4Sphere);
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
  currentSeparator->addChild(pMaterial);

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
  currentSeparator->addChild(sTrans);
  
  SoCube *g4Square = new SoCube();
  g4Square->width =  sSize;
  g4Square->height = sSize;
  g4Square->depth =  sSize;
  currentSeparator->addChild(g4Square);
}

//
// Method for handling G4Polyhedron objects for drawing solids.
//
void G4OpenInventorSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {
#define MAXPOLYH	32767
  //
  // Assume all facets are convex quadrilaterals.
  //
  SbVec3f* polyVerts   = new SbVec3f[MAXPOLYH*4];
  G4int*   numPolyVert = new G4int [MAXPOLYH];
  SbVec3f* polyNorms   = new SbVec3f[MAXPOLYH];

  //
  // Loop through the G4Facets...
  //
  G4int vertIdx = 0, polyIdx = 0;

  G4bool notLastFace;
  G4Normal3D SurfaceUnitNormal;

  do {
  
    //
    // First, find surface normal for the facet...
    //
    notLastFace = polyhedron.GetNextUnitNormal (SurfaceUnitNormal);
    
    //
    // Second, set the normal for all vertices in the facet...
    //
    polyNorms[polyIdx].setValue(SurfaceUnitNormal.x(),
				SurfaceUnitNormal.y(),
				SurfaceUnitNormal.z());
    
    //
    // Loop through the four edges of each G4Facet...
    //
    G4int faceIdx = 0;
    G4bool notLastEdge;
    G4Point3D vertex;
    G4int edgeFlag=1;
    do {
      notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
      vertex.transform (*fpObjectTransformation);
      
      polyVerts[vertIdx].setValue(vertex.x(),
				  vertex.y(),
				  vertex.z());

      vertIdx++;
      faceIdx++;
    } while (notLastEdge);

    numPolyVert[polyIdx] = faceIdx;
    polyIdx++;

  } while ((notLastFace)&&(polyIdx < MAXPOLYH));

  //
  // Store Norms
  //
  SoNormal *pNorms = new SoNormal;
  pNorms->vector.setValues(0,polyIdx,polyNorms);
  currentSeparator->addChild(pNorms);
  
  SoNormalBinding *pNormalBinding = new SoNormalBinding;
  pNormalBinding->value = SoNormalBinding::PER_FACE;
  currentSeparator->addChild(pNormalBinding);

  //
  // Define Material
  //
  SoMaterial *pMaterial = new SoMaterial;
  const G4Colour& c = GetColour (polyhedron);
  pMaterial->diffuseColor.setValue(c.GetRed (), c.GetGreen (), c.GetBlue ());
  currentSeparator->addChild(pMaterial);

  //
  // Store Vertex Coordinates  
  //
  SoCoordinate3 *polyCoords = new SoCoordinate3;
  polyCoords->point.setValues(0,vertIdx,polyVerts);
  currentSeparator->addChild(polyCoords);
  
  //
  // Define FaceSet
  //
  SoFaceSet *pFaceSet = new SoFaceSet;
#ifdef INVENTOR2_0
  pFaceSet->numVertices.setValues(0,polyIdx,(const long *)numPolyVert);
#else 
  pFaceSet->numVertices.setValues(0,polyIdx,numPolyVert);
#endif
  currentSeparator->addChild(pFaceSet);  

  //
  // Clean-up
  //
  delete [] polyVerts;
  delete [] numPolyVert;
  delete [] polyNorms;
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
  
  currentSeparator->addChild(surfSep);

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
    // set the destination to "transientRoot"
    //  
    currentSeparator=transientRoot;
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
    currentSeparator->addChild(new SoResetTransform);
    currentSeparator->addChild(xform);
  }
}

void G4OpenInventorSceneHandler::EndPrimitives () {
  G4VSceneHandler::EndPrimitives ();
}

void G4OpenInventorSceneHandler::EndModeling () {
}

void G4OpenInventorSceneHandler::ClearStore () {
  G4VSceneHandler::ClearStore();
  staticRoot->removeAllChildren();
  transientRoot->removeAllChildren();
}

void G4OpenInventorSceneHandler::ClearTransientStore () {
  transientRoot->removeAllChildren();
}

void G4OpenInventorSceneHandler::RequestPrimitives (const G4VSolid& solid) {
    // Stop-gap solution for display List re-use.
    // A proper implementation would use geometry hierarchy.
    //
    G4VSceneHandler::RequestPrimitives (solid);
}

G4int G4OpenInventorSceneHandler::fSceneIdCount = 0;

G4int G4OpenInventorSceneHandler::fSceneCount = 0;

//
// These methods add primitives using the HEPVis classes
//
void G4OpenInventorSceneHandler::AddThis (const G4Box & Box) {
  SoCube *g4Box = new SoCube();
  g4Box->width =2*Box.GetXHalfLength();
  g4Box->height=2*Box.GetYHalfLength();
  g4Box->depth =2*Box.GetZHalfLength(); 
  currentSeparator->addChild(g4Box);
}
void G4OpenInventorSceneHandler::AddThis (const G4Tubs & Tubs) {
  SoG4Tubs *g4Tubs = new SoG4Tubs();
  g4Tubs->pRMin = Tubs.GetRMin();
  g4Tubs->pRMax = Tubs.GetRMax();
  g4Tubs->pDz = Tubs.GetDz();
  g4Tubs->pSPhi = Tubs.GetSPhi();
  g4Tubs->pDPhi = Tubs.GetDPhi();
  currentSeparator->addChild(g4Tubs);
}
void G4OpenInventorSceneHandler::AddThis (const G4Cons &cons) {
  SoG4Cons *g4Cons = new SoG4Cons();
  g4Cons->fRmin1 = cons.GetRmin1();
  g4Cons->fRmin2 = cons.GetRmin2();
  g4Cons->fRmax1 = cons.GetRmax1();
  g4Cons->fRmax2 = cons.GetRmax2();
  g4Cons->fDz    = cons.GetDz();
  g4Cons->fSPhi  = cons.GetSPhi();
  g4Cons->fDPhi  = cons.GetDPhi();
  currentSeparator->addChild(g4Cons);
}

void G4OpenInventorSceneHandler::AddThis (const G4Trap &trap) {
  G4ThreeVector SymAxis=trap.GetSymAxis();
  
  SoG4Trap *g4Trap = new SoG4Trap();
  g4Trap->pDz  = trap.GetZHalfLength();
  g4Trap->pPhi = atan2(SymAxis(kYAxis),SymAxis(kXAxis));
  g4Trap->pTheta = acos(SymAxis(kZAxis));
  g4Trap->pDy1 = trap.GetYHalfLength1();
  g4Trap->pDx1 = trap.GetXHalfLength1();
  g4Trap->pDx2 = trap.GetXHalfLength2();
  g4Trap->pDy2 = trap.GetYHalfLength2();
  g4Trap->pDx3 = trap.GetXHalfLength3();
  g4Trap->pDx4 = trap.GetXHalfLength4();
  currentSeparator->addChild(g4Trap);
}

void G4OpenInventorSceneHandler::AddThis (const G4Trd &trd) {
  SoG4Trd *g4Trd = new SoG4Trd();
  g4Trd->fDx1 = trd.GetXHalfLength1();
  g4Trd->fDx2 = trd.GetXHalfLength2();
  g4Trd->fDy1 = trd.GetYHalfLength1();
  g4Trd->fDy2 = trd.GetYHalfLength2();
  g4Trd->fDz  = trd.GetZHalfLength();
  currentSeparator->addChild(g4Trd);
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
  // we create a detector tree kit.  This has two separators.  The solid itself
  // goes in the preview separator, the full separator is forseen for daughters.
  // This separator is not only created--it is also put in a dictionary for 
  // future use by the leaf part.

  // For leaf parts, we first locate the mother volume and find its separator 
  // through the dictionary.

  // The private member currentSeparator is set to the proper location on in the
  // scene database so that when the solid is actually added (in addthis), it is
  //  put in the right place.

  // First find the color attributes.
  //
  const G4VisAttributes* pVisAttribs =
  fpViewer -> GetApplicableVisAttributes (&visAttribs);
  const G4Colour& g4Col =  pVisAttribs -> GetColour ();
  const double red=g4Col.GetRed ();
  const double green=g4Col.GetGreen ();
  const double blue=g4Col.GetBlue ();
               
  
  //
  // This block of code is executed for non-leaf parts:
  //
  if (fpCurrentLV->GetNoDaughters()!=0 ||
      fpCurrentPV->IsReplicated()) {

    //
    // Make the detector tree kit:
    //
    SoG4DetectorTreeKit* g4DetectorTreeKit = new SoG4DetectorTreeKit();  

    SoSeparator* previewSeparator   =  
      (SoSeparator*) g4DetectorTreeKit->getPart("previewSeparator",TRUE);
    SoSeparator* fullSeparator =  
      (SoSeparator*) g4DetectorTreeKit->getPart("fullSeparator",   TRUE);
    previewSeparator->renderCaching=SoSeparator::OFF;
    fullSeparator->renderCaching=SoSeparator::OFF;
    SoMaterial* matColor = (SoMaterial*) 
      g4DetectorTreeKit->getPart("appearance.material", TRUE);
    matColor->diffuseColor.setValue(red,green,blue);
   
    //
    // Add the full separator to the dictionary; it is indexed by the 
    // address of the physical volume!
    //
    SeparatorMap[fpCurrentPV]=fullSeparator;

    //
    // Find out where to add this volume.  This means locating its mother.  
    // If no mother can be found, it goes under root.
    //
    G4VPhysicalVolume* MotherVolume = fpCurrentPV->GetMother();
    while (MotherVolume) {
      if (SeparatorMap.find(MotherVolume) != SeparatorMap.end()) {
	SeparatorMap[MotherVolume]->addChild(g4DetectorTreeKit);
        break;
      }
      else {
        G4cerr << "OIScene non-leaf protocol error.  Mother volume " << 
	          MotherVolume->GetName() << " missing." << G4endl;
        G4cerr << "                         Daughter volume was: "
	     << fpCurrentPV->GetName()
	     << G4endl;
      }
      MotherVolume=MotherVolume->GetMother();
    }
    if (!MotherVolume) staticRoot->addChild(g4DetectorTreeKit);
    currentSeparator = previewSeparator;
  }
  //
  // This block of code is executed for leaf parts.
  //
  else {
    //
    // Locate the mother volume and find it's corresponding full separator
    //
    currentSeparator = NULL;
    G4VPhysicalVolume* MotherVolume = fpCurrentPV->GetMother();
    while (MotherVolume) {
      if (SeparatorMap.find(MotherVolume) != SeparatorMap.end()) {
	currentSeparator=SeparatorMap[MotherVolume];
        break;
      }
      else {
	G4cerr << "OIScene leaf protocol error.  Mother volume " << 
	          MotherVolume->GetName() <<  " missing." << G4endl;
	G4cerr << "                         Daughter volume was: "
	     << fpCurrentPV->GetName()
	     << G4endl;
      }
      MotherVolume = MotherVolume->GetMother();
    }
    //
    // If the mother volume has no full separator, then the solid and its 
    // attributes will go under "staticRoot"
    //
    if (!currentSeparator) {
      SoMaterial* newColor = new SoMaterial();
      newColor->diffuseColor.setValue(red,green,blue);
      staticRoot->addChild(newColor);
      currentSeparator=staticRoot;
    }
    else {
      SoMaterial* newColor = new SoMaterial();
      newColor->diffuseColor.setValue(red,green,blue);
      currentSeparator->addChild(newColor);    
    }
  }
  //
  // Set up the geometrical transformation for the coming 
  //
  G4OpenInventorTransform3D oiTran (objectTransformation);
  SoSFMatrix* oiMat = oiTran.GetOIMatrix();
  SoMatrixTransform* xform = new SoMatrixTransform;
  xform->matrix.setValue(oiMat->getValue());
  currentSeparator->addChild(new SoResetTransform);
  currentSeparator->addChild(xform);
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
	if        ( size = mark.GetWorldSize()  ) {

		// get mark radius in 3D units
		size = 0.5 * mark.GetWorldSize()  ; 

	} else if ( size = mark.GetScreenSize() ) {

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
#endif
