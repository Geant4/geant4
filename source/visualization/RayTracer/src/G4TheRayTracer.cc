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
//
//


#include "G4TheRayTracer.hh"
#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4RTMessenger.hh"
#include "G4RayShooter.hh"
#include "G4VFigureFileMaker.hh"
#include "G4RTTrackingAction.hh"
#include "G4RTSteppingAction.hh"
#include "G4RayTrajectory.hh"
#include "G4RayTrajectoryPoint.hh"
#include "G4RTJpegMaker.hh"
#include "G4RTSimpleScanner.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"
#include "G4TransportationManager.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VVisManager.hh"

#define G4warn G4cout

G4VFigureFileMaker * G4TheRayTracer::theFigMaker = 0;
G4VRTScanner * G4TheRayTracer::theScanner = 0;

G4TheRayTracer::G4TheRayTracer(G4VFigureFileMaker* figMaker,
			       G4VRTScanner* scanner)
{
  theFigMaker = figMaker;
  if(!theFigMaker) theFigMaker = new G4RTJpegMaker;
  theScanner = scanner;
  if(!theScanner) theScanner = new G4RTSimpleScanner;
  theRayShooter = new G4RayShooter();
  theUserEventAction = 0;
  theUserStackingAction = 0;
  theUserTrackingAction = 0;
  theUserSteppingAction = 0;
  theRayTracerEventAction = 0;
  theRayTracerStackingAction = 0;
  theRayTracerTrackingAction = 0;
  theRayTracerSteppingAction = 0;
  colorR = 0;
  colorG = 0;
  colorB = 0;

  theMessenger = G4RTMessenger::GetInstance(this);
  theEventManager = G4EventManager::GetEventManager();

  nColumn = 640;
  nRow = 640;

  eyePosition = G4ThreeVector(1.*m,1.*m,1.*m);
  targetPosition = G4ThreeVector(0.,0.,0.);
  lightDirection = G4ThreeVector(-0.1,-0.2,-0.3).unit();
  up = G4ThreeVector(0,1,0);
  viewSpan = 5.0*deg;
  headAngle = 0.;
  attenuationLength = 1.0*m;

  distortionOn = false;
  antialiasingOn = false;

  backgroundColour = G4Colour(1.,1.,1.);
}

G4TheRayTracer::~G4TheRayTracer()
{
  delete theRayShooter;
  if(theRayTracerTrackingAction) delete theRayTracerTrackingAction;
  if(theRayTracerSteppingAction) delete theRayTracerSteppingAction;
  delete theMessenger;
  delete theScanner;
  delete theFigMaker;
}

void G4TheRayTracer::Trace(const G4String& fileName)
{
  G4StateManager* theStateMan = G4StateManager::GetStateManager();
  G4ApplicationState currentState = theStateMan->GetCurrentState();
  if(currentState!=G4State_Idle)
  {
    G4warn << "Illegal application state - Trace() ignored." << G4endl;
    return;
  }

  if(!theFigMaker)
  {
    G4warn << "Figure file maker class is not specified - Trace() ignored." << G4endl;
    return;
  }

  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int storeTrajectory = UI->GetCurrentIntValue("/tracking/storeTrajectory");
  if(storeTrajectory==0) UI->ApplyCommand("/tracking/storeTrajectory 1");


  G4ThreeVector tmpVec = targetPosition - eyePosition;
  eyeDirection = tmpVec.unit();
  colorR = new unsigned char[nColumn*nRow];
  colorG = new unsigned char[nColumn*nRow];
  colorB = new unsigned char[nColumn*nRow];

  StoreUserActions();
  G4bool succeeded = CreateBitMap();
  if(succeeded)
  { CreateFigureFile(fileName); }
  else
  { G4warn << "Could not create figure file" << G4endl;
    G4warn << "You might set the eye position outside of the world volume" << G4endl; }
  RestoreUserActions();

  if(storeTrajectory==0) UI->ApplyCommand("/tracking/storeTrajectory 0");

  delete [] colorR;
  delete [] colorG;
  delete [] colorB;
}

void G4TheRayTracer::StoreUserActions()
{ 
  theUserEventAction = theEventManager->GetUserEventAction();
  theUserStackingAction = theEventManager->GetUserStackingAction();
  theUserTrackingAction = theEventManager->GetUserTrackingAction();
  theUserSteppingAction = theEventManager->GetUserSteppingAction();

  if(!theRayTracerTrackingAction) theRayTracerTrackingAction = new G4RTTrackingAction();
  if(!theRayTracerSteppingAction) theRayTracerSteppingAction = new G4RTSteppingAction();

  theEventManager->SetUserAction(theRayTracerEventAction);
  theEventManager->SetUserAction(theRayTracerStackingAction);
  theEventManager->SetUserAction(theRayTracerTrackingAction);
  theEventManager->SetUserAction(theRayTracerSteppingAction);

  G4SDManager* theSDMan = G4SDManager::GetSDMpointerIfExist();
  if(theSDMan)
  { theSDMan->Activate("/",false); }

  G4GeometryManager* theGeomMan = G4GeometryManager::GetInstance();
  theGeomMan->OpenGeometry();
  theGeomMan->CloseGeometry(true);
}

void G4TheRayTracer::RestoreUserActions()
{
  theEventManager->SetUserAction(theUserEventAction);
  theEventManager->SetUserAction(theUserStackingAction);
  theEventManager->SetUserAction(theUserTrackingAction);
  theEventManager->SetUserAction(theUserSteppingAction);

  G4SDManager* theSDMan = G4SDManager::GetSDMpointerIfExist();
  if(theSDMan)
  { theSDMan->Activate("/",true); }
}

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Geantino.hh"

G4bool G4TheRayTracer::CreateBitMap()
{
  G4int iEvent = 0;
  G4double stepAngle = viewSpan/100.;
  G4double viewSpanX = stepAngle*nColumn;
  G4double viewSpanY = stepAngle*nRow;
  G4bool succeeded;

  G4VVisManager* visMan = G4VVisManager::GetConcreteInstance();
  visMan->IgnoreStateChanges(true);

// Confirm process(es) of Geantino is initialized
  G4VPhysicalVolume* pWorld =
	G4TransportationManager::GetTransportationManager()->
	GetNavigatorForTracking()->GetWorldVolume();
  G4RegionStore::GetInstance()->UpdateMaterialList(pWorld);
  G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable(pWorld);
  G4ProcessVector* pVector
    = G4Geantino::GeantinoDefinition()->GetProcessManager()->GetProcessList();
  for (G4int j=0; j < (G4int)pVector->size(); ++j) {
      (*pVector)[j]->BuildPhysicsTable(*(G4Geantino::GeantinoDefinition()));
  }

// Close geometry and set the application state
  G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
  geomManager->OpenGeometry();
  geomManager->CloseGeometry(1,0);
  
  G4ThreeVector center(0,0,0);
  G4Navigator* navigator =
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  navigator->LocateGlobalPointAndSetup(center,0,false);

  G4StateManager* theStateMan = G4StateManager::GetStateManager();
  theStateMan->SetNewState(G4State_GeomClosed); 

// Event loop
  theScanner->Initialize(nRow,nColumn);
  G4int iRow, iColumn;
  while (theScanner->Coords(iRow,iColumn)) {
      G4int iCoord = iRow * nColumn + iColumn;
      G4double dRow = 0, dColumn = 0;  // Antialiasing increments.
      G4Event* anEvent = new G4Event(iEvent++);
      G4double angleX = -(viewSpanX/2. - (iColumn+dColumn)*stepAngle);
      G4double angleY = viewSpanY/2. - (iRow+dRow)*stepAngle;
      G4ThreeVector rayDirection;
      if(distortionOn)
      {
	rayDirection = G4ThreeVector(-std::tan(angleX)/std::cos(angleY),std::tan(angleY)/std::cos(angleX),1.0);
      }
      else
      {
	rayDirection = G4ThreeVector(-std::tan(angleX),std::tan(angleY),1.0);
      }
      G4double cp = std::cos(eyeDirection.phi());
      G4double sp = std::sqrt(1.-cp*cp);
      G4double ct = std::cos(eyeDirection.theta());
      G4double st = std::sqrt(1.-ct*ct);
      G4double gamma = std::atan2(ct*cp*up.x()+ct*sp*up.y()-st*up.z(), -sp*up.x()+cp*up.y());
      rayDirection.rotateZ(-gamma);
      rayDirection.rotateZ(headAngle);
      rayDirection.rotateUz(eyeDirection);
      G4ThreeVector rayPosition(eyePosition);
      G4bool interceptable = true;
      // Check if rayPosition is in the world.
      EInside whereisit =
	pWorld->GetLogicalVolume()->GetSolid()->Inside(rayPosition);
      if (whereisit != kInside) {
	// It's outside the world, so move it inside.
	G4double outsideDistance =
	  pWorld->GetLogicalVolume()->GetSolid()->
	  DistanceToIn(rayPosition,rayDirection);  
	if (outsideDistance != kInfinity) {
	  // Borrowing from geometry, where 1e-8 < epsilon < 1e-3, in
	  // absolute/internal length units, is used for ensuring good
	  // behaviour, choose to add 0.001 to ensure rayPosition is
	  // definitely inside the world volume (JA 16/9/2005)...
	  rayPosition = rayPosition+(outsideDistance+0.001)*rayDirection;
	}
	else {
	  interceptable = false;
	}
      }
      if (interceptable) {
	theRayShooter->Shoot(anEvent,rayPosition,rayDirection.unit());
	theEventManager->ProcessOneEvent(anEvent);
	succeeded = GenerateColour(anEvent);
	colorR[iCoord] = (unsigned char)(G4int(255*rayColour.GetRed()));
	colorG[iCoord] = (unsigned char)(G4int(255*rayColour.GetGreen()));
	colorB[iCoord] = (unsigned char)(G4int(255*rayColour.GetBlue()));
      } else {  // Ray does not intercept world at all.
	// Store background colour...
	colorR[iCoord] = (unsigned char)(G4int(255*backgroundColour.GetRed()));
	colorG[iCoord] = (unsigned char)(G4int(255*backgroundColour.GetGreen()));
	colorB[iCoord] = (unsigned char)(G4int(255*backgroundColour.GetBlue()));
	succeeded = true;
      }

      theScanner->Draw(colorR[iCoord],colorG[iCoord],colorB[iCoord]);

      delete anEvent;
      if(!succeeded) return false;
  }

  theStateMan->SetNewState(G4State_Idle); 
  visMan->IgnoreStateChanges(false);
  return true;
}

void G4TheRayTracer::CreateFigureFile(const G4String& fileName)
{
  //G4cout << nColumn << " " << nRow << G4endl;
  theFigMaker->CreateFigureFile(fileName,nColumn,nRow,colorR,colorG,colorB);
}

G4bool G4TheRayTracer::GenerateColour(G4Event* anEvent)
{
  G4TrajectoryContainer * trajectoryContainer = anEvent->GetTrajectoryContainer();
  
  G4RayTrajectory* trajectory = (G4RayTrajectory*)( (*trajectoryContainer)[0] );
  if(!trajectory) return false;

  G4int nPoint = trajectory->GetPointEntries();
  if(nPoint==0) return false;

  G4Colour initialColour(backgroundColour);
  if( trajectory->GetPointC(nPoint-1)->GetPostStepAtt() )
  { initialColour = GetSurfaceColour(trajectory->GetPointC(nPoint-1)); }
  rayColour = Attenuate(trajectory->GetPointC(nPoint-1),initialColour);

  for(G4int i=nPoint-2;i>=0;--i)
  {
    G4Colour surfaceColour = GetSurfaceColour(trajectory->GetPointC(i));
    G4double weight = 1.0 - surfaceColour.GetAlpha();
    G4Colour mixedColour = GetMixedColour(rayColour,surfaceColour,weight);
    rayColour = Attenuate(trajectory->GetPointC(i),mixedColour);
  }
    
  return true;
}

G4Colour G4TheRayTracer::GetMixedColour
(const G4Colour& surfCol,const G4Colour& transCol,G4double weight)
{
  G4double red   = weight*surfCol.GetRed() + (1.-weight)*transCol.GetRed();
  G4double green = weight*surfCol.GetGreen() + (1.-weight)*transCol.GetGreen();
  G4double blue  = weight*surfCol.GetBlue() + (1.-weight)*transCol.GetBlue();
  G4double alpha = weight*surfCol.GetAlpha() + (1.-weight)*transCol.GetAlpha();
  return G4Colour(red,green,blue,alpha);
}

G4Colour G4TheRayTracer::GetSurfaceColour(G4RayTrajectoryPoint* point)
{
  const G4VisAttributes* preAtt = point->GetPreStepAtt();
  const G4VisAttributes* postAtt = point->GetPostStepAtt();

  G4bool preVis = ValidColour(preAtt);
  G4bool postVis = ValidColour(postAtt);

  G4Colour transparent(1.,1.,1.,0.);

  if(!preVis&&!postVis) return transparent;

  G4ThreeVector normal = point->GetSurfaceNormal();

  G4Colour preCol(1.,1.,1.);
  G4Colour postCol(1.,1.,1.);

  if(preVis)
  {
    const G4Colour& preAttColour = preAtt->GetColour();
    G4double brill = (1.0-(-lightDirection).dot(normal))/2.0;
    G4double red   = preAttColour.GetRed();
    G4double green = preAttColour.GetGreen();
    G4double blue  = preAttColour.GetBlue();
    preCol = G4Colour
      (red*brill,green*brill,blue*brill,preAttColour.GetAlpha());
  }
  else
  { preCol = transparent; }

  if(postVis)
  {
    const G4Colour& postAttColour = postAtt->GetColour();
    G4double brill = (1.0-(-lightDirection).dot(-normal))/2.0;
    G4double red   = postAttColour.GetRed();
    G4double green = postAttColour.GetGreen();
    G4double blue  = postAttColour.GetBlue();
    postCol = G4Colour
      (red*brill,green*brill,blue*brill,postAttColour.GetAlpha());
  }
  else
  { postCol = transparent; }
    
  if(!preVis) return postCol;
  if(!postVis) return preCol;

  G4double weight = 0.5;
  return GetMixedColour(preCol,postCol,weight);
}

G4Colour G4TheRayTracer::Attenuate
(G4RayTrajectoryPoint* point,const G4Colour& sourceCol)
{
  const G4VisAttributes* preAtt = point->GetPreStepAtt();

  G4bool visible = ValidColour(preAtt);
  if(!visible) return sourceCol;

  G4Colour objCol = preAtt->GetColour();
  G4double stepRed = objCol.GetRed();
  G4double stepGreen = objCol.GetGreen();
  G4double stepBlue = objCol.GetBlue();
  G4double stepAlpha = objCol.GetAlpha();
  G4double stepLength = point->GetStepLength();

  G4double attenuationFuctor;
  if(stepAlpha > 0.9999999){ stepAlpha = 0.9999999; } // patch to the next line
    attenuationFuctor = -stepAlpha/(1.0-stepAlpha)*stepLength/attenuationLength;
 
  G4double KtRed = std::exp((1.0-stepRed)*attenuationFuctor);
  G4double KtGreen = std::exp((1.0-stepGreen)*attenuationFuctor);
  G4double KtBlue = std::exp((1.0-stepBlue)*attenuationFuctor);
  if(KtRed>1.0){KtRed=1.0;}
  if(KtGreen>1.0){KtGreen=1.0;}
  if(KtBlue>1.0){KtBlue=1.0;}
  return G4Colour(sourceCol.GetRed()*KtRed,
    sourceCol.GetGreen()*KtGreen,sourceCol.GetBlue()*KtBlue);
}

G4bool G4TheRayTracer::ValidColour(const G4VisAttributes* visAtt)
{
  G4bool val = true;
  if(!visAtt)
  { val = false; }
  else if(!(visAtt->IsVisible()))
  { val = false; }
  else if(visAtt->IsForceDrawingStyle()
    &&(visAtt->GetForcedDrawingStyle()==G4VisAttributes::wireframe))
  { val = false; }
  return val;
}

