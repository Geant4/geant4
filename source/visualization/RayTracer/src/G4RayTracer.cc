// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayTracer.cc,v 1.6 2000-06-07 02:52:45 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//


#include "G4RayTracer.hh"
#include "G4RayTracerFeatures.hh"
#include "G4RayTracerSceneHandler.hh"
#include "G4RayTracerViewer.hh"
#include "G4EventManager.hh"
#include "G4RTMessenger.hh"
#include "G4RayShooter.hh"
#include "G4VFigureFileMaker.hh"
#include "G4RTTrackingAction.hh"
#include "G4RTSteppingAction.hh"
#include "G4RayTrajectory.hh"
#include "G4RayTrajectoryPoint.hh"
#include "G4RTJpegMaker.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"
#include "G4TransportationManager.hh"

G4RayTracer::G4RayTracer(G4VFigureFileMaker* figMaker)
:G4VGraphicsSystem("RayTracer","RayTracer",RAYTRACER_FEATURES,
                     G4VGraphicsSystem::threeD)
{
  theFigMaker = figMaker;
  if(!theFigMaker) theFigMaker = new G4RTJpegMaker();
  theRayShooter = new G4RayShooter();
  theRayTracerEventAction = 0;
  theRayTracerStackingAction = 0;
  theRayTracerTrackingAction = new G4RTTrackingAction();
  theRayTracerSteppingAction = new G4RTSteppingAction();
  theMessenger = new G4RTMessenger(this,theRayTracerSteppingAction);
  theEventManager = G4EventManager::GetEventManager();

  nColumn = 640;
  nRow = 640;

  eyePosition = G4ThreeVector(1.*m,1.*m,1.*m);
  targetPosition = G4ThreeVector(0.,0.,0.);
  lightDirection = G4ThreeVector(-0.1,-0.2,-0.3).unit();
  viewSpan = 5.0*deg;
  headAngle = 270.*deg; 
  attenuationLength = 1.0*m;

  distortionOn = false;
}

G4RayTracer::~G4RayTracer()
{
  delete theRayShooter;
  delete theRayTracerTrackingAction;
  delete theRayTracerSteppingAction;
  delete theMessenger;
}

void G4RayTracer::Trace(G4String fileName)
{
  G4StateManager* theStateMan = G4StateManager::GetStateManager();
  G4ApplicationState currentState = theStateMan->GetCurrentState();
  if(currentState!=Idle)
  {
    G4cerr << "Illegal application state - Trace() ignored." << G4endl;
    return;
  }

  if(!theFigMaker)
  {
    G4cerr << "Figure file maker class is not specified - Trace() ignored." << G4endl;
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
  { G4cerr << "Could not create figure file" << G4endl;
    G4cerr << "You might set the eye position outside of the world volume" << G4endl; }
  RestoreUserActions();

  if(storeTrajectory==0) UI->ApplyCommand("/tracking/storeTrajectory 0");

  delete [] colorR;
  delete [] colorG;
  delete [] colorB;
}

G4VSceneHandler* G4RayTracer::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4RayTracerSceneHandler (*this, name);
  G4cout << G4RayTracerSceneHandler::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << G4endl;
  return pScene;
}

G4VViewer* G4RayTracer::CreateViewer (G4VSceneHandler& sceneHandler,
				      const G4String& name) {
  G4VViewer* pView = new G4RayTracerViewer (sceneHandler, name);
  return pView;
}

void G4RayTracer::StoreUserActions()
{ 
  theUserEventAction = theEventManager->GetUserEventAction();
  theUserStackingAction = theEventManager->GetUserStackingAction();
  theUserTrackingAction = theEventManager->GetUserTrackingAction();
  theUserSteppingAction = theEventManager->GetUserSteppingAction();

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

void G4RayTracer::RestoreUserActions()
{
  theEventManager->SetUserAction(theUserEventAction);
  theEventManager->SetUserAction(theUserStackingAction);
  theEventManager->SetUserAction(theUserTrackingAction);
  theEventManager->SetUserAction(theUserSteppingAction);

  G4SDManager* theSDMan = G4SDManager::GetSDMpointerIfExist();
  if(theSDMan)
  { theSDMan->Activate("/",true); }
}

G4bool G4RayTracer::CreateBitMap()
{
  G4int iEvent = 0;
  G4double stepAngle = viewSpan/100.;
  G4double viewSpanX = stepAngle*nColumn;
  G4double viewSpanY = stepAngle*nRow;
  G4bool succeeded;
  for(int iRow=0;iRow<nRow;iRow++)
  {
    for(int iColumn=0;iColumn<nColumn;iColumn++)
    {
      G4Event* anEvent = new G4Event(iEvent++);
      
      G4double angleX = -(viewSpanX/2. - iColumn*stepAngle);
      G4double angleY = viewSpanY/2. - iRow*stepAngle;
      G4ThreeVector rayDirection;
      if(distortionOn)
      {
        rayDirection = G4ThreeVector(tan(angleX)/cos(angleY),-tan(angleY)/cos(angleX),1.0);
      }
      else
      {
        rayDirection = G4ThreeVector(tan(angleX),-tan(angleY),1.0);
      }
      rayDirection.rotateZ(headAngle);
      rayDirection.rotateUz(eyeDirection);
      G4ThreeVector rayPosition(eyePosition);
      G4bool interceptable = true;
      // Check if rayPosition is in the world.
      G4VPhysicalVolume* pWorld =
	G4TransportationManager::GetTransportationManager()->
	GetNavigatorForTracking()->GetWorldVolume ();
      EInside whereisit =
	pWorld->GetLogicalVolume()->GetSolid()->Inside(rayPosition);
      if (whereisit != kInside) {
	// It's outside the world, so move it inside.
	G4double outsideDistance =
	  pWorld->GetLogicalVolume()->GetSolid()->
	  DistanceToIn(rayPosition,rayDirection);  
	if (outsideDistance != kInfinity) {
	  rayPosition = rayPosition+outsideDistance*rayDirection;
	}
	else {
	  interceptable = false;
	}
      }
      if (interceptable) {
	theRayShooter->Shoot(anEvent,rayPosition,rayDirection);
	theEventManager->ProcessOneEvent(anEvent);
	succeeded = GenerateColour(anEvent);
  //G4cout << iColumn << " " << iRow << " " << anEvent->GetEventID() << G4endl;
      }
      else {  // Ray does not intercept world at all.
	// Generate background colour...
	G4int iEvent = anEvent->GetEventID();
	colorR[iEvent] = 0;
	colorG[iEvent] = 0;
	colorB[iEvent] = 0;
	succeeded = true;
      }
      delete anEvent;
      if(!succeeded) return false;
    }
  }
  return true;
}

void G4RayTracer::CreateFigureFile(G4String fileName)
{
  //G4cout << nColumn << " " << nRow << G4endl;
  theFigMaker->CreateFigureFile(fileName,nColumn,nRow,colorR,colorG,colorB);
}

G4bool G4RayTracer::GenerateColour(G4Event* anEvent)
{
  G4TrajectoryContainer * trajectoryContainer = anEvent->GetTrajectoryContainer();
  
  G4RayTrajectory* trajectory = (G4RayTrajectory*)( (*trajectoryContainer)[0] );
  if(!trajectory) return false;

  G4int nPoint = trajectory->GetPointEntries();
  if(nPoint==0) return false;

  G4Colour rayColour;
  G4Colour initialColour(1.,1.,1.);
  if( trajectory->GetPointC(nPoint-1)->GetPostStepAtt() )
  { initialColour = GetSurfaceColour(trajectory->GetPointC(nPoint-1)); }
  rayColour = Attenuate(trajectory->GetPointC(nPoint-1),initialColour);

  for(int i=nPoint-2;i>=0;i--)
  {
    G4Colour surfaceColour = GetSurfaceColour(trajectory->GetPointC(i));
    G4double weight = 1.0 - surfaceColour.GetAlpha();
    G4Colour mixedColour = GetMixedColour(rayColour,surfaceColour,weight);
    rayColour = Attenuate(trajectory->GetPointC(i),mixedColour);
  }
    
  G4int iEvent = anEvent->GetEventID();
  colorR[iEvent] = (unsigned char)(int(255*rayColour.GetRed()));
  colorG[iEvent] = (unsigned char)(int(255*rayColour.GetGreen()));
  colorB[iEvent] = (unsigned char)(int(255*rayColour.GetBlue()));
  return true;
}

G4Colour G4RayTracer::GetMixedColour(G4Colour surfCol,G4Colour transCol,G4double weight)
{
  G4double r = weight*surfCol.GetRed() + (1.-weight)*transCol.GetRed();
  G4double g = weight*surfCol.GetGreen() + (1.-weight)*transCol.GetGreen();
  G4double b = weight*surfCol.GetBlue() + (1.-weight)*transCol.GetBlue();
  G4double a = weight*surfCol.GetAlpha() + (1.-weight)*transCol.GetAlpha();
  return G4Colour(r,g,b,a);
}

G4Colour G4RayTracer::GetSurfaceColour(G4RayTrajectoryPoint* point)
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
    G4double brill = (1.0-(-lightDirection).dot(normal))/2.0;
    G4double r = preAtt->GetColour().GetRed();
    G4double g = preAtt->GetColour().GetGreen();
    G4double b = preAtt->GetColour().GetBlue();
    preCol = G4Colour(r*brill,g*brill,b*brill,preAtt->GetColour().GetAlpha());
  }
  else
  { preCol = transparent; }

  if(postVis)
  {
    G4double brill = (1.0-(-lightDirection).dot(-normal))/2.0;
    G4double r = postAtt->GetColour().GetRed();
    G4double g = postAtt->GetColour().GetGreen();
    G4double b = postAtt->GetColour().GetBlue();
    postCol = G4Colour(r*brill,g*brill,b*brill,postAtt->GetColour().GetAlpha());
  }
  else
  { postCol = transparent; }
    
  if(!preVis) return postCol;
  if(!postVis) return preCol;

  G4double weight = 0.5;
  return GetMixedColour(preCol,postCol,weight);
}

G4Colour G4RayTracer::Attenuate(G4RayTrajectoryPoint* point, G4Colour sourceCol)
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

  G4double attenuationFuctor = -stepAlpha/(1.0-stepAlpha)*stepLength/attenuationLength;
  G4double KtRed = exp((1.0-stepRed)*attenuationFuctor);
  G4double KtGreen = exp((1.0-stepGreen)*attenuationFuctor);
  G4double KtBlue = exp((1.0-stepBlue)*attenuationFuctor);
  if(KtRed>1.0){KtRed=1.0;}
  if(KtGreen>1.0){KtGreen=1.0;}
  if(KtBlue>1.0){KtBlue=1.0;}
  return G4Colour(sourceCol.GetRed()*KtRed,
    sourceCol.GetGreen()*KtGreen,sourceCol.GetBlue()*KtBlue);
}

G4bool G4RayTracer::ValidColour(const G4VisAttributes* visAtt)
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

