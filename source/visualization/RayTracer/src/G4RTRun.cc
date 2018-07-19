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
// $Id: G4RTRun.cc 66264 2012-12-14 10:17:44Z allison $
//
//
//

///////////////////////
//G4RTRun.cc
///////////////////////

#include "G4RTRun.hh"
#include "G4TheMTRayTracer.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"

#include "G4RayTrajectory.hh"
#include "G4RayTrajectoryPoint.hh"

G4RTRun::G4RTRun()
{
  colorMap = new G4THitsMap<G4Colour>("G4RTRun","ColorMap");

  backgroundColour = G4TheMTRayTracer::theInstance->backgroundColour;
  lightDirection = G4TheMTRayTracer::theInstance->lightDirection;
  attenuationLength = G4TheMTRayTracer::theInstance->attenuationLength;
}

G4RTRun::~G4RTRun()
{ 
  colorMap->clear();
  delete colorMap;
}

void G4RTRun::RecordEvent(const G4Event* evt)
{
  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  if(!trajectoryContainer) return;
  G4RayTrajectory* trajectory = static_cast<G4RayTrajectory*>( (*trajectoryContainer)[0] );
  if(!trajectory) return;

  G4int nPoint = trajectory->GetPointEntries();
  if(nPoint==0) return;

  G4int evId = evt->GetEventID();
  G4Colour initialCol(backgroundColour);
  if( trajectory->GetPointC(nPoint-1)->GetPostStepAtt() )
  { initialCol = GetSurfaceColour(trajectory->GetPointC(nPoint-1)); }
  G4Colour rayColour = Attenuate(trajectory->GetPointC(nPoint-1),initialCol);

  for(int i=nPoint-2;i>=0;i--)
  {
    G4Colour surfaceCol = GetSurfaceColour(trajectory->GetPointC(i));
    G4double weight = 1.0 - surfaceCol.GetAlpha();
    G4Colour mixedCol = GetMixedColour(rayColour,surfaceCol,weight);
    rayColour = Attenuate(trajectory->GetPointC(i),mixedCol);
  }

  colorMap->add(evId,rayColour);
}

void G4RTRun::Merge(const G4Run* aLocalRun)
{
  const G4RTRun* theLocalRun = static_cast<const G4RTRun*>(aLocalRun);
  if(theLocalRun) *(colorMap) += *(theLocalRun->colorMap);
  G4Run::Merge(aLocalRun);
}

G4Colour G4RTRun::GetSurfaceColour(G4RayTrajectoryPoint* point)
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
    G4double red   = preAtt->GetColour().GetRed();
    G4double green = preAtt->GetColour().GetGreen();
    G4double blue  = preAtt->GetColour().GetBlue();
    preCol = G4Colour
      (red*brill,green*brill,blue*brill,preAtt->GetColour().GetAlpha());
  }
  else
  { preCol = transparent; }

  if(postVis)
  {
    G4double brill = (1.0-(-lightDirection).dot(-normal))/2.0;
    G4double red   = postAtt->GetColour().GetRed();
    G4double green = postAtt->GetColour().GetGreen();
    G4double blue  = postAtt->GetColour().GetBlue();
    postCol = G4Colour
      (red*brill,green*brill,blue*brill,postAtt->GetColour().GetAlpha());
  }
  else
  { postCol = transparent; }

  if(!preVis) return postCol;
  if(!postVis) return preCol;

  G4double weight = 0.5;
  return GetMixedColour(preCol,postCol,weight);
}

G4Colour G4RTRun::GetMixedColour(G4Colour surfCol,G4Colour transCol,G4double weight)
{
  G4double red   = weight*surfCol.GetRed() + (1.-weight)*transCol.GetRed();
  G4double green = weight*surfCol.GetGreen() + (1.-weight)*transCol.GetGreen();
  G4double blue  = weight*surfCol.GetBlue() + (1.-weight)*transCol.GetBlue();
  G4double alpha = weight*surfCol.GetAlpha() + (1.-weight)*transCol.GetAlpha();
  return G4Colour(red,green,blue,alpha);
}

G4Colour G4RTRun::Attenuate(G4RayTrajectoryPoint* point, G4Colour sourceCol)
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

G4bool G4RTRun::ValidColour(const G4VisAttributes* visAtt)
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


