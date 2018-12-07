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
/// \file Par01/src/Par01EnergySpot.cc
/// \brief Implementation of the Par01EnergySpot class
//
//
//
#include "Par01EnergySpot.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4VVisManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01EnergySpot::Par01EnergySpot()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01EnergySpot::Par01EnergySpot(const G4ThreeVector& point, G4double E)
{
  fPoint = point;
  fEnergy = E;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01EnergySpot::~Par01EnergySpot()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01EnergySpot::Draw(G4Colour *color)
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager)
    {
      G4Polyline polyline;
      G4Colour colour(1.,.5,.5);
      if (color != 0) colour = *color;
      polyline.SetVisAttributes(colour);
      G4ThreeVector pp(fPoint);
      // Draw a "home made" marker:
      // Will be better by using a real Marker:
      pp.setZ(pp.z()+1*cm);
      polyline.push_back(pp);
      pp.setZ(pp.z()-2*cm);
      polyline.push_back(pp);
      pp = fPoint;
      polyline.push_back(pp);
      pp.setX(pp.x()+1*cm);
      polyline.push_back(pp);
      pp.setX(pp.x()-2*cm);
      polyline.push_back(pp);
      pp = fPoint;
      polyline.push_back(pp);
      pp.setY(pp.y()+1*cm);
      polyline.push_back(pp);
      pp.setY(pp.y()-2*cm);
      polyline.push_back(pp);
      pVVisManager -> Draw(polyline);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01EnergySpot::Print()
{
  G4cout << " Par01EnergySpot {E = " << fEnergy << "; Position = " << fPoint << " }"<< G4endl;
}
