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
// John Allison  26th August 1998.
// Model which knows how to draw GEANT4 hits.

#include "G4HitsModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4Event.hh"

G4HitsModel::~G4HitsModel () {}

G4HitsModel::G4HitsModel ():
  fpCurrentHit(nullptr)
{
  fType = "G4HitsModel";
  fGlobalTag = "G4HitsModel for all hits.";
  fGlobalDescription = fGlobalTag;
}

void G4HitsModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
  const G4Event* event = fpMP->GetEvent();
  if (event) {
    G4HCofThisEvent* HCE = event->GetHCofThisEvent();
    if (HCE) {
      G4int nHC = (G4int)HCE->GetCapacity();
      for (G4int iHC = 0; iHC < nHC; ++iHC) {
        G4VHitsCollection* HC = HCE -> GetHC (iHC);
        if (HC) {
          for(std::size_t iHit = 0; iHit < HC->GetSize(); ++iHit) {
            fpCurrentHit = HC->GetHit(iHit);
            if (fpCurrentHit) sceneHandler.AddCompound(*fpCurrentHit);
          }
        }
      }
    }
  }
}
