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
// Model which knows how to draw GEANT4 digis (digitisations).

#include "G4DigiModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4Event.hh"

G4DigiModel::~G4DigiModel () {}

G4DigiModel::G4DigiModel ():
  fpCurrentDigi(nullptr)
{
  fType = "G4DigiModel";
  fGlobalTag = "G4DigiModel for all digis.";
  fGlobalDescription = fGlobalTag;
}

void G4DigiModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
  const G4Event* event = fpMP->GetEvent();
  if (event) {
    G4DCofThisEvent* DCE = event->GetDCofThisEvent();
    if (DCE) {
      G4int nDC = (G4int)DCE->GetCapacity();
      for (G4int iDC = 0; iDC < nDC; ++iDC) {
        G4VDigiCollection* DC = DCE->GetDC(iDC);
        if (DC) {
          for(std::size_t iDigi = 0; iDigi < DC->GetSize(); ++iDigi) {
            fpCurrentDigi = DC->GetDigi(iDigi);
            if (fpCurrentDigi) sceneHandler.AddCompound(*fpCurrentDigi);
          }
        }
      }
    }
  }
}
