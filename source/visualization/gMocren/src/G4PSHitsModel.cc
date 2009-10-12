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
// $Id: G4PSHitsModel.cc,v 1.2 2009-10-12 11:25:15 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Created:  Mar. 31, 2009  Akinori Kimura
// Model which draws the primtive scorer hits.
//

#include "G4PSHitsModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4GMocrenFileSceneHandler.hh"

G4PSHitsModel::~G4PSHitsModel () {}

G4PSHitsModel::G4PSHitsModel ():
  kpCurrentHit(0)
{
  fGlobalTag = "G4PSHitsModel for all hits.";
  fGlobalDescription = fGlobalTag;
}

void G4PSHitsModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
  const G4Event* event = fpMP->GetEvent();
  G4GMocrenFileSceneHandler * gHandler = dynamic_cast<G4GMocrenFileSceneHandler*>(&sceneHandler);
  if(gHandler) {
    if (event) {
      G4HCofThisEvent * HCE = event->GetHCofThisEvent();
      if (HCE) {
	G4int nHC = HCE -> GetCapacity ();
	for (int iHC = 0; iHC < nHC; iHC++) {
	  G4THitsMap<G4double> * hits = dynamic_cast<G4THitsMap<G4double> *>(HCE -> GetHC (iHC));
	  if(hits) gHandler->AddCompound (*hits);
	  
	}
      }
    }
  }
}
