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
// $Id: G4HitsModel.cc,v 1.14 2006-05-22 07:26:24 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th August 1998.
// Model which knows how to draw GEANT4 hits.

#include "G4HitsModel.hh"

#include "G4VGraphicsScene.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

G4HitsModel::~G4HitsModel () {}

G4HitsModel::G4HitsModel () {
  fGlobalTag = "G4HitsModel for all hits.";
  fGlobalDescription = fGlobalTag;
}

void G4HitsModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler) {
  const G4Event* event =
    G4EventManager::GetEventManager()->GetConstCurrentEvent();
  if (event) {
    G4HCofThisEvent* HCE = event -> GetHCofThisEvent ();
    if (HCE) {
      G4int nHC = HCE -> GetCapacity ();
      for (int iHC = 0; iHC < nHC; iHC++) {
	G4VHitsCollection* HC = HCE -> GetHC (iHC);
	if (HC) {
	  for(size_t iHit = 0; iHit < HC->GetSize(); ++iHit) {
	    G4VHit* hit = HC -> GetHit (iHit);
	    sceneHandler.AddCompound (*hit);
	  }
	}
      }
    }
  }
}
