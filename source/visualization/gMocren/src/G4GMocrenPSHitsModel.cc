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
// $Id: G4GMocrenPSHitsModel.cc,v 1.1 2009-10-12 10:24:23 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Created:  Mar. 31, 2009  Akinori Kimura  
//

#include "G4GMocrenPSHitsModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"

G4GMocrenPSHitsModel::~G4GMocrenPSHitsModel () {}

G4GMocrenPSHitsModel::G4GMocrenPSHitsModel ():
  kpCurrentHit(0)
{
  fGlobalTag = "G4GMocrenPSHitsModel for all hits.";
  fGlobalDescription = fGlobalTag;
}

#include "G4GMocrenFileSceneHandler.hh"
const bool DEBUG = false;//true;
void G4GMocrenPSHitsModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
  if(DEBUG) G4cout << "G4GMocrenPSHitsModel::DescribeYourselfTo() "
		   << (unsigned long)fpMP << G4endl;
  G4GMocrenFileSceneHandler * ghandler = dynamic_cast<G4GMocrenFileSceneHandler*>(&sceneHandler);
  if(ghandler) {
    const G4Event* event = fpMP->GetEvent();
    if (event) {
      G4HCofThisEvent * HCE = event->GetHCofThisEvent();
      if (HCE) {
	G4int nHC = HCE -> GetCapacity ();
	if(DEBUG) G4cout << "            # of HC: " << nHC << G4endl;
	for (int iHC = 0; iHC < nHC; iHC++) {
	  G4THitsMap<G4double> * hits = dynamic_cast<G4THitsMap<G4double> *>(HCE -> GetHC (iHC));
	  if(hits) ghandler->AddCompound (*hits);

	}
      }
    }
  }
}
