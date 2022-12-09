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
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelEventAction  ------
//           by R.Giannitrapani, F. Longo & G.Santin (13 nov 2000)
//
//- inclusion of Digits by F.Longo & R.Giannitrapani (24 oct 2001)
//  20.11.01 G.Santin: new analysis management, modified according to GammaRayTelAnalysis
// ************************************************************

#ifndef GammaRayTelEventAction_h
#define GammaRayTelEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "GammaRayTelRunAction.hh"

class GammaRayTelEventAction: public G4UserEventAction {
public:
	explicit GammaRayTelEventAction(GammaRayTelRunAction *runAction = nullptr);

	~GammaRayTelEventAction() override;

	void BeginOfEventAction(const G4Event *event) override;

	void EndOfEventAction(const G4Event *event) override;

	void SetDrawFlag(G4String value) {
		drawFlag = value;
	};

private:
    G4int anticoincidenceCollectionID{-1};

    G4int calorimeterCollectionID{-1};

	G4int trackerCollectionID{-1};

	G4String drawFlag{"all"};

	GammaRayTelRunAction *theRunAction;
};
#endif
