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
///////////////////////////////////////////////////////
// File name:     G4hInverseIonisation
//
// Author:        Laurent Desorgher
//
// Creation date: 15.02.2009
//
///////////////////////////////////////////////////////
#include "G4hInverseIonisation.hh"
#include "G4VEmAdjointModel.hh"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4hInverseIonisation::G4hInverseIonisation(G4bool whichScatCase,G4String process_name,G4AdjointhIonisationModel* aEmAdjointModel):
				G4VAdjointReverseReaction(process_name,whichScatCase)
{theAdjointEMModel = aEmAdjointModel;
 theAdjointEMModel->SetSecondPartOfSameType(false);
 SetIntegralMode(true); 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4hInverseIonisation::~G4hInverseIonisation(){
}
