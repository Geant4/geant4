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
// Geant4 header : G4LowEGammaNuclearModel
// Created:  15 May 2019
// Author  V.Ivanchenko
//
// Class Description
// Sampling of gamma nuclear interaction at low energy 
// Class Description - End
//

#ifndef G4LowEGammaNuclearModel_h
#define G4LowEGammaNuclearModel_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4LorentzVector.hh"

class G4PreCompoundModel;

class G4LowEGammaNuclearModel : public G4HadronicInteraction
{
public:

  explicit G4LowEGammaNuclearModel();

  ~G4LowEGammaNuclearModel() override;
 
  G4HadFinalState* ApplyYourself(const G4HadProjectile & aTrack, 
				 G4Nucleus & targetNucleus) final;

  void InitialiseModel() final;

private:

  G4LowEGammaNuclearModel & operator=
  (const G4LowEGammaNuclearModel &right);
  G4LowEGammaNuclearModel(const G4LowEGammaNuclearModel&);

  G4PreCompoundModel* fPreco;
  G4LorentzVector lab4mom;

  G4int secID;  // Creator model ID for the secondaries created by this model
};

#endif
