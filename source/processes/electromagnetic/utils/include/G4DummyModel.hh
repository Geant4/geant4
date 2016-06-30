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
// $Id: G4DummyModel.hh 96139 2016-03-17 14:10:31Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4DummyModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.03.2006
//
// Modifications:
//
// Class Description:
//
// EM model doing nothing

// -------------------------------------------------------------------
//

#ifndef G4DummyModel_h
#define G4DummyModel_h 1

#include "globals.hh"
#include "G4VMscModel.hh"

class G4DummyModel :  public G4VMscModel
{

public:

  explicit G4DummyModel(const G4String& nam = "DummyModel");

  virtual ~G4DummyModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) final;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double tmax) final;

private:

  //  hide assignment operator
  G4DummyModel & operator=(const  G4DummyModel &right) = delete;
  G4DummyModel(const  G4DummyModel&) = delete;
};

#endif

