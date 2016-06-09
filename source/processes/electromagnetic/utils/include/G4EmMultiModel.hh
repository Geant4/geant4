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
// $Id: G4EmMultiModel.hh,v 1.6 2007/05/22 17:31:57 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmMultiModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.05.2004
//
// Modifications:
// 15-04-05 optimize internal interface (V.Ivanchenko)
//
// Class Description:
//
// Energy loss model using several G4VEmModels

// -------------------------------------------------------------------
//

#ifndef G4EmMultiModel_h
#define G4EmMultiModel_h 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include <vector>

class G4Region;
class G4PhysicsTable;
class G4DynamicParticle;

class G4EmMultiModel :  public G4VEmModel
{

public:

  G4EmMultiModel(const G4String& nam = "MultiModel");

  virtual ~G4EmMultiModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*);


  G4double ComputeDEDX(const G4MaterialCutsCouple*,
		       const G4ParticleDefinition*,
		       G4double kineticEnergy,
		       G4double cutEnergy);

  G4double CrossSection(const G4MaterialCutsCouple*,
			const G4ParticleDefinition*,
			G4double kineticEnergy,
			G4double cutEnergy,
			G4double maxEnergy);

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double tmax);

  void DefineForRegion(const G4Region*);

  void AddModel(G4VEmModel*, G4double tmin, G4double tmax);

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
    				    G4double kineticEnergy);

private: 

  //  hide assignment operator
  G4EmMultiModel & operator=(const  G4EmMultiModel &right);
  G4EmMultiModel(const  G4EmMultiModel&);

  G4int                         nModels;
  std::vector<G4VEmModel*>      model;
  G4DataVector                  tsecmin;
  G4DataVector                  cross_section;

};

#endif

