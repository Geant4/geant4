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
// $Id: G4WentzelVIRelModel.hh 104307 2017-05-24 09:01:45Z gcosmo $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4WentzelVIRelModel
//
// Author:        V.Ivanchenko 
//
// Creation date: 08.06.2012 from G4WentzelVIModel
//
// Modifications:
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// G.Wentzel, Z. Phys. 40 (1927) 590.
// H.W.Lewis, Phys Rev 78 (1950) 526.
// J.M. Fernandez-Varea et al., NIM B73 (1993) 447.
// L.Urban, CERN-OPEN-2006-077.

// -------------------------------------------------------------------
//

#ifndef G4WentzelVIRelModel_h
#define G4WentzelVIRelModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WentzelVIModel.hh"
#include "G4Threading.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4NistManager;

class G4WentzelVIRelModel : public G4WentzelVIModel
{

public:

  explicit G4WentzelVIRelModel();

  virtual ~G4WentzelVIRelModel();

  
  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) override;

  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
					      G4double KineticEnergy,
					      G4double AtomicNumber,
					      G4double AtomicWeight=0., 
					      G4double cut = DBL_MAX,
					      G4double emax= DBL_MAX) override;

  virtual void DefineMaterial(const G4MaterialCutsCouple* cup);

private:

  void ComputeEffectiveMass();

  //  hide assignment operator
  G4WentzelVIRelModel & operator=(const  G4WentzelVIRelModel &right) = delete;
  G4WentzelVIRelModel(const  G4WentzelVIRelModel&) = delete;

  static std::vector<G4double> effMass;
  G4NistManager* fNistManager;

#ifdef G4MULTITHREADED
  static G4Mutex WentzelVIRelModelMutex;
#endif
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

