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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VEmFluctuationModel
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
// 28-12-02 add method Dispersion (V.Ivanchenko)
// 07-02-03 change signature (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
//
//
// Class Description: 
//
// Abstract class for interface to simualtion of energy loss fluctuations

// -------------------------------------------------------------------
//

#ifndef G4VEmFluctuationModel_h
#define G4VEmFluctuationModel_h 1


#include "globals.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4Material;

class G4VEmFluctuationModel 
{

public:

  G4VEmFluctuationModel(const G4String& nam): name(nam) {};

  virtual ~G4VEmFluctuationModel() {};

  virtual G4double SampleFluctuations(const G4Material*,
                                  const G4DynamicParticle*,
				        G4double& tmax,
                                        G4double& length,
                                        G4double& meanLoss) = 0;

  virtual G4double Dispersion(const G4Material*,
                              const G4DynamicParticle*,
				        G4double& tmax,
                                        G4double& length) = 0;

  virtual void Initialise(const G4ParticleDefinition*) = 0;
  
  G4String GetName() const {return name;};

protected:

private:

  // hide assignment operator
  G4VEmFluctuationModel & operator=(const  G4VEmFluctuationModel &right);
  G4VEmFluctuationModel(const  G4VEmFluctuationModel&);

  const G4String   name;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

