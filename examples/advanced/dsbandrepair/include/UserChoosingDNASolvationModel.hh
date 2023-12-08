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
/// \file UserChoosingDNASolvationModel.hh
/// \brief Definition of the UserChoosingDNASolvationModel class

#ifndef UserChoosingDNASolvationModel_h
#define UserChoosingDNASolvationModel_h 1
#include "G4DNAOneStepThermalizationModel.hh"

template<typename MODEL=DNA::Penetration::Meesungnoen2002>
class UserTDNAOneStepThermalizationModel: public G4TDNAOneStepThermalizationModel<MODEL>
{
public:
    typedef MODEL Model;
    UserTDNAOneStepThermalizationModel(const G4ParticleDefinition* p = 0,
                                        const G4String& nam ="DNAOneStepThermalizationModel");
    ~UserTDNAOneStepThermalizationModel() override = default;

    void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy) override;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UserChoosingDNASolvationModel
{
public:
    static G4VEmModel* UserCreate(const G4String& penetrationModel);
  
    static G4VEmModel* UserGetMacroDefinedModel();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif