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
// Contact authors: S. Meylan, C. Villagrasa
//
// email: sylvain.meylan@symalgo-tech.com, carmen.villagrasa@irsn.fr

#ifndef G4DNADUMMYMODEL_HH
#define G4DNADUMMYMODEL_HH

#include "G4VDNAModel.hh"
#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4ParticleChangeForGamma.hh"

class G4DNADummyModel : public G4VDNAModel
{
public:
    G4DNADummyModel(const G4String& applyToMaterial,
                    const G4ParticleDefinition* p,
                    const G4String& nam,
                    G4VEmModel* emModel);
    ~G4DNADummyModel();

    virtual void Initialise(const G4ParticleDefinition* particle, const G4DataVector& = *(new G4DataVector()), G4ParticleChangeForGamma* changeForGamme=nullptr);

    virtual G4double CrossSectionPerVolume(const G4Material* material,
                                           const G4String& materialName,
                                           const G4ParticleDefinition* p,
                                           G4double ekin,
                                           G4double emin,
                                           G4double emax);

    virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                   const G4MaterialCutsCouple*,
                                   const G4String& materialName,
                                   const G4DynamicParticle*,
                                   G4ParticleChangeForGamma *particleChangeForGamma,
                                   G4double tmin,
                                   G4double tmax);

    const G4VEmModel* GetEmModel() const {return fpEmModel;}
    G4VEmModel* GetEmModel() {return fpEmModel;}

private:
    G4VEmModel* fpEmModel;
    const G4ParticleDefinition* fpParticleDef;
    const std::vector<double>* fMaterialMolPerVol;

    G4double GetNumMoleculePerVolumeUnitForMaterial(const G4Material *mat);
};

#endif // G4DNADUMMYMODEL_HH
