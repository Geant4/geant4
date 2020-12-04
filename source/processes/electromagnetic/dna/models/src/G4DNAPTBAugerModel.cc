//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                            *
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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)
//

#include "G4DNAPTBAugerModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

#include "G4Material.hh"

using namespace std;

G4DNAPTBAugerModel::G4DNAPTBAugerModel(const G4String& modelAugerName): modelName(modelAugerName)
{
    // To inform the user that the Auger model is enabled
    G4cout << modelName <<" is constructed" << G4endl;
}

G4DNAPTBAugerModel::~G4DNAPTBAugerModel()
{
    if( verboseLevel>0 ) G4cout << modelName <<" is deleted" << G4endl;
}

void G4DNAPTBAugerModel::Initialise()
{
    verboseLevel = 0;

    if( verboseLevel>0 )
    {
        G4cout << "PTB Auger model is initialised " << G4endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBAugerModel::ComputeAugerEffect(std::vector<G4DynamicParticle*>* fvect, const G4String& materialNameIni, G4double bindingEnergy)
{
    // Rename material if modified NIST material
    // This is needed when material is obtained from G4MaterialCutsCouple
    G4String materialName = materialNameIni;
    if(materialName.find("_MODIFIED")){
        materialName = materialName.substr(0,materialName.size()-9);
    }

    // check if there is a k-shell ionisation and find the ionised atom
    G4int atomId(0);

    atomId = DetermineIonisedAtom(atomId, materialName, bindingEnergy);

    if(atomId!=0)
    {
        G4double kineticEnergy = CalculAugerEnergyFor(atomId);

        if(kineticEnergy<0)
        {
            G4cerr<<"**************************"<<G4endl;
            G4cerr<<"FatalError. Auger kineticEnergy: "<<kineticEnergy<<G4endl;
            exit(EXIT_FAILURE);
        }

        if(atomId==1 || atomId==2 || atomId==3)
        {
            GenerateAugerWithRandomDirection(fvect, kineticEnergy);
        }
        else if(atomId==4)
        {
            GenerateAugerWithRandomDirection(fvect, kineticEnergy);
            GenerateAugerWithRandomDirection(fvect, kineticEnergy);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4DNAPTBAugerModel::DetermineIonisedAtom(G4int atomId, const G4String& materialName, G4double bindingEnergy)
{
    if(materialName=="THF" || materialName=="backbone_THF"){
        if(bindingEnergy==305.07){
            atomId=1; //"carbon";
        }
        else if(bindingEnergy==557.94){
            atomId=2; //"oxygen";
        }
    }
    else if(materialName=="PY" || materialName=="PU"
            || materialName=="cytosine_PY" || materialName=="thymine_PY"
            || materialName=="adenine_PU" || materialName=="guanine_PU"
            )
    {
        if(bindingEnergy==307.52){
            atomId=1; //"carbon";
        }
        else if(bindingEnergy==423.44){
            atomId=4; //"nitrogen";
        }
    }
    else if(materialName=="TMP"|| materialName=="backbone_TMP"){
        if(bindingEnergy==209.59 || bindingEnergy==152.4)
            atomId=3; //"carbonTMP";
    }

    return atomId;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBAugerModel::CalculAugerEnergyFor(G4int atomId)
{
    G4double kineticEnergy;

    if(atomId==2) // oxygen
    {
        kineticEnergy = 495*eV;
    }
    else
    {
        G4double f1, f2, f3, g1, g2, Y;

        Y = G4UniformRand();

        if(atomId == 1){ // carbon
            f1 = -7.331e-2;
            f2 = -3.306e-5;
            f3 = 2.433e0;
            g1 = 4.838e-1;
            g2 = 3.886e0;
        }
        else if(atomId == 4){ // nitrogen
            f1 = -7.518e-2;
            f2 = 1.178e-4;
            f3 = 2.600e0;
            g1 = 4.639e-1;
            g2 = 3.770e0;
        }
        else// if(atomId == 3) // carbon_TMP
        {
            f1 = -5.700e-2;
            f2 = 1.200e-4;
            f3 = 2.425e0;
            g1 = 5.200e-1;
            g2 = 2.560e0;
        }

        kineticEnergy = pow(10, f1*pow( abs( log10(Y) ) , g1) + f2*pow( abs( log10(Y) ) , g2) + f3 )*eV;
    }

    return kineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBAugerModel::SetCutForAugerElectrons(G4double cut)
{
  minElectronEnergy = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBAugerModel::GenerateAugerWithRandomDirection(std::vector<G4DynamicParticle*>* fvect, G4double kineticEnergy)
{
      // Isotropic angular distribution for the outcoming e-
      G4double newcosTh = 1.-2.*G4UniformRand();
      G4double  newsinTh = std::sqrt(1.-newcosTh*newcosTh);
      G4double newPhi = twopi*G4UniformRand();
      
      G4double xDir =  newsinTh*std::sin(newPhi);
      G4double yDir = newsinTh*std::cos(newPhi);
      G4double zDir = newcosTh;
      
      G4ThreeVector ElectronDirection(xDir,yDir,zDir);

      // generation of new particle
      G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(), ElectronDirection, kineticEnergy) ;
      fvect->push_back(dp);
}
