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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4hParametrisedLossModel
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
//
// Creation date: 20 July 2000
//
// Modifications:
// 20/07/2000  V.Ivanchenko First implementation
// 18/08/2000  V.Ivanchenko TRIM85 model is added
// 03/10/2000  V.Ivanchenko CodeWizard clean up
// 10/05/2001  V.Ivanchenko Clean up againist Linux compilation with -Wall
// 30/12/2003  V.Ivanchenko SRIM2003 model is added
// 07/05/2004  V.Ivanchenko Fix Graphite problem, add QAO model
//
// Class Description:
//
// Low energy protons/ions electronic stopping power parametrisation
//
// Class Description: End
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hParametrisedLossModel.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4hZiegler1985p.hh"
#include "G4hICRU49p.hh"
#include "G4hICRU49He.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ElementVector.hh"
#include "G4Material.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hParametrisedLossModel::G4hParametrisedLossModel(const G4String& name)
  :G4VLowEnergyModel(name), modelName(name)
{
  InitializeMe();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hParametrisedLossModel::InitializeMe()
{
  expStopPower125 = 0.;

  theZieglerFactor = eV*cm2*1.0e-15 ;

  // Registration of parametrisation models
  G4String blank  = G4String(" ") ;
  G4String ir49p  = G4String("ICRU_R49p") ;
  G4String ir49He = G4String("ICRU_R49He") ;
  G4String zi85p  = G4String("Ziegler1985p") ;
  if(zi85p == modelName) {
      eStopingPowerTable = new G4hZiegler1985p();
      highEnergyLimit = 100.0*MeV;
      lowEnergyLimit  = 1.0*keV;

  } else if(ir49p == modelName || blank == modelName) {
      eStopingPowerTable = new G4hICRU49p();
      highEnergyLimit = 2.0*MeV;
      lowEnergyLimit  = 1.0*keV;

  } else if(ir49He == modelName) {
      eStopingPowerTable = new G4hICRU49He();
      highEnergyLimit = 10.0*MeV/4.0;
      lowEnergyLimit  = 1.0*keV/4.0;

  } else {
      eStopingPowerTable = new G4hICRU49p();
      highEnergyLimit = 2.0*MeV;
      lowEnergyLimit  = 1.0*keV;
      G4cout << "G4hParametrisedLossModel Warning: <" << modelName
             << "> is unknown - default <"
             << ir49p << ">" << " is used for Electronic Stopping"
             << G4endl;
      modelName = ir49p;
  }
      /*
      G4cout << "G4hParametrisedLossModel: the model <"
             << modelName << ">" << " is used for Electronic Stopping"
             << G4endl;
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hParametrisedLossModel::~G4hParametrisedLossModel()
{
  delete eStopingPowerTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hParametrisedLossModel::TheValue(const G4DynamicParticle* particle,
					    const G4Material* material)
{
  G4double scaledEnergy = (particle->GetKineticEnergy())
                        * proton_mass_c2/(particle->GetMass());
  G4double factor = theZieglerFactor;
  if (scaledEnergy < lowEnergyLimit) {
    if (modelName != "QAO") factor *= std::sqrt(scaledEnergy/lowEnergyLimit);
    scaledEnergy = lowEnergyLimit;
  }
  G4double eloss = StoppingPower(material,scaledEnergy) * factor;

  return eloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hParametrisedLossModel::TheValue(const G4ParticleDefinition* aParticle,
					    const G4Material* material,
					          G4double kineticEnergy)
{
  G4double scaledEnergy = kineticEnergy
                        * proton_mass_c2/(aParticle->GetPDGMass());

  G4double factor = theZieglerFactor;
  if (scaledEnergy < lowEnergyLimit) {
    if (modelName != "QAO") factor *= std::sqrt(scaledEnergy/lowEnergyLimit);
    scaledEnergy = lowEnergyLimit;
  }
  G4double eloss = StoppingPower(material,scaledEnergy) * factor;

  return eloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hParametrisedLossModel::LowEnergyLimit(const G4ParticleDefinition* ,
						  const G4Material*) const
{
  return lowEnergyLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hParametrisedLossModel::HighEnergyLimit(const G4ParticleDefinition* ,
						   const G4Material*) const
{
  return highEnergyLimit;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hParametrisedLossModel::LowEnergyLimit(const G4ParticleDefinition* ) const
{
  return lowEnergyLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hParametrisedLossModel::HighEnergyLimit(const G4ParticleDefinition* ) const
{
  return highEnergyLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4hParametrisedLossModel::IsInCharge(const G4DynamicParticle* ,
					    const G4Material*) const
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4hParametrisedLossModel::IsInCharge(const G4ParticleDefinition* ,
					    const G4Material*) const
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hParametrisedLossModel::StoppingPower(const G4Material* material,
						       G4double kineticEnergy)
{
  G4double eloss = 0.0;

  const std::size_t numberOfElements = material->GetNumberOfElements() ;
  const G4double* theAtomicNumDensityVector =
    material->GetAtomicNumDensityVector() ;


  // compound material with parametrisation
  if( (eStopingPowerTable->HasMaterial(material)) ) {

    eloss = eStopingPowerTable->StoppingPower(material, kineticEnergy);
    if ("QAO" != modelName) {
      eloss *=  material->GetTotNbOfAtomsPerVolume();
      if(1 < numberOfElements) {
        G4int nAtoms = 0;

        const G4int* theAtomsVector = material->GetAtomsVector();
        for (std::size_t iel=0; iel<numberOfElements; ++iel) {
          nAtoms += theAtomsVector[iel];
        }
        eloss /= nAtoms;
      }
    }

  // pure material
  } else if(1 == numberOfElements) {

    G4double z = material->GetZ();
    eloss = (eStopingPowerTable->ElectronicStoppingPower(z, kineticEnergy))
                               * (material->GetTotNbOfAtomsPerVolume()) ;

  // Experimental data exist only for kinetic energy 125 keV
  } else if( MolecIsInZiegler1988(material)) {

  // Cycle over elements - calculation based on Bragg's rule
    G4double eloss125 = 0.0 ;
    const G4ElementVector* theElementVector =
                           material->GetElementVector() ;


    //  loop for the elements in the material
    for (std::size_t i=0; i<numberOfElements; ++i) {
      const G4Element* element = (*theElementVector)[i] ;
      G4double z = element->GetZ() ;
      eloss    +=(eStopingPowerTable->ElectronicStoppingPower(z,kineticEnergy))
                                    * theAtomicNumDensityVector[i] ;
      eloss125 +=(eStopingPowerTable->ElectronicStoppingPower(z,125.0*keV))
                                    * theAtomicNumDensityVector[i] ;
    }

    // Chemical factor is taken into account
    eloss *= ChemicalFactor(kineticEnergy, eloss125) ;

  // Brugg's rule calculation
  } else {
    const G4ElementVector* theElementVector =
                           material->GetElementVector() ;

    //  loop for the elements in the material
    for (std::size_t i=0; i<numberOfElements; ++i)
    {
      const G4Element* element = (*theElementVector)[i] ;
      G4double z = element->GetZ() ;
      eloss   += (eStopingPowerTable->ElectronicStoppingPower(z,kineticEnergy))
                                   * theAtomicNumDensityVector[i];
    }
  }
  return eloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4hParametrisedLossModel::MolecIsInZiegler1988(
                                 const G4Material* material)
{
  // The list of molecules from
  // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
  // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.

  G4String myFormula = G4String(" ") ;
  const G4String chFormula = material->GetChemicalFormula() ;
  if (myFormula == chFormula ) return false ;

  //  There are no evidence for difference of stopping power depended on
  //  phase of the compound except for water. The stopping power of the
  //  water in gas phase can be predicted using Bragg's rule.
  //
  //  No chemical factor for water-gas

  myFormula = G4String("H_2O") ;
  const G4State theState = material->GetState() ;
  if( theState == kStateGas && myFormula == chFormula) return false ;

  const std::size_t numberOfMolecula = 53 ;

  // The coffecient from Table.4 of Ziegler & Manoyan
  static const G4double HeEff = 2.8735 ;

  static const G4String name[numberOfMolecula] = {
    "H_2O",     "C_2H_4O",    "C_3H_6O",  "C_2H_2",             "C_H_3OH",
    "C_2H_5OH",  "C_3H_7OH",   "C_3H_4",   "NH_3",               "C_14H_10",
    "C_6H_6",    "C_4H_10",    "C_4H_6",   "C_4H_8O",            "CCl_4",
    "CF_4",      "C_6H_8",     "C_6H_12",  "C_6H_10O",           "C_6H_10",
    "C_8H_16",   "C_5H_10",    "C_5H_8",   "C_3H_6-Cyclopropane","C_2H_4F_2",
    "C_2H_2F_2", "C_4H_8O_2",  "C_2H_6",   "C_2F_6",             "C_2H_6O",
    "C_3H_6O",   "C_4H_10O",   "C_2H_4",   "C_2H_4O",            "C_2H_4S",
    "SH_2",      "CH_4",       "CCLF_3",   "CCl_2F_2",           "CHCl_2F",
    "(CH_3)_2S", "N_2O",       "C_5H_10O", "C_8H_6",          "(CH_2)_N",
    "(C_3H_6)_N","(C_8H_8)_N", "C_3H_8",   "C_3H_6-Propylene",   "C_3H_6O",
    "C_3H_6S",   "C_4H_4S",    "C_7H_8"
  };

  static const G4double expStopping[numberOfMolecula] = {
     66.1,  190.4, 258.7,  42.2, 141.5,
    210.9,  279.6, 198.8,  31.0, 267.5,
    122.8,  311.4, 260.3, 328.9, 391.3,
    206.6,  374.0, 422.0, 432.0, 398.0,
    554.0,  353.0, 326.0,  74.6, 220.5,
    197.4,  362.0, 170.0, 330.5, 211.3,
    262.3,  349.6,  51.3, 187.0, 236.9,
    121.9,   35.8, 247.0, 292.6, 268.0,
    262.3,   49.0, 398.9, 444.0,  22.91,
     68.0,  155.0,  84.0,  74.2, 254.7,
    306.8,  324.4, 420.0
  } ;

  static const G4double expCharge[numberOfMolecula] = {
    HeEff, HeEff, HeEff,   1.0, HeEff,
    HeEff, HeEff, HeEff,   1.0,   1.0,
      1.0, HeEff, HeEff, HeEff, HeEff,
    HeEff, HeEff, HeEff, HeEff, HeEff,
    HeEff, HeEff, HeEff,   1.0, HeEff,
    HeEff, HeEff, HeEff, HeEff, HeEff,
    HeEff, HeEff,   1.0, HeEff, HeEff,
    HeEff,   1.0, HeEff, HeEff, HeEff,
    HeEff,   1.0, HeEff, HeEff,   1.0,
      1.0,   1.0,   1.0,   1.0, HeEff,
    HeEff, HeEff, HeEff
  } ;

  static const G4double numberOfAtomsPerMolecula[numberOfMolecula] = {
    3.0,  7.0, 10.0,  4.0,  6.0,
    9.0, 12.0,  7.0,  4.0, 24.0,
    12.0, 14.0, 10.0, 13.0,  5.0,
    5.0, 14.0, 18.0, 17.0, 17.0,
    24.0, 15.0, 13.0,  9.0,  8.0,
    6.0, 14.0,  8.0,  8.0,  9.0,
    10.0, 15.0,  6.0,  7.0,  7.0,
    3.0,  5.0,  5.0,  5.0,  5.0,
    9.0,  3.0, 16.0, 14.0,  3.0,
    9.0, 16.0, 11.0,  9.0, 10.0,
    10.0,  9.0, 15.0
  } ;

  // Search for the compaund in the table
  for (std::size_t i=0; i<numberOfMolecula; ++i)
    {
      if(chFormula == name[i]) {
        G4double exp125 = expStopping[i] *
	                  (material->GetTotNbOfAtomsPerVolume()) /
	                  (expCharge[i] * numberOfAtomsPerMolecula[i]) ;
        SetExpStopPower125(exp125) ;
        return true ;
      }
    }

  return false ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hParametrisedLossModel::ChemicalFactor(
                            G4double kineticEnergy, G4double eloss125) const
{
  // Approximation of Chemical Factor according to
  // J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
  // Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.

  G4double gamma    = 1.0 + kineticEnergy/proton_mass_c2 ;
  G4double gamma25  = 1.0 + 25.0*keV /proton_mass_c2 ;
  G4double gamma125 = 1.0 + 125.0*keV/proton_mass_c2 ;
  G4double beta     = std::sqrt(1.0 - 1.0/(gamma*gamma)) ;
  G4double beta25   = std::sqrt(1.0 - 1.0/(gamma25*gamma25)) ;
  G4double beta125  = std::sqrt(1.0 - 1.0/(gamma125*gamma125)) ;

  G4double factor = 1.0 + (expStopPower125/eloss125 - 1.0) *
                   (1.0 + G4Exp( 1.48 * ( beta125/beta25 - 7.0 ) ) ) /
                   (1.0 + G4Exp( 1.48 * ( beta/beta25    - 7.0 ) ) ) ;

  return factor ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
