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
// Created on 2016/04/08
//
// Authors: D. Sakata, S. Incerti
//
// This class perform electric excitation for electron transportation,
// based on Dirac B-Spline R-Matrix Model and scaled experimental data.
// See following reference paper 
// Phys.Rev.A77,062711(2008) and Phys.Rev.A78,042713(2008)	

#ifndef G4DNADiracRMatrixExcitationModel_h
#define G4DNADiracRMatrixExcitationModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VAtomDeexcitation.hh"

#include "G4LogLogInterpolation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4NistManager.hh"

#include "G4DNACrossSectionDataSet.hh"

class G4DNADiracRMatrixExcitationModel: public G4VEmModel
{

public:

  G4DNADiracRMatrixExcitationModel(const G4ParticleDefinition* p = 0,
                        const G4String& nam = "DNADiracRMatrixExcitationModel");

  virtual ~G4DNADiracRMatrixExcitationModel();

  virtual void Initialise(const G4ParticleDefinition*,
                          const G4DataVector& = *(new G4DataVector()));

  virtual G4double CrossSectionPerVolume(const G4Material* material,
                                         const G4ParticleDefinition* p,
                                         G4double ekin,
                                         G4double emin,
                                         G4double emax);

  virtual G4double GetExtendedTotalCrossSection  (const G4Material* material,
                                          const G4ParticleDefinition*,
                                          G4double kineticEnergy);
  
  virtual G4double GetExtendedPartialCrossSection(const G4Material* material,
                                          G4int level,
                                          const G4ParticleDefinition*,
                                          G4double kineticEnergy);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy);

  inline  void SelectStationary(G4bool input);

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:

  const G4double paramFuncTCS_5dto6s1[3]={-3e-50      , 9.46358e-16,  1.4237  }; // y = [0]+[1]/pow(x-[2],2)
  const G4double paramFuncTCS_5dto6s2[3]={-3e-50      , 4.24498e-15, -0.674543}; // y = [0]+[1]/pow(x-[2],2)
  const G4double paramFuncTCS_6sto6p1[3]={ 1.50018e-26, 2.459e-15  ,-40.8088  }; // y = [0]+[1]*log(x-[2])/(x-[2])
  const G4double paramFuncTCS_6sto6p2[3]={ 1.26684e-25, 3.97221e-15,-55.6954  }; // y = [0]+[1]*log(x-[2])/(x-[2])
  const G4int    ShellEnumAu       [4]={19    , 20   ,21    , 21   }; 
  //      5d3/2  ,6s1/2  ,6s1/2   //from EADL
  const G4double BindingEnergyAu   [4]={12.16 ,10.46 , 8.3  ,  8.3 }; 
  // [eV] 5d3/2  ,6s1/2  ,6s1/2   //from EADL
  const G4double ExcitationEnergyAu[4]={ 2.66 , 1.14 , 4.63 ,  5.11}; 
  // [eV] 5dto6s1,6sto6p1,6sto6p2

  G4double fLowEnergyLimit=0.;
  G4double fExperimentalEnergyLimit=0.;
  G4double fHighEnergyLimit=0.;

  G4bool   isInitialised=false;
  G4bool   statCode=false;
  G4int    verboseLevel=0;

  G4String                     fTableFile="";
  G4DNACrossSectionDataSet*    fTableData=nullptr;
  const std::vector<G4double>* fpMaterialDensity=nullptr;
  const  G4ParticleDefinition* fParticleDefinition=nullptr;
  G4VAtomDeexcitation*         fAtomDeexcitation=nullptr;

  G4int RandomSelect(const G4Material* material,
                     const G4ParticleDefinition*,
                     G4double kineticEnergy);
  
   
  G4DNADiracRMatrixExcitationModel & operator
                           =(const  G4DNADiracRMatrixExcitationModel &right);
  G4DNADiracRMatrixExcitationModel(const  G4DNADiracRMatrixExcitationModel&);

};

inline void G4DNADiracRMatrixExcitationModel::SelectStationary(G4bool input)
{
  statCode = input;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
