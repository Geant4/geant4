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
// Created on 2016/05/02
//
// Authors: D Sakata, S. Incerti
//
// This class perform electric excitation for electron transportation in gold,
// based on Dirac B-Spline R-Matrix method with scaled experimental data 
// for low energy.
// See following reference paper 
// Phys.Rev.A77,062711(2008) and Phys.Rev.A78,042713(2008)	

#include "G4DNADiracRMatrixExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4Gamma.hh"
#include "G4RandomDirection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNADiracRMatrixExcitationModel::G4DNADiracRMatrixExcitationModel
(const G4ParticleDefinition*,const G4String& nam) :
    G4VEmModel(nam), isInitialised(false), fTableData(0)
{
  fpMaterialDensity       = 0;
  fHighEnergyLimit        = 0;
  fExperimentalEnergyLimit= 0;
  fLowEnergyLimit         = 0;
  fParticleDefinition     = 0;

  verboseLevel = 0;

  if (verboseLevel > 0)
  {
    G4cout << "Dirac R-matrix excitation model is constructed " << G4endl;
  }
  
  fParticleChangeForGamma = 0;
  statCode                = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNADiracRMatrixExcitationModel::~G4DNADiracRMatrixExcitationModel()
{
  if (fTableData) delete fTableData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNADiracRMatrixExcitationModel::Initialise
(const G4ParticleDefinition* particle,const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
  {
    G4cout << 
           "Calling G4DNADiracRMatrixExcitationModel::Initialise()" 
           << G4endl;
  }

  fParticleDefinition = particle;
  
  if(particle->GetParticleName() == "e-")
  {
    fTableFile = "dna/sigma_excitation_e_diracrmatrix_Z79";
    fLowEnergyLimit           =  10  *  eV;
    fExperimentalEnergyLimit  =  577.*  eV;
    fHighEnergyLimit          =  1.0 * GeV;
  }
  else
  { 
   G4Exception("G4DNADiracRMatrixExcitationModel::Initialise","em0001",
        FatalException,"Not defined for other particles than electrons.");
   return;
  }

  G4double scaleFactor = 1. * cm * cm;
  fTableData = new G4DNACrossSectionDataSet
              (new G4LogLogInterpolation,eV,scaleFactor );
  fTableData->LoadData(fTableFile);

  if( verboseLevel>0 )
  {
    G4cout << "Dirac R-matrix excitation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "<< HighEnergyLimit() / keV << " keV "
    << " for "<< particle->GetParticleName()
    << G4endl;
  }

  if (isInitialised){return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNADiracRMatrixExcitationModel::CrossSectionPerVolume
                         (const G4Material* material,
                          const G4ParticleDefinition* particleDefinition,
                          G4double ekin,
                          G4double,
                          G4double)
{
  if (verboseLevel > 3)
  {
    G4cout << 
        "Calling CrossSectionPerVolume() of G4DNADiracRMatrixExcitationModel"
           << G4endl;
  }

  G4double atomicNDensity = material->GetAtomicNumDensityVector()[0];

  // Protection: for single element
  if(material->GetNumberOfElements()>1) return 0.; 

  G4double z              = material->GetZ();

  // Protection: for Gold
  if(z!=79){return 0.;}

  G4double sigma=0.;

  if(atomicNDensity!= 0.0)
  {
    if (ekin >= fLowEnergyLimit && ekin < fExperimentalEnergyLimit)
    {
      sigma = fTableData->FindValue(ekin);
    }
    else if ((fExperimentalEnergyLimit <= ekin) && (ekin < fHighEnergyLimit))
    {    
      sigma = GetExtendedTotalCrossSection(material,particleDefinition,ekin);
    }

    if (verboseLevel > 2)
    {
      G4cout<<"__________________________________" << G4endl;
      G4cout<<"=== G4DNADiracRMatrixExcitationModel - XS INFO START"<<G4endl;
      G4cout<<"=== Kinetic energy (eV)=" << ekin/eV << " particle : " 
            <<particleDefinition->GetParticleName() << G4endl;
      G4cout<<"=== Cross section per atom for Z="<<z<<" is (cm^2)" 
            <<sigma/cm/cm << G4endl;
      G4cout<<"=== Cross section per atom for Z="<<z<<" is (cm^-1)=" 
            <<sigma*atomicNDensity/(1./cm) << G4endl;
      G4cout<<"=== G4DNADiracRMatrixExcitationModel - XS INFO END"<<G4endl;
    }
  } 

  return sigma*atomicNDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNADiracRMatrixExcitationModel::SampleSecondaries
                         (std::vector<G4DynamicParticle*>* /*fvect*/,
                          const G4MaterialCutsCouple* couple,
                          const G4DynamicParticle* aDynamicParticle,
                          G4double,G4double)
{

  if (verboseLevel > 3)
  {
    G4cout << 
          "Calling SampleSecondaries() of G4DNADiracRMatrixExcitationModel"
           << G4endl;
  }

  G4ParticleDefinition* particle = aDynamicParticle->GetDefinition();
  G4double k                     = aDynamicParticle->GetKineticEnergy();

  G4int    level                 = RandomSelect(couple->GetMaterial(),particle,
                                                k);
  G4double excitationEnergy      = ExcitationEnergyAu[level]*eV;
  G4double newEnergy             = k - excitationEnergy;

  if (newEnergy > 0)
  {
    //Energy Loss
    fParticleChangeForGamma->ProposeMomentumDirection 
                  (aDynamicParticle->GetMomentumDirection());
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
    if(!statCode) fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
    else          fParticleChangeForGamma->SetProposedKineticEnergy(k);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNADiracRMatrixExcitationModel::GetExtendedTotalCrossSection
                          (const G4Material* material,
                           const G4ParticleDefinition* particle,
                           G4double kineticEnergy)
{
  G4double value=0;
  
  size_t N=fTableData->NumberOfComponents();
  
  for(int i=0;i<(int)N;i++){
   value = value+GetExtendedPartialCrossSection(material,i,particle,
                                                kineticEnergy);
  }
  
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNADiracRMatrixExcitationModel::GetExtendedPartialCrossSection
                         (const G4Material*,
                          G4int level,
                          const G4ParticleDefinition* particle,
                          G4double kineticEnergy)
{
  G4double value=0; 
  
  if(particle->GetParticleName()=="e-"){
  
    if(level==0){
      // y = [0]+[1]/pow(x-2,2)
      value = paramFuncTCS_5dto6s1[0]+paramFuncTCS_5dto6s1[1]
              /std::pow(kineticEnergy/eV-paramFuncTCS_5dto6s1[2],2);
    }
    else if(level==1){
      // y = [0]+[1]/pow(x-2,2)
      value = paramFuncTCS_5dto6s2[0]+paramFuncTCS_5dto6s2[1]
              /std::pow(kineticEnergy/eV-paramFuncTCS_5dto6s2[2],2);
    }
    else if(level==2){
      // y = [0]+[1]*log(x-2)/(x-[2])
      value = paramFuncTCS_6sto6p1[0]+paramFuncTCS_6sto6p1[1]
              *G4Log(kineticEnergy/eV-paramFuncTCS_6sto6p1[2])
              /(kineticEnergy/eV-paramFuncTCS_6sto6p1[2]);
    }
    else if(level==3){
      // y = [0]+[1]*log(x-2)/(x-[2])
      value = paramFuncTCS_6sto6p2[0]+paramFuncTCS_6sto6p2[1]
              *G4Log(kineticEnergy/eV-paramFuncTCS_6sto6p2[2])
              /(kineticEnergy/eV-paramFuncTCS_6sto6p2[2]);
    }
  }

  return value*cm*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNADiracRMatrixExcitationModel::RandomSelect
                         (const G4Material* material,
                          const G4ParticleDefinition* particle,
                          G4double kineticEnergy)
{
  G4double value = 0.;

  std::size_t NOfComp = fTableData->NumberOfComponents();
 
  auto valuesBuffer = new G4double[NOfComp];
  
  const G4int n = (G4int)fTableData->NumberOfComponents();
  
  G4int i(n);

  while (i > 0)
  {
    --i;
    if 
    ((fLowEnergyLimit<=kineticEnergy)&&(kineticEnergy<fExperimentalEnergyLimit))
    {
      valuesBuffer[i] = fTableData->GetComponent(i)->FindValue(kineticEnergy);
    }
    else if 
    ((fExperimentalEnergyLimit<=kineticEnergy)&&(kineticEnergy<fHighEnergyLimit))
    {
      valuesBuffer[i] 
            = GetExtendedPartialCrossSection(material,i,particle,kineticEnergy);
    }
    value += valuesBuffer[i];
  }
  value *= G4UniformRand();
  i = n;
  while (i > 0)
  {
    --i;
    if (valuesBuffer[i] > value)
    {
      delete[] valuesBuffer;
      return i;
    }
    value -= valuesBuffer[i];
  }
  if (valuesBuffer) delete[] valuesBuffer;
  return 9999;
}

