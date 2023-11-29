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
// This class perform transmission term of volume plasmon excitation,
// based on Quinn Model, see Phys. Rev. vol 126, number 4 (1962)

#include "G4DNAQuinnPlasmonExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAQuinnPlasmonExcitationModel::G4DNAQuinnPlasmonExcitationModel
                        (const G4ParticleDefinition*,
                         const G4String& nam):
G4VEmModel(nam), isInitialised(false)
{
  fpMaterialDensity       = 0;
  fLowEnergyLimit         =  10  *  eV;
  fHighEnergyLimit        =  1.0 * GeV;

  for(G4int i=0;i<100;i++) nValenceElectron[i]=0;

  verboseLevel = 0;

  if (verboseLevel > 0)
  {
    G4cout << "Quinn plasmon excitation model is constructed " << G4endl;
  }
  fParticleChangeForGamma = 0;
  statCode                = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAQuinnPlasmonExcitationModel::~G4DNAQuinnPlasmonExcitationModel()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAQuinnPlasmonExcitationModel::Initialise
                        (const G4ParticleDefinition* particle,
                         const G4DataVector& /*cuts*/)
{
  for(G4int i=0;i<100;i++) nValenceElectron[i]=0;

  if (verboseLevel > 3)
  {
    G4cout << 
           "Calling G4DNAQuinnPlasmonExcitationModel::Initialise()" 
           << G4endl;
  }

  if(particle == G4Electron::ElectronDefinition())
  {
    fLowEnergyLimit           =  10  *  eV;
    fHighEnergyLimit          =  1.0 * GeV;
  }
  else
  { 
    G4Exception("G4DNAQuinnPlasmonExcitationModel::Initialise","em0001",
    FatalException,"Not defined for other particles than electrons.");
   return;
  }

  // Get Number of valence electrons
  G4ProductionCutsTable* theCoupleTable = 
      G4ProductionCutsTable::GetProductionCutsTable();
      
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
  
  for(G4int i=0;i<numOfCouples;i++){
    
    const G4MaterialCutsCouple* couple = 
          theCoupleTable->GetMaterialCutsCouple(i);
    
    const G4Material* material = couple->GetMaterial();
    
    const G4ElementVector* theElementVector =material->GetElementVector();
    
    std::size_t nelm = material->GetNumberOfElements();
    if (nelm==1) // Protection: only for single element
    {
      G4int z = G4lrint((*theElementVector)[0]->GetZ());
      if(z<=100)
      {
        nValenceElectron[z]  = GetNValenceElectron(z);
      }
      else
      {
        G4Exception("G4DNAQuinnPlasmonExcitationModel::Initialise","em0002",
          FatalException,"The model is not applied for z>100");
      }
    }
    //for(G4int j=0;j<nelm;j++){
    //    G4int z=G4lrint((*theElementVector)[j]->GetZ());
    //    if(z<=100){nValenceElectron[z]  = GetNValenceElectron(z);}
    //}
  }

  if( verboseLevel>0 )
  {
    G4cout << "Quinn plasmon excitation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / keV << " keV for "
    << particle->GetParticleName()
    << G4endl;
  }

  if (isInitialised){return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAQuinnPlasmonExcitationModel::CrossSectionPerVolume
                         (const G4Material* material,
                          const G4ParticleDefinition* particleDefinition,
                          G4double ekin,
                          G4double,
                          G4double)
{
  if (verboseLevel > 3)
  {
    G4cout << 
         "Calling CrossSectionPerVolume() of G4DNAQuinnPlasmonExcitationModel"
           << G4endl;
  }

  // Protection: only for single element
  if(material->GetNumberOfElements()>1) return 0.;
  G4double z = material->GetZ();

  // Protection: only for Gold
  if (z!=79){return 0.;}
  

  G4double sigma           = 0;
  G4double atomicNDensity = material->GetAtomicNumDensityVector()[0];

  if(atomicNDensity!= 0.0)
  {
    if (ekin >= fLowEnergyLimit && ekin < fHighEnergyLimit)
    {
      sigma = GetCrossSection(material,particleDefinition,ekin);
    }

    if (verboseLevel > 2)
    {
      G4cout<<"__________________________________" << G4endl;
      G4cout<<"=== G4DNAQuinnPlasmonExcitationModel - XS INFO START"<<G4endl;
      G4cout<<"=== Kinetic energy (eV)=" << ekin/eV << " particle : " 
            <<particleDefinition->GetParticleName() << G4endl;
      G4cout<<"=== Cross section per atom for Z="<<z<<" is (cm^2)" 
            <<sigma/cm/cm << G4endl;
      G4cout<<"=== Cross section per atom for Z="<<z<<" is (cm^-1)=" 
            <<sigma*atomicNDensity/(1./cm) << G4endl;
      G4cout<<"=== G4DNAQuinnPlasmonExcitationModel - XS INFO END" << G4endl;
    }
  } 

  return sigma*atomicNDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAQuinnPlasmonExcitationModel::SampleSecondaries
                         (std::vector<G4DynamicParticle*>* /*fvect*/,
                          const G4MaterialCutsCouple* couple,
                          const G4DynamicParticle* aDynamicParticle,
                          G4double,G4double)
{

  if (verboseLevel > 3)
  {
    G4cout << 
           "Calling SampleSecondaries() of G4DNAQuinnPlasmonExcitationModel"
           << G4endl;
  }

  const G4Material *material = couple->GetMaterial();

  G4ParticleDefinition* particle = aDynamicParticle->GetDefinition();
  
  G4double k                     = aDynamicParticle->GetKineticEnergy();

  if(particle == G4Electron::ElectronDefinition())
  {
    G4double  e      = 1.;
    G4int     z      = material->GetZ();
    G4int     Nve    = 0;

    //TODO: have to be change to realistic!!
    if(z<100) Nve    = nValenceElectron[z];

    G4double  A      = material->GetA()/g/mole;
    G4double  Dens   = material->GetDensity()/g*cm*cm*cm;
    G4double  veDens = Dens*CLHEP::Avogadro*Nve/A;

    G4double omega_p = std::sqrt(veDens*std::pow(e,2)/
                      (CLHEP::epsilon0/(1./cm)*CLHEP::electron_mass_c2
                      /(CLHEP::c_squared/cm/cm)));
    
    G4double excitationEnergy   = CLHEP::hbar_Planck*omega_p;
    G4double newEnergy          = k - excitationEnergy;


    if (newEnergy > 0)
    {
      fParticleChangeForGamma->
            ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
      
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
      
      if(!statCode)
      {
           fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
      }
      else
      {
           fParticleChangeForGamma->SetProposedKineticEnergy(k);
      
      }
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAQuinnPlasmonExcitationModel::GetCrossSection
                         (const G4Material* material,
                          const G4ParticleDefinition* particle,
                          G4double kineticEnergy)
{
  G4double value=0; 

  if(particle == G4Electron::ElectronDefinition())
  {
     G4double  e      = 1.;
     G4int     z      = material->GetZ();
     G4int     Nve    = 0;
     if(z<100) Nve    = nValenceElectron[z];
     G4double  A      = material->GetA()/g/mole;
     G4double  Dens   = material->GetDensity()/g*cm*cm*cm;
     G4double  veDens = Dens*CLHEP::Avogadro*Nve/A;

     G4double omega_p = std::sqrt(veDens*std::pow(e,2)
                        /(CLHEP::epsilon0/(1./cm)*CLHEP::electron_mass_c2/
                        (CLHEP::c_squared/cm/cm)));
     
     G4double fEnergy = std::pow(CLHEP::h_Planck,2)/(8*CLHEP::electron_mass_c2)*
                        std::pow(3*veDens/CLHEP::pi,2./3.)/e
                        *(CLHEP::c_squared/cm/cm);
     
     G4double p0      = sqrt(2*CLHEP::electron_mass_c2
                        /(CLHEP::c_squared/cm/cm)*fEnergy);
     
     G4double p       = sqrt(2*CLHEP::electron_mass_c2
                        /(CLHEP::c_squared/cm/cm)*kineticEnergy);

     G4double mfp     = 2*CLHEP::Bohr_radius/cm*kineticEnergy
                        /(CLHEP::hbar_Planck*omega_p)/
                        (G4Log((std::pow(std::pow(p0,2)
                        +2*CLHEP::electron_mass_c2/
                        (CLHEP::c_squared/cm/cm)*omega_p
                        *CLHEP::hbar_Planck,1./2.)-p0)
                        /(p-std::pow(std::pow(p,2)-2*CLHEP::electron_mass_c2/
                        (CLHEP::c_squared/cm/cm)*omega_p
                        *CLHEP::hbar_Planck,1./2.))));

     G4double excitationEnergy   = CLHEP::hbar_Planck*omega_p;
     
     if((0<mfp)&&(0<veDens)&&(excitationEnergy<kineticEnergy)){
       value            = 1./(veDens*mfp);
     }
  }
  return value*cm*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNAQuinnPlasmonExcitationModel::GetNValenceElectron(G4int z)
{

  G4int Nve=0;
  
  // Current limitation to gold
  if (z!=79){return 0.;}

  if (verboseLevel > 3)
  {
    G4cout << 
           "Calling GetNValenceElectron() of G4DNAQuinnPlasmonExcitationModel"
           << G4endl;
  }
 
  const char *datadir=0;
  
  if(!datadir)
  {
     datadir = G4FindDataDir("G4LEDATA");
     if(!datadir)
     {
       G4Exception("G4DNAQuinnPlasmonExcitationModel::GetNValenceElectron()"
                   ,"em0002",FatalException,
                   "Enviroment variable G4LEDATA not defined");
       return 0;
     }
  }

  std::ostringstream targetfile;
  targetfile.str("");
  targetfile.clear(stringstream::goodbit);
  targetfile << datadir <<"/dna/atomicstate_Z"<< z <<".dat";
  std::ifstream fin(targetfile.str().c_str());

  if(!fin)
  {
    G4cout<< " Error : "<< targetfile.str() <<" is not found "<<endl;
    G4Exception("G4DNAQuinnPlasmonExcitationModel::GetNValenceElectron()"
                ,"em0003",FatalException,
                "There is no target file");
    return 0;
  }
  
  string buff0,buff1,buff2,buff3,buff4,buff5,buff6;
  fin >> buff0 >>buff1>>buff2>>buff3>>buff4>>buff5>>buff6;
  
  while(true){
    fin >> buff0 >>buff1>>buff2>>buff3>>buff4>>buff5>>buff6;
    if(!fin.eof())
    {
      Nve = stoi(buff3);
    }
    else
    {
      break;
    }
  }
  return Nve;
}

