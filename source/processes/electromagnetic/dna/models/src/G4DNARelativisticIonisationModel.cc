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
// $Id: G4DNARelativisticIonisationModel.cc $
//
// Created on 2016/05/12
//
// Authors: D Sakata, S. Incerti
//
// This class perform ionisation for electron transportation in gold,
// based on Relativistic Binary Encounter Bethe-Vriens(RBEBV) model.
// See following reference paper, 
// M. Guerra et al, J. Phys. B: At. Mol. Opt. Phys. 48, 185202 (2015)
// =======================================================================
// Limitation of secondaries by GEANT4 atomic de-excitation:
// The cross section and energy of secondary production is based on 
// EADL database. If there are no tabele for several orbitals, this class
// will not provide secondaries for the orbitals.
// For gold(Au), this class provide secondaries for inner 18 orbitals
// but don't provide for outer 3 orbitals due to EADL databese limitation. 
// =======================================================================

#include "G4DNARelativisticIonisationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4AtomicShell.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4Gamma.hh"
#include "G4RandomDirection.hh"

#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARelativisticIonisationModel::G4DNARelativisticIonisationModel(
                const G4ParticleDefinition*,
                const G4String& nam) :
G4VEmModel(nam), isInitialised(false),statCode(false),fasterCode(true)
{
  fHighEnergyLimit        = 0;
  fLowEnergyLimit         = 0;

  verboseLevel = 0;

  SetDeexcitationFlag(true);
  fAtomDeexcitation       = 0;
  fMaterialDensity        = 0;
  fParticleDefinition     = 0;
  fParticleChangeForGamma = 0;

  if (verboseLevel > 0)
  {
    G4cout << "Relativistic Ionisation Model is constructed " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARelativisticIonisationModel::~G4DNARelativisticIonisationModel()
{
  // Cross section
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARelativisticIonisationModel::Initialise(const G4ParticleDefinition* particle,
                                           const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
  {
    G4cout << 
    "Calling G4DNARelativisticIonisationModel::Initialise()"
    << G4endl;
  }


  if(fParticleDefinition != 0 && fParticleDefinition != particle)
  {
    G4Exception("G4DNARelativisticIonisationModel::Initialise","em0001",
        FatalException,"Model already initialized for another particle type.");
  }

  fParticleDefinition = particle;
  G4ParticleDefinition *electronDef = G4Electron::ElectronDefinition();
  if(particle == electronDef)
  {
    fLowEnergyLimit           =  10  *  eV;
    fHighEnergyLimit          =  1.0 * GeV;

    std::ostringstream eFullFileNameZ;

    const char *path = G4FindDataDir("G4LEDATA");
    if (!path)
    {
      G4Exception("G4DNARelativisticIonisationModel::Initialise","em0006",
        FatalException,"G4LEDATA environment variable not set.");
      return;
    }


    G4ProductionCutsTable *coupletable 
                  = G4ProductionCutsTable::GetProductionCutsTable();
    G4int Ncouple = (G4int)coupletable ->GetTableSize();
    for(G4int i=0; i<Ncouple; ++i)
    {
      const G4MaterialCutsCouple* couple   
                           = coupletable->GetMaterialCutsCouple(i);
      const G4Material          * material = couple     ->GetMaterial();
      {
        // Protection: only for single element
        if(material->GetNumberOfElements()>1) continue; 

        G4int Z = material->GetZ();
        // Protection: only for GOLD
        if(Z!=79) continue;

        iState    [Z].clear();
        iShell    [Z].clear();
        iSubShell [Z].clear();
        Nelectrons[Z].clear();
        Ebinding  [Z].clear();
        Ekinetic  [Z].clear();
        LoadAtomicStates(Z,path);

        /////////////Load cumulated DCS////////////////
        eVecEZ.clear();
        eVecEjeEZ.clear();
        eProbaShellMapZ.clear();
        eDiffCrossSectionDataZ.clear();

        eFullFileNameZ.str("");
        eFullFileNameZ.clear(stringstream::goodbit);
              
        eFullFileNameZ
          << path
          << "/dna/sigmadiff_cumulated_ionisation_e_RBEBV_Z"
          << Z << ".dat";
        std::ifstream eDiffCrossSectionZ(eFullFileNameZ.str().c_str());
        if (!eDiffCrossSectionZ)
          G4Exception("G4DNARelativisticIonisationModel::Initialise","em0003",
          FatalException,
          "Missing data file for cumulated DCS");

        eVecEZ[Z].push_back(0.);
        while(!eDiffCrossSectionZ.eof())
        {
          G4double tDummy;
          G4double eDummy;
          eDiffCrossSectionZ>>tDummy>>eDummy;
          if (tDummy != eVecEZ[Z].back())
          {
           eVecEZ[Z].push_back(tDummy);
           eVecEjeEZ[Z][tDummy].push_back(0.);
          }

          for(G4int istate=0;istate<(G4int)iState[Z].size();istate++)
          {
            eDiffCrossSectionZ>>
            eDiffCrossSectionDataZ[Z][istate][tDummy][eDummy];
            eEjectedEnergyDataZ[Z][istate][tDummy]
            [eDiffCrossSectionDataZ[Z][istate][tDummy][eDummy]] 
            = eDummy;
            eProbaShellMapZ[Z][istate][tDummy].push_back(
            eDiffCrossSectionDataZ[Z][istate][tDummy][eDummy]);
          }

          if (eDummy != eVecEjeEZ[Z][tDummy].back()){
                             eVecEjeEZ[Z][tDummy].push_back(eDummy);
          }
        }
      }
    }
  }
  else
  { 
    G4cout<< 
    "Error : No particle Definition is found in G4DNARelativisticIonisationModel"
    <<G4endl;
    return;
  }

  if( verboseLevel>0 )
  {
    G4cout << "Relativistic Ionisation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / keV << " keV for "
    << particle->GetParticleName()
    << G4endl;
  }

  // Initialise gold density pointer
  fMaterialDensity = G4DNAMolecularMaterial::Instance()
                    ->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_Au"));

  fAtomDeexcitation       = G4LossTableManager::Instance()->AtomDeexcitation();
  fParticleChangeForGamma = GetParticleChangeForGamma();

  if (isInitialised){return;}
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNARelativisticIonisationModel::CrossSectionPerVolume(
                const G4Material* material,
                const G4ParticleDefinition* particleDefinition,
                G4double ekin,
                G4double,
                G4double)
{
  if (verboseLevel > 3)
  {
    G4cout << 
       "Calling CrossSectionPerVolume() of G4DNARelativisticIonisationModel"
    << G4endl;
  }

  if(particleDefinition != fParticleDefinition) return 0;

  // Calculate total cross section for model
  G4double sigma=0;

  if(material->GetNumberOfElements()>1) return 0.; // Protection for Molecules
  G4double atomicNDensity = material->GetAtomicNumDensityVector()[0];
  G4double z              = material->GetZ();

  if(atomicNDensity!= 0.0)
  {
    if (ekin >= fLowEnergyLimit && ekin < fHighEnergyLimit)
    {    
      sigma = GetTotalCrossSection(material,particleDefinition,ekin);
    }

    if (verboseLevel > 2)
    {
      G4cout << "__________________________________" << G4endl;
      G4cout << "=== G4DNARelativisticIonisationModel - XS INFO START" <<G4endl;
      G4cout << "=== Kinetic energy (eV)=" << ekin/eV << " particle : " 
             << particleDefinition->GetParticleName() << G4endl;
      G4cout << "=== Cross section per atom for Z="<<z<<" is (cm^2)" 
             << sigma/cm/cm << G4endl;
      G4cout << "=== Cross section per atom for Z="<<z<<" is (cm^-1)=" 
             << sigma*atomicNDensity/(1./cm) << G4endl;
      G4cout << "=== G4DNARelativisticIonisationModel - XS INFO END" << G4endl;
    }
  } 
  return sigma*atomicNDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARelativisticIonisationModel::SampleSecondaries(
                std::vector<G4DynamicParticle*>* fvect,
                const G4MaterialCutsCouple* couple,
                const G4DynamicParticle* particle,
                G4double,G4double)
{
  if (verboseLevel > 3)
  {
    G4cout << 
      "Calling SampleSecondaries() of G4DNARelativisticIonisationModel"
    << G4endl;
  }


  G4ParticleDefinition* particleDef = particle->GetDefinition();
  G4double k                        = particle->GetKineticEnergy();
  G4double ejectedE                 = 0.*eV;

  if(fLowEnergyLimit <= k && k<fHighEnergyLimit)
  {
    G4ThreeVector primaryDir        = particle   ->GetMomentumDirection(); 

    G4double      particleMass      = particleDef->GetPDGMass();
    G4double      totalEnergy       = k+particleMass;
    G4double      pSquare           = k*(totalEnergy+particleMass);
    G4double      totalMomentum     = std::sqrt(pSquare);

    const G4Material *material = couple->GetMaterial();
    G4int             z        = material->GetZ();
    G4int             level    = RandomSelect(material,particleDef,k);

    if(k<Ebinding[z].at(level)) return;

    G4int NumSecParticlesInit =0;
    G4int NumSecParticlesFinal=0;

    if(fAtomDeexcitation){
      G4AtomicShellEnumerator as = G4AtomicShellEnumerator(level);
      const G4AtomicShell *shell = fAtomDeexcitation->GetAtomicShell(z,as);
      NumSecParticlesInit  = (G4int)fvect->size();
      fAtomDeexcitation->GenerateParticles(fvect,shell,z,0,0);
      NumSecParticlesFinal = (G4int)fvect->size();
    }

    ejectedE                
            = GetEjectedElectronEnergy   (material,particleDef,k,level);
    G4ThreeVector ejectedDir
            = GetEjectedElectronDirection(particleDef,k,ejectedE);
    ejectedDir.rotateUz(primaryDir);

    G4double scatteredE = k - Ebinding[z].at(level) - ejectedE;

    if(particleDef == G4Electron::ElectronDefinition()){
      G4double secondaryTotMomentum  
            = std::sqrt(ejectedE*(ejectedE+2*CLHEP::electron_mass_c2));
      G4double finalMomentumX    
            = totalMomentum*primaryDir.x()- secondaryTotMomentum*ejectedDir.x();
      G4double finalMomentumY    
            = totalMomentum*primaryDir.y()- secondaryTotMomentum*ejectedDir.y();
      G4double finalMomentumZ    
            = totalMomentum*primaryDir.z()- secondaryTotMomentum*ejectedDir.z();
      
      G4ThreeVector scatteredDir(finalMomentumX,finalMomentumY,finalMomentumZ);
      fParticleChangeForGamma->ProposeMomentumDirection(scatteredDir.unit());

    }
    else 
    {
      fParticleChangeForGamma->ProposeMomentumDirection(primaryDir);
    }

    //G4double deexSecEnergy=0.;
    G4double restEproduction = Ebinding[z].at(level);
    for(G4int iparticle=NumSecParticlesInit;
                   iparticle<NumSecParticlesFinal;iparticle++)
    {
      //deexSecEnergy = deexSecEnergy + (*fvect)[iparticle]->GetKineticEnergy();
      G4double Edeex = (*fvect)[iparticle]->GetKineticEnergy();
      if(restEproduction>=Edeex){
        restEproduction -= Edeex;
      }
      else{
        delete (*fvect)[iparticle];
        (*fvect)[iparticle]=0;
      }
    }
    if(restEproduction < 0.0){
      G4Exception("G4DNARelativisticIonisationModel::SampleSecondaries()",
                  "em0008",FatalException,"Negative local energy deposit");
    }

    if(!statCode)
    {
      if(scatteredE>0){
        fParticleChangeForGamma->SetProposedKineticEnergy (scatteredE);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(restEproduction);
        //fParticleChangeForGamma
        //->ProposeLocalEnergyDeposit(k-scatteredE-ejectedE-deexSecEnergy);
      }
    }
    else
    {
      fParticleChangeForGamma->SetProposedKineticEnergy (k);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredE);
    }
    
    if(ejectedE>0){
      G4DynamicParticle* ejectedelectron 
             = new G4DynamicParticle(G4Electron::Electron(),ejectedDir,ejectedE);
      fvect->push_back(ejectedelectron);
    }
  }
}

void G4DNARelativisticIonisationModel::LoadAtomicStates(
                                          G4int z,const char* path)
{

  if (verboseLevel > 3)
  {
    G4cout << 
      "Calling LoadAtomicStates() of G4DNARelativisticIonisationModel"
    << G4endl;
  }
  const char *datadir =  path;
  if(!datadir)
  {
     datadir = G4FindDataDir("G4LEDATA");
     if(!datadir)
     {
       G4Exception("G4DNARelativisticIonisationModel::LoadAtomicStates()",
            "em0002",FatalException,"Enviroment variable G4LEDATA not defined");
                   
       return;
     }
  }
  std::ostringstream targetfile;
  targetfile << datadir <<"/dna/atomicstate_Z"<< z <<".dat";
  std::ifstream fin(targetfile.str().c_str());
  if(!fin)
  {
    G4cout<< " Error : "<< targetfile.str() <<" is not found "<<G4endl;
    G4Exception("G4DNARelativisticIonisationModel::LoadAtomicStates()","em0002",
                FatalException,"There is no target file");
    return;
  }

  G4String buff0,buff1,buff2,buff3,buff4,buff5,buff6;
  fin >> buff0 >>buff1>>buff2>>buff3>>buff4>>buff5>>buff6;
  G4int iline=0;
  while(true){
    fin >> buff0 >>buff1>>buff2>>buff3>>buff4>>buff5>>buff6;
    if(!fin.eof())
    {
      iState    [z].push_back(stoi(buff0));
      iShell    [z].push_back(stoi(buff1));
      iSubShell [z].push_back(stoi(buff2));
      Nelectrons[z].push_back(stoi(buff3));
      Ebinding  [z].push_back(stod(buff4));
      if(stod(buff5)==0.)
      {// if there is no kinetic energy in the file, kinetic energy 
       // for Bhor atomic model will be calculated: !!! I's not realistic!!!
        G4double radius   = std::pow(iShell[z].at(iline),2)
                  *std::pow(CLHEP::hbar_Planck,2)*(4*CLHEP::pi*CLHEP::epsilon0)
                  /CLHEP::electron_mass_c2;
        G4double momentum = iShell[z].at(iline)*CLHEP::hbar_Planck/radius;
        Ekinetic[z].push_back(std::pow(momentum,2)/(2*CLHEP::electron_mass_c2));
      }
      else
      {
        Ekinetic  [z].push_back(stod(buff5));
      }
      iline++;
    }
    else
    {
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARelativisticIonisationModel::GetTotalCrossSection(
                const G4Material* material,
                const G4ParticleDefinition* particle,
                G4double kineticEnergy)
{
  G4double value=0;
  G4int  z = material->GetZ();
  if(z!=79){ return 0.;}
  else {
    std::size_t N=iState[z].size();
    for(G4int i=0; i<(G4int)N; ++i){
      value = value+GetPartialCrossSection(material,i,particle,kineticEnergy);
    }
    return value;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARelativisticIonisationModel::GetPartialCrossSection(
                const G4Material* material,
                G4int level,
                const G4ParticleDefinition* particle,
                G4double kineticEnergy)
{
  G4double value   = 0; 
  G4double constRy =13.6057E-6;//MeV

  G4ParticleDefinition *electronDef = G4Electron::ElectronDefinition();
  G4int z = material->GetZ();
  if(particle==electronDef){

    G4double t       = kineticEnergy        /Ebinding[z].at(level);
    G4double tdash   = kineticEnergy        /CLHEP::electron_mass_c2;
    G4double udash   = Ekinetic[z].at(level)/CLHEP::electron_mass_c2;
    G4double bdash   = Ebinding[z].at(level)/CLHEP::electron_mass_c2;
    G4double beta_t2 = 1.-1./std::pow(1.+tdash,2);
    G4double beta_u2 = 1.-1./std::pow(1.+udash,2);
    G4double beta_b2 = 1.-1./std::pow(1.+bdash,2);
    G4double alpha   = std::sqrt(2*constRy/CLHEP::electron_mass_c2);
    G4double phi     = std::cos(std::sqrt(std::pow(alpha,2)
                       /(beta_t2+beta_b2))*G4Log(beta_t2/beta_b2));
    G4double constS  = 4*CLHEP::pi*std::pow(CLHEP::Bohr_radius,2)
                     *Nelectrons[z].at(level)*std::pow(alpha,4);

    if(Ebinding[z].at(level)<=kineticEnergy)
    {
      value =constS/((beta_t2+(beta_u2+beta_b2)/iShell[z].at(level))*2.*bdash)
            *(1./2.*(G4Log(beta_t2/(1.-beta_t2))-beta_t2-G4Log(2.*bdash))
            *(1.-1./std::pow(t,2.))
            +1.-1./t-G4Log(t)/(t+1.)*(1.+2.*tdash)/(std::pow(1.+tdash/2.,2.))
            *phi+std::pow(bdash,2)/(std::pow(1+tdash/2.,2))*(t-1)/2.);
    }
    
  }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARelativisticIonisationModel::GetDifferentialCrossSection(
                const G4Material* material,
                const G4ParticleDefinition* particle,
                G4double kineticEnergy,
                G4double secondaryEnergy,
                G4int  level)
{
  G4double value=0.; 
  G4double constRy =13.6057E-6;//MeV

  G4int z = material->GetZ();

  G4ParticleDefinition *electronDef = G4Electron::ElectronDefinition();
  if(particle==electronDef){
    G4double w       = secondaryEnergy      /Ebinding[z].at(level);
    G4double t       = kineticEnergy        /Ebinding[z].at(level);
    G4double tdash   = kineticEnergy        /CLHEP::electron_mass_c2;
    G4double udash   = Ekinetic[z].at(level)/CLHEP::electron_mass_c2;
    G4double bdash   = Ebinding[z].at(level)/CLHEP::electron_mass_c2;
    G4double beta_t2 = 1.-1./std::pow(1.+tdash,2);
    G4double beta_u2 = 1.-1./std::pow(1.+udash,2);
    G4double beta_b2 = 1.-1./std::pow(1.+bdash,2);
    G4double alpha   = std::sqrt(2*constRy/CLHEP::electron_mass_c2);
    G4double phi     = std::cos(std::sqrt(std::pow(alpha,2)/(beta_t2+beta_b2))
                     *G4Log(beta_t2/beta_b2));
    G4double constS  = 4*CLHEP::pi*std::pow(CLHEP::Bohr_radius,2)
                     *Nelectrons[z].at(level)*std::pow(alpha,4);

    if(secondaryEnergy<=((kineticEnergy-Ebinding[z].at(level))/2.))
    {
      value = constS/((beta_t2+(beta_u2+beta_b2)/iShell[z].at(level))*2.*bdash)
              *(-phi/(t+1.)*(1./std::pow(w+1.,1.)+1./std::pow(t-w,1.))
              *(1.+2*tdash)/std::pow(1.+tdash/2.,2.)
              +1./std::pow(w+1.,2.)+1./std::pow(t-w,2.)
              +std::pow(bdash,2)/std::pow(1+tdash/2.,2)
              +(1./std::pow(w+1.,3.)+1./std::pow(t-w,3.))
              *(G4Log(beta_t2/(1.-beta_t2))-beta_t2-G4Log(2*bdash)));
    }
  }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNARelativisticIonisationModel::RandomSelect(
                const G4Material* material,
                const G4ParticleDefinition* particle,
                G4double kineticEnergy)
{
  G4double value = 0.;
  G4int z = material->GetZ();
  std::size_t numberOfShell = iShell[z].size();
  auto valuesBuffer = new G4double[numberOfShell];
  const G4int n = (G4int)iShell[z].size();
  G4int i(n);

  while (i > 0)
  {
    i--;
    if((fLowEnergyLimit<=kineticEnergy)&&(kineticEnergy<fHighEnergyLimit))
    {
      valuesBuffer[i]=GetPartialCrossSection(material,i,particle,kineticEnergy);
    }
    value += valuesBuffer[i];
  }

  value *= G4UniformRand();
  i = n;

  while (i > 0)
  {
    i--;

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double  G4DNARelativisticIonisationModel::GetEjectedElectronEnergy(
                       const G4Material* material,
                       const G4ParticleDefinition* particle,
                       G4double energy, G4int ishell)
{
  G4double secondaryEnergy=0;

  G4ParticleDefinition *electronDef = G4Electron::ElectronDefinition();
  G4int z = material->GetZ();
  if(!fasterCode){ // for 2D rejection method
    if(particle==electronDef){
      G4double maximumsecondaryEnergy = (energy-Ebinding[z].at(ishell))/2.;
      if(maximumsecondaryEnergy<0.) return 0.;
      G4double maximumCrossSection=-999.;

      maximumCrossSection 
            = GetDifferentialCrossSection(material,particle,energy,0.,ishell);
      do{
        secondaryEnergy = G4UniformRand()* maximumsecondaryEnergy;
      }while(G4UniformRand()*maximumCrossSection > 
           GetDifferentialCrossSection(
                        material,particle,energy,secondaryEnergy,ishell));
    }
  }
  else { // for cumulative method using cumulated DCS file

    G4double valueE1  =0.;
    G4double valueE2  =0.;
    G4double valueXS21=0.;
    G4double valueXS22=0.;
    G4double valueXS11=0.;
    G4double valueXS12=0.;
    G4double ejeE21   =0.;
    G4double ejeE22   =0.;
    G4double ejeE11   =0.;
    G4double ejeE12   =0.;
    G4double random = G4UniformRand();

    if (particle == G4Electron::ElectronDefinition())
    {
      if((eVecEZ[z].at(0)<=energy)&&(energy<eVecEZ[z].back()))
      {
        std::vector<G4double>::iterator k2 
                = std::upper_bound(eVecEZ[z].begin(),eVecEZ[z].end(), energy);
        std::vector<G4double>::iterator k1 = k2-1;

        if (     random < eProbaShellMapZ[z][ishell][(*k1)].back()
              && random < eProbaShellMapZ[z][ishell][(*k2)].back() )
        {
          std::vector<G4double>::iterator xs12 =
             std::upper_bound(eProbaShellMapZ[z][ishell][(*k1)].begin(),
                              eProbaShellMapZ[z][ishell][(*k1)].end(), random);
          std::vector<G4double>::iterator xs11 = xs12-1;

          std::vector<G4double>::iterator xs22 =
             std::upper_bound(eProbaShellMapZ[z][ishell][(*k2)].begin(),
                              eProbaShellMapZ[z][ishell][(*k2)].end(), random);
          std::vector<G4double>::iterator xs21 = xs22-1;

          valueE1  =*k1;
          valueE2  =*k2;
          valueXS21 =*xs21;
          valueXS22 =*xs22;
          valueXS12 =*xs12;
          valueXS11 =*xs11;

          ejeE11 = eEjectedEnergyDataZ[z][ishell][valueE1][valueXS11];
          ejeE12 = eEjectedEnergyDataZ[z][ishell][valueE1][valueXS12];
          ejeE21 = eEjectedEnergyDataZ[z][ishell][valueE2][valueXS21];
          ejeE22 = eEjectedEnergyDataZ[z][ishell][valueE2][valueXS22];

          secondaryEnergy  = QuadInterpolator(  valueXS11, valueXS12,
                                                valueXS21, valueXS22,
                                                ejeE11 , ejeE12 ,
                                                ejeE21 , ejeE22 ,
                                                valueE1, valueE2,
                                                energy, random );
        }
      }
    }
  }

  if(secondaryEnergy<0) secondaryEnergy=0;
  return secondaryEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector  G4DNARelativisticIonisationModel::GetEjectedElectronDirection(
                             const G4ParticleDefinition* ,
                             G4double energy,G4double secondaryenergy)
{
  G4double phi       = 2*CLHEP::pi*G4UniformRand();
  G4double sintheta  = std::sqrt((1.-secondaryenergy/energy) 
                       / (1.+secondaryenergy/(2*CLHEP::electron_mass_c2)));

  G4double dirX = sintheta*std::cos(phi);
  G4double dirY = sintheta*std::sin(phi);
  G4double dirZ = std::sqrt(1.-sintheta*sintheta);

  G4ThreeVector vec(dirX,dirY,dirZ);
  return vec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARelativisticIonisationModel::Interpolate(    G4double e1,
                                                           G4double e2,
                                                           G4double e,
                                                           G4double xs1,
                                                           G4double xs2)
{

    G4double value = 0.;

    if((xs1!=0)&&(e1!=0)){
      // Log-log interpolation by default
      G4double a = (std::log10(xs2)-std::log10(xs1)) 
                  / (std::log10(e2)-std::log10(e1));
      G4double b = std::log10(xs2) - a*std::log10(e2);
      G4double sigma = a*std::log10(e) + b;
      value = (std::pow(10.,sigma));
    }
    else{
      // Lin-Lin interpolation 
      G4double d1 = xs1;
      G4double d2 = xs2;
      value = (d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
    }

    return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARelativisticIonisationModel::QuadInterpolator(
                G4double e11, G4double e12,
                G4double e21, G4double e22,
                G4double xs11, G4double xs12,
                G4double xs21, G4double xs22,
                G4double t1, G4double t2,
                G4double t, G4double e)
{
    G4double interpolatedvalue1 = Interpolate(e11, e12, e, xs11, xs12);
    G4double interpolatedvalue2 = Interpolate(e21, e22, e, xs21, xs22);
    G4double value 
             = Interpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
    return value;
}

