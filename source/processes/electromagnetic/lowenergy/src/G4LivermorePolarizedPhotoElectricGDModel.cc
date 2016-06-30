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
// $Id: G4LivermorePolarizedPhotoElectricGDModel.hh,v 1.2 2010-11-23 16:42:15 flongo Exp $
//
// Authors: G.Depaola & F.Longo
//

#include "G4LivermorePolarizedPhotoElectricGDModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4ParticleChangeForGamma.hh" 

#include <vector>

G4LPhysicsFreeVector*  G4LivermorePolarizedPhotoElectricGDModel::fCrossSection[] = {nullptr};
G4LPhysicsFreeVector*  G4LivermorePolarizedPhotoElectricGDModel::fCrossSectionLE[] = {nullptr};
std::vector<G4double>* G4LivermorePolarizedPhotoElectricGDModel::fParam[] = {0};
G4int                  G4LivermorePolarizedPhotoElectricGDModel::fNShells[] = {0};
G4int                  G4LivermorePolarizedPhotoElectricGDModel::fNShellsUsed[] = {0};
G4ElementData*         G4LivermorePolarizedPhotoElectricGDModel::fShellCrossSection = nullptr;
G4Material*            G4LivermorePolarizedPhotoElectricGDModel::fWater = nullptr;
G4double               G4LivermorePolarizedPhotoElectricGDModel::fWaterEnergyLimit = 0.0;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedPhotoElectricGDModel::G4LivermorePolarizedPhotoElectricGDModel(
   const G4String& nam)
  :G4VEmModel(nam),fParticleChange(nullptr),maxZ(99),
   nShellLimit(100), fDeexcitationActive(false), isInitialised(false), 
   fAtomDeexcitation(nullptr)
{
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  theGamma    = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();

  SetDeexcitationFlag(true);
  fSandiaCof.resize(4,0.0);
  fCurrSection = 0.0;
  
  if (verboseLevel > 0) {
    G4cout << "Livermore Polarized PhotoElectric is constructed " 
	   << " nShellLimit "
	   << nShellLimit << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedPhotoElectricGDModel::~G4LivermorePolarizedPhotoElectricGDModel()
{  
  if(IsMaster()) {
    delete fShellCrossSection;
    for(G4int i=0; i<maxZ; ++i) { 
      delete fParam[i];
      fParam[i] = 0;
      delete fCrossSection[i];
      fCrossSection[i] = 0;
      delete fCrossSectionLE[i];
      fCrossSectionLE[i] = 0;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePolarizedPhotoElectricGDModel::Initialise(
                    const G4ParticleDefinition*,
		    const G4DataVector&)
{
  if (verboseLevel > 2) {
    G4cout << "Calling G4LivermorePolarizedPhotoElectricGDModel::Initialise()" << G4endl;
  }

  if(IsMaster()) {
    
    if(!fWater) { 
      fWater = G4Material::GetMaterial("G4_WATER", false); 
      if(fWater) { fWaterEnergyLimit = 13.6*eV; }
    }
    
    if(!fShellCrossSection) { fShellCrossSection = new G4ElementData(); }
    
    char* path = getenv("G4LEDATA");
    
    G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = theCoupleTable->GetTableSize();
    
    for(G4int i=0; i<numOfCouples; ++i) {
      const G4MaterialCutsCouple* couple = 
	theCoupleTable->GetMaterialCutsCouple(i);
      const G4Material* material = couple->GetMaterial();
      const G4ElementVector* theElementVector = material->GetElementVector();
      G4int nelm = material->GetNumberOfElements();
      
      for (G4int j=0; j<nelm; ++j) {        
	G4int Z = (G4int)(*theElementVector)[j]->GetZ();
	if(Z < 1)          { Z = 1; }
	else if(Z > maxZ)  { Z = maxZ; }
	if(!fCrossSection[Z]) { ReadData(Z, path); }
      }
    }
  }  
  //  
  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files for LivermorePhotoElectric model" 
	   << G4endl;
  }
  if(!isInitialised) {
    isInitialised = true;
    fParticleChange = GetParticleChangeForGamma();
    
    fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  }
  fDeexcitationActive = false;
  if(fAtomDeexcitation) { 
    fDeexcitationActive = fAtomDeexcitation->IsFluoActive(); 
  }
  
  if (verboseLevel > 0) { 
    G4cout << "LivermorePolarizedPhotoElectric model is initialized " << G4endl
	   << G4endl;
  }
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4double G4LivermorePolarizedPhotoElectricGDModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double ZZ, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "G4LivermorePolarizedPhotoElectricGDModel::ComputeCrossSectionPerAtom():" 
	   << " Z= " << ZZ << "  R(keV)= " << GammaEnergy/keV << G4endl;
  }
  G4double cs = 0.0;
  G4int Z = G4lrint(ZZ);
  if(Z < 1 || Z >= maxZ) { return cs; }
  
  // if element was not initialised
  // do initialisation safely for MT mode
  if(!fCrossSection[Z]) {
    InitialiseForElement(0, Z);
    if(!fCrossSection[Z]) { return cs; }
  }
  
  G4int idx = fNShells[Z]*6 - 4;
  if (GammaEnergy < (*(fParam[Z]))[idx-1]) { GammaEnergy = (*(fParam[Z]))[idx-1]; }
  
  G4double x1 = 1.0/GammaEnergy;
  G4double x2 = x1*x1;
  G4double x3 = x2*x1;
  
  // parameterisation
  if(GammaEnergy >= (*(fParam[Z]))[0]) {
    G4double x4 = x2*x2;
    cs = x1*((*(fParam[Z]))[idx] + x1*(*(fParam[Z]))[idx+1]
	     + x2*(*(fParam[Z]))[idx+2] + x3*(*(fParam[Z]))[idx+3] 
	     + x4*(*(fParam[Z]))[idx+4]);
    // high energy part
  } else if (GammaEnergy >= (*(fParam[Z]))[1]) {
    cs = x3*(fCrossSection[Z])->Value(GammaEnergy);
    
    // low energy part
  } else {
    cs = x3*(fCrossSectionLE[Z])->Value(GammaEnergy);
  }
  if (verboseLevel > 1) { 
    G4cout << "LivermorePolarizedPhotoElectricGDModel: E(keV)= " << GammaEnergy/keV
           << " Z= " << Z << " cross(barn)= " << cs/barn << G4endl;
  }
  return cs;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedPhotoElectricGDModel::SampleSecondaries(
                                std::vector<G4DynamicParticle*>* fvect,
				const G4MaterialCutsCouple* couple,
				const G4DynamicParticle* aDynamicGamma,
				G4double,
				G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling SampleSecondaries() of G4LivermorePolarizedPhotoElectricGDModel" 
	   << G4endl;
  }

   G4double photonEnergy = aDynamicGamma->GetKineticEnergy();
  if (verboseLevel > 3) {
    G4cout << "G4LivermorePolarizedPhotoElectricGDModel::SampleSecondaries() Egamma(keV)= "
	   << photonEnergy/keV << G4endl;
  }
  
  G4ThreeVector gammaPolarization0 = aDynamicGamma->GetPolarization();  
  G4ThreeVector photonDirection = aDynamicGamma->GetMomentumDirection();
  
  // kill incident photon
  fParticleChange->ProposeTrackStatus(fStopAndKill);   
  fParticleChange->SetProposedKineticEnergy(0.);
  
  // low-energy photo-effect in water - full absorption

  const G4Material* material = couple->GetMaterial();
  if(fWater && (material == fWater || 
		material->GetBaseMaterial() == fWater)) {
    if(photonEnergy <= fWaterEnergyLimit) { 
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);      
      return;
    }
  }   

  // Protection: a polarisation parallel to the
  // direction causes problems;
  // in that case find a random polarization
  
  // Make sure that the polarization vector is perpendicular to the
  // gamma direction. If not
  
  if(!(gammaPolarization0.isOrthogonal(photonDirection, 1e-6))||(gammaPolarization0.mag()==0))
    { // only for testing now
      gammaPolarization0 = GetRandomPolarization(photonDirection);
    }
  else
    {
      if ( gammaPolarization0.howOrthogonal(photonDirection) != 0)
	{
	  gammaPolarization0 = GetPerpendicularPolarization(photonDirection, gammaPolarization0);
	}
    }
  
  // End of Protection
  
  //  G4double E0_m = photonEnergy / electron_mass_c2 ;

  // Shell 

  // Select randomly one element in the current material
  //G4cout << "Select random atom Egamma(keV)= " << photonEnergy/keV << G4endl;
  const G4Element* elm = SelectRandomAtom(material, theGamma, photonEnergy);
  G4int Z = G4lrint(elm->GetZ());


  // Select the ionised shell in the current atom according to shell 
  //   cross sections
  // G4cout << "Select random shell Z= " << Z << G4endl;
  
  if(Z >= maxZ) { Z = maxZ-1; }
  
  // element was not initialised gamma should be absorbed
  if(!fCrossSection[Z]) {
    fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
    return;
  }
  
  // shell index
  size_t shellIdx = 0;
  size_t nn = fNShellsUsed[Z];
  
  if(nn > 1) {
    if(photonEnergy >= (*(fParam[Z]))[0]) {
      G4double x1 = 1.0/photonEnergy;
      G4double x2 = x1*x1;
      G4double x3 = x2*x1;
      G4double x4 = x3*x1;
      G4int idx   = nn*6 - 4;
      // when do sampling common factors are not taken into account
      // so cross section is not real
      G4double cs0 = G4UniformRand()*((*(fParam[Z]))[idx] 
				      + x1*(*(fParam[Z]))[idx+1]
				      + x2*(*(fParam[Z]))[idx+2] 
				      + x3*(*(fParam[Z]))[idx+3] 
				      + x4*(*(fParam[Z]))[idx+4]);
      for(shellIdx=0; shellIdx<nn; ++shellIdx) {
	idx = shellIdx*6 + 2;
	if(photonEnergy > (*(fParam[Z]))[idx-1]) {
	  G4double cs = (*(fParam[Z]))[idx] + x1*(*(fParam[Z]))[idx+1] 
	    + x2*(*(fParam[Z]))[idx+2] + x3*(*(fParam[Z]))[idx+3] 
	    + x4*(*(fParam[Z]))[idx+4];
	  if(cs >= cs0) { break; }
	}
      }
      if(shellIdx >= nn) { shellIdx = nn-1; }
      
    } else {
      
      // when do sampling common factors are not taken into account
      // so cross section is not real
      G4double cs = G4UniformRand();
      
      if(photonEnergy >= (*(fParam[Z]))[1]) {
	cs *= (fCrossSection[Z])->Value(photonEnergy);
      } else {
	cs *= (fCrossSectionLE[Z])->Value(photonEnergy);
      }
      
      for(size_t j=0; j<nn; ++j) {
	shellIdx = (size_t)fShellCrossSection->GetComponentID(Z, j);
	if(photonEnergy > (*(fParam[Z]))[6*shellIdx+1]) {
	  cs -= fShellCrossSection->GetValueForComponent(Z, j, photonEnergy);
	}
	if(cs <= 0.0 || j+1 == nn) { break; }
      }
    }
  }
  
  G4double bindingEnergy = (*(fParam[Z]))[shellIdx*6 + 1];
  //G4cout << "Z= " << Z << " shellIdx= " << shellIdx 
  //       << " nShells= " << fNShells[Z] 
  //       << " Ebind(keV)= " << bindingEnergy/keV 
  //       << " Egamma(keV)= " << photonEnergy/keV << G4endl;
  
  const G4AtomicShell* shell = 0;
  
  // no de-excitation from the last shell
  if(fDeexcitationActive && shellIdx + 1 < nn) {
    G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIdx);
    shell = fAtomDeexcitation->GetAtomicShell(Z, as);
  }
  
  // If binding energy of the selected shell is larger than photon energy
  //    do not generate secondaries
  if(photonEnergy < bindingEnergy) {
    fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
    return;
  }


  // Electron 
  
  G4double eKineticEnergy = photonEnergy - bindingEnergy;
  G4double edep = bindingEnergy;
  
  G4double costheta = SetCosTheta(eKineticEnergy);
  G4double sintheta = sqrt(1. - costheta*costheta);
  G4double phi = SetPhi(photonEnergy,eKineticEnergy,costheta);
  G4double dirX = sintheta*cos(phi);
  G4double dirY = sintheta*sin(phi);
  G4double dirZ = costheta;
  G4ThreeVector electronDirection(dirX, dirY, dirZ);
  SystemOfRefChange(photonDirection, electronDirection, gammaPolarization0);
  G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(), 
						       electronDirection, 
						       eKineticEnergy);
  fvect->push_back(electron);
  
  // Deexcitation
  //  Sample deexcitation
  if(shell) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      G4int nbefore = fvect->size();
      
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      G4int nafter = fvect->size();
      if(nafter > nbefore) {
	G4double esec = 0.0;
	for (G4int j=nbefore; j<nafter; ++j) {
	  
	  G4double e = ((*fvect)[j])->GetKineticEnergy();
	  if(esec + e > edep) {
	    // correct energy in order to have energy balance
            e = edep - esec;
	    ((*fvect)[j])->SetKineticEnergy(e);
	    esec += e;
	    // delete the rest of secondaries (should not happens)
            for (G4int jj=nafter-1; jj>j; --jj) { 
              delete (*fvect)[jj]; 
              fvect->pop_back(); 
            }
	    break;        
	  }
	  esec += e; 
	}
	edep -= esec;
      }
    }
  }
  // energy balance - excitation energy left
  if(edep > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedPhotoElectricGDModel::ReadData(G4int Z, const char* path)
{
  if (verboseLevel > 1) 
    {
      G4cout << "Calling ReadData() of G4LivermorePolarizedPhotoElectricGDModel" 
	     << G4endl;
    }
  
  if(fCrossSection[Z]) { return; }
  
  const char* datadir = path;
  
  if(!datadir) 
    {
      datadir = getenv("G4LEDATA");
      if(!datadir) 
	{
	  G4Exception("G4LivermorePhotoElectricModel::ReadData()",
		      "em0006",FatalException,
		      "Environment variable G4LEDATA not defined");
	  return;
	}
    }
  
  // spline for photoeffect total x-section above K-shell
  fCrossSection[Z] = new G4LPhysicsFreeVector();
  fCrossSection[Z]->SetSpline(true);
  
  std::ostringstream ost;
  ost << datadir << "/livermore/phot/pe-cs-" << Z <<".dat";
  std::ifstream fin(ost.str().c_str());
  if( !fin.is_open()) {
    G4ExceptionDescription ed;
    ed << "G4LivermorePolarizedPhotoElectricGDModel data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4LivermorePolarizedPhotoElectricGDModel::ReadData()",
		"em0003",FatalException,
		ed,"G4LEDATA version should be G4EMLOW6.32 or later.");
    return;
  } else {
    if(verboseLevel > 3) { G4cout << "File " << ost.str().c_str() 
				  << " is opened by G4LivermorePolarizedPhotoElectricGDModel" << G4endl;}
    fCrossSection[Z]->Retrieve(fin, true);
    fCrossSection[Z]->ScaleVector(MeV, barn);
    fin.close();
  }
  
  fParam[Z] = new std::vector<G4double>;
  
  // read fit parameters
  G4int n1 = 0;
  G4int n2 = 0;
  G4double x;
  std::ostringstream ost1;
  ost1 << datadir << "/livermore/phot/pe-" << Z <<".dat";
  std::ifstream fin1(ost1.str().c_str());
  if( !fin1.is_open()) {
    G4ExceptionDescription ed;
    ed << "G4LivermorePolarizedPhotoElectricGDModel data file <" << ost1.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4LivermorePolarizedPhotoElectricGDModel::ReadData()",
		"em0003",FatalException,
		ed,"G4LEDATA version should be G4EMLOW6.32 or later.");
    return;
  } else {
    if(verboseLevel > 3) { 
      G4cout << "File " << ost1.str().c_str()
	     << " is opened by G4LivermorePolarizedPhotoElectricGDModel" << G4endl;
    }
    fin1 >> n1;
    if(fin1.fail()) { return; }
    if(0 > n1 || n1 >= INT_MAX) { n1 = 0; }
    
    fin1 >> n2;
    if(fin1.fail()) { return; }
    if(0 > n2 || n2 >= INT_MAX) { n2 = 0; }
    
    fin1 >> x;
    if(fin1.fail()) { return; }
    
    fNShells[Z] = n1;
    fParam[Z]->reserve(6*n1+1);
    fParam[Z]->push_back(x*MeV);
    for(G4int i=0; i<n1; ++i) {
      for(G4int j=0; j<6; ++j) {
        fin1 >> x;
	if(0 == j) { x *= MeV; }
	else       { x *= barn; }
        fParam[Z]->push_back(x);
      }
    }
    fin1.close();
  }
  // there is a possibility to used only main shells
  if(nShellLimit < n2) { n2 = nShellLimit; }
  fShellCrossSection->InitialiseForComponent(Z, n2);
  fNShellsUsed[Z] = n2;
  
  if(1 < n2) {
    std::ostringstream ost2;
    ost2 << datadir << "/livermore/phot/pe-ss-cs-" << Z <<".dat";
    std::ifstream fin2(ost2.str().c_str());
    if( !fin2.is_open()) {
      G4ExceptionDescription ed;
      ed << "G4LivermorePolarizedPhotoElectricGDModel data file <" << ost2.str().c_str()
	 << "> is not opened!" << G4endl;
      G4Exception("G4LivermorePolarizedPhotoElectricGDModel::ReadData()",
		  "em0003",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW6.32 or later.");
      return;
    } else {
      if(verboseLevel > 3) { 
	G4cout << "File " << ost2.str().c_str()
	       << " is opened by G4LivermorePolarizedPhotoElectricGDModel" << G4endl;
      }
      
      G4int n3, n4;
      G4double y;
      for(G4int i=0; i<n2; ++i) {
	fin2 >> x >> y >> n3 >> n4;
	G4LPhysicsFreeVector* v = new G4LPhysicsFreeVector(n3, x, y);
	for(G4int j=0; j<n3; ++j) {
	  fin2 >> x >> y;
	  v->PutValues(j, x*MeV, y*barn);
	}
	fShellCrossSection->AddComponent(Z, n4, v);
      }
      fin2.close();
    }
  }
  
  // no spline for photoeffect total x-section below K-shell
  if(1 < fNShells[Z]) {
    fCrossSectionLE[Z] = new G4LPhysicsFreeVector();
    
    std::ostringstream ost3;
    ost3 << datadir << "/livermore/phot/pe-le-cs-" << Z <<".dat";
    std::ifstream fin3(ost3.str().c_str());
    if( !fin3.is_open()) {
      G4ExceptionDescription ed;
      ed << "G4LivermorePolarizedPhotoElectricGDModel data file <" << ost3.str().c_str()
         << "> is not opened!" << G4endl;
      G4Exception("G4LivermorePolarizedPhotoElectricGDModel::ReadData()",
		  "em0003",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW6.32 or later.");
      return;
    } else {
      if(verboseLevel > 3) { 
        G4cout << "File " << ost3.str().c_str() 
	       << " is opened by G4LivermorePolarizedPhotoElectricGDModel" << G4endl;
      }
      fCrossSectionLE[Z]->Retrieve(fin3, true);
      fCrossSectionLE[Z]->ScaleVector(MeV, barn);
      fin3.close();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4LivermorePolarizedPhotoElectricGDModel::SetCosTheta(G4double energyE)
{
  G4double rand1,rand2,onemcost,greject;
  G4double masarep = 510.99906*keV;
  
  G4double gamma = 1. + energyE/masarep;
  G4double gamma2 = gamma*gamma;
  
  G4double beta = sqrt((gamma2 - 1.)/gamma2);
  
  G4double alfa = 1./beta - 1.;
  
  G4double g1 = 0.5*beta*gamma*(gamma-1.)*(gamma-2.);
  
  G4double alfap2 = alfa+2.;
  
  G4double grejectmax = 2.*(g1+1./alfa);
  
  do
    {
      rand1 = G4UniformRand();
      onemcost = 2.*alfa*(2.*rand1 + alfap2 * sqrt(rand1))/
	(alfap2*alfap2 - 4.*rand1);
      greject = (2. - onemcost)*(g1+1./(alfa+onemcost));
      rand2 = G4UniformRand();
    }
  while (rand2*grejectmax > greject);
  G4double cosTheta = 1. - onemcost;
  return cosTheta;
} 

 
G4double G4LivermorePolarizedPhotoElectricGDModel::SetPhi(G4double Ph_energy,
							  G4double E_energy,
							  G4double costheta)
{
  G4double epsilon = E_energy/electron_mass_c2;
  G4double k = Ph_energy/electron_mass_c2;
  G4double gamma = 1. + epsilon;
  G4double gamma2 = gamma*gamma;
  G4double beta = sqrt((gamma2 - 1.)/gamma2);
  
  G4double d = (2./(k*gamma*(1-beta*costheta))-1)*(1/k);
  
  G4double norm_factor = 1 +2*d;
  
  G4double rnd1; 
  G4double rnd2;
  G4double phi, phiprob;
  
  do
    {
      rnd1 =G4UniformRand();
      rnd2 =G4UniformRand();
      phi = rnd1*twopi;
      phiprob = 1 +2*d*cos(phi)*cos(phi);
    }
  while (rnd2*norm_factor > phiprob);
  return phi;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LivermorePolarizedPhotoElectricGDModel::SetPerpendicularVector(G4ThreeVector& a)
{
  G4double dx = a.x();
  G4double dy = a.y();
  G4double dz = a.z();
  G4double x = dx < 0.0 ? -dx : dx;
  G4double y = dy < 0.0 ? -dy : dy;
  G4double z = dz < 0.0 ? -dz : dz;
  if (x < y) {
    return x < z ? G4ThreeVector(-dy,dx,0) : G4ThreeVector(0,-dz,dy);
  }else{
    return y < z ? G4ThreeVector(dz,0,-dx) : G4ThreeVector(-dy,dx,0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LivermorePolarizedPhotoElectricGDModel::GetRandomPolarization(G4ThreeVector& direction0)
{
  G4ThreeVector d0 = direction0.unit();
  G4ThreeVector a1 = SetPerpendicularVector(d0); //different orthogonal
  G4ThreeVector a0 = a1.unit(); // unit vector
  
  G4double rand1 = G4UniformRand();
  
  G4double angle = twopi*rand1; // random polar angle
  G4ThreeVector b0 = d0.cross(a0); // cross product
  
  G4ThreeVector c;
  
  c.setX(std::cos(angle)*(a0.x())+std::sin(angle)*b0.x());
  c.setY(std::cos(angle)*(a0.y())+std::sin(angle)*b0.y());
  c.setZ(std::cos(angle)*(a0.z())+std::sin(angle)*b0.z());
  
  G4ThreeVector c0 = c.unit();
  
  return c0;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LivermorePolarizedPhotoElectricGDModel::GetPerpendicularPolarization
(const G4ThreeVector& gammaDirection, const G4ThreeVector& gammaPolarization) const
{
  
  // 
  // The polarization of a photon is always perpendicular to its momentum direction.
  // Therefore this function removes those vector component of gammaPolarization, which
  // points in direction of gammaDirection
  //

  // Mathematically we search the projection of the vector a on the plane E, where n is the
     // plains normal vector.
     // The basic equation can be found in each geometry book (e.g. Bronstein):
     // p = a - (a o n)/(n o n)*n
   
  return gammaPolarization - gammaPolarization.dot(gammaDirection)/gammaDirection.dot(gammaDirection) * gammaDirection;  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 
void G4LivermorePolarizedPhotoElectricGDModel::SystemOfRefChange
(G4ThreeVector& direction0,G4ThreeVector& direction1,
 G4ThreeVector& polarization0)
{
  // direction0 is the original photon direction ---> z
  // polarization0 is the original photon polarization ---> x
  // need to specify y axis in the real reference frame ---> y 
  G4ThreeVector Axis_Z0 = direction0.unit();
  G4ThreeVector Axis_X0 = polarization0.unit();
  G4ThreeVector Axis_Y0 = (Axis_Z0.cross(Axis_X0)).unit(); // to be confirmed;
  
  G4double direction_x = direction1.getX();
  G4double direction_y = direction1.getY();
  G4double direction_z = direction1.getZ();
  
  direction1 = (direction_x*Axis_X0 + direction_y*Axis_Y0 +  direction_z*Axis_Z0).unit();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AutoLock.hh"
namespace { G4Mutex LivermorePolarizedPhotoElectricGDModelMutex = G4MUTEX_INITIALIZER; }

void G4LivermorePolarizedPhotoElectricGDModel::InitialiseForElement(
								    const G4ParticleDefinition*, G4int Z)
{
  G4AutoLock l(&LivermorePolarizedPhotoElectricGDModelMutex);
  //  G4cout << "G4LivermorePhotoElectricModel::InitialiseForElement Z= " 
  //   << Z << G4endl;
  if(!fCrossSection[Z]) { ReadData(Z); }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 

