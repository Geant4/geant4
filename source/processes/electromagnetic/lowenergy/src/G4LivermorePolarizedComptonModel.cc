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
// $Id: G4LivermorePolarizedComptonModel.cc 95950 2016-03-03 10:42:48Z gcosmo $
//
// Authors: G.Depaola & F.Longo
//
// History:
// --------
// 02 May 2009   S Incerti as V. Ivanchenko proposed in G4LivermoreComptonModel.cc
//
// Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - remove GetMeanFreePath method and table
//                  - added protection against numerical problem in energy sampling 
//                  - use G4ElementSelector

#include "G4LivermorePolarizedComptonModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4ShellData.hh"
#include "G4DopplerProfile.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4LogLogInterpolation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4int G4LivermorePolarizedComptonModel::maxZ = 99;
G4LPhysicsFreeVector* G4LivermorePolarizedComptonModel::data[] = {0};
G4ShellData*       G4LivermorePolarizedComptonModel::shellData = 0;
G4DopplerProfile*  G4LivermorePolarizedComptonModel::profileData = 0; 
G4CompositeEMDataSet* G4LivermorePolarizedComptonModel::scatterFunctionData = 0;

//static const G4double ln10 = G4Log(10.);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedComptonModel::G4LivermorePolarizedComptonModel(const G4ParticleDefinition*, const G4String& nam)
  :G4VEmModel(nam),isInitialised(false)
{ 
  verboseLevel= 1;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if( verboseLevel>1 )  
    G4cout << "Livermore Polarized Compton is constructed " << G4endl;
        
  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
  
  fParticleChange = 0;
  fAtomDeexcitation = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedComptonModel::~G4LivermorePolarizedComptonModel()
{  
  if(IsMaster()) {
    delete shellData;
    shellData = 0;
    delete profileData;
    profileData = 0;
    delete scatterFunctionData;
    scatterFunctionData = 0;
    for(G4int i=0; i<maxZ; ++i) {
      if(data[i]) { 
	delete data[i];
	data[i] = 0;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedComptonModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& cuts)
{
  if (verboseLevel > 1)
    G4cout << "Calling G4LivermorePolarizedComptonModel::Initialise()" << G4endl;

  // Initialise element selector 

  if(IsMaster()) {
    
    // Access to elements 

    char* path = getenv("G4LEDATA");

    G4ProductionCutsTable* theCoupleTable = 
      G4ProductionCutsTable::GetProductionCutsTable();

    G4int numOfCouples = theCoupleTable->GetTableSize();
    
    for(G4int i=0; i<numOfCouples; ++i) {
      const G4Material* material = 
	theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
      const G4ElementVector* theElementVector = material->GetElementVector();
      G4int nelm = material->GetNumberOfElements();
      
      for (G4int j=0; j<nelm; ++j) {
	G4int Z = G4lrint((*theElementVector)[j]->GetZ());
	if(Z < 1)        { Z = 1; }
	else if(Z > maxZ){ Z = maxZ; }
	
	if( (!data[Z]) ) { ReadData(Z, path); }
      }
    }
    
    // For Doppler broadening
    if(!shellData) {
      shellData = new G4ShellData(); 
      shellData->SetOccupancyData();
      G4String file = "/doppler/shell-doppler";
      shellData->LoadData(file);
    }
    if(!profileData) { profileData = new G4DopplerProfile(); }

    // Scattering Function 
    
    if(!scatterFunctionData)
      {
	
	G4VDataSetAlgorithm* scatterInterpolation = new G4LogLogInterpolation;
	G4String scatterFile = "comp/ce-sf-";
	scatterFunctionData = new G4CompositeEMDataSet(scatterInterpolation, 1., 1.);
	scatterFunctionData->LoadData(scatterFile);
      }
    
    InitialiseElementSelectors(particle, cuts);
  }
 
  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files" << G4endl;
  }
  
  if( verboseLevel>1 ) { 
    G4cout << "G4LivermoreComptonModel is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }
  //  
  if(isInitialised) { return; }
  
  fParticleChange = GetParticleChangeForGamma();
  fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();
  isInitialised = true;
}


void G4LivermorePolarizedComptonModel::InitialiseLocal(const G4ParticleDefinition*,
					      G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedComptonModel::ReadData(size_t Z, const char* path)
{
  if (verboseLevel > 1) 
    {
      G4cout << "G4LivermorePolarizedComptonModel::ReadData()" 
	     << G4endl;
    }
  if(data[Z]) { return; }  
  const char* datadir = path;
  if(!datadir) 
    {
      datadir = getenv("G4LEDATA");
      if(!datadir) 
	{
	  G4Exception("G4LivermorePolarizedComptonModel::ReadData()",
		      "em0006",FatalException,
		      "Environment variable G4LEDATA not defined");
	  return;
	}
    }
  
  data[Z] = new G4LPhysicsFreeVector();
  
  // Activation of spline interpolation
  data[Z]->SetSpline(false);
  
  std::ostringstream ost;
  ost << datadir << "/livermore/comp/ce-cs-" << Z <<".dat";
  std::ifstream fin(ost.str().c_str());
  
  if( !fin.is_open()) 
    {
      G4ExceptionDescription ed;
      ed << "G4LivermorePolarizedComptonModel data file <" << ost.str().c_str()
	 << "> is not opened!" << G4endl;
      G4Exception("G4LivermoreComptonModel::ReadData()",
		  "em0003",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW6.34 or later");
      return;
    } else {
    if(verboseLevel > 3) {
      G4cout << "File " << ost.str() 
	     << " is opened by G4LivermorePolarizedComptonModel" << G4endl;
    }   
    data[Z]->Retrieve(fin, true);
    data[Z]->ScaleVector(MeV, MeV*barn);
  }   
  fin.close();
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedComptonModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4LivermorePolarizedComptonModel" << G4endl;

  G4double cs = 0.0; 
  
  if (GammaEnergy < LowEnergyLimit()) 
    return 0.0;

  G4int intZ = G4lrint(Z);
  if(intZ < 1 || intZ > maxZ) { return cs; } 
  
  G4LPhysicsFreeVector* pv = data[intZ];
  
  // if element was not initialised
  // do initialisation safely for MT mode
  if(!pv) 
    {
      InitialiseForElement(0, intZ);
      pv = data[intZ];
      if(!pv) { return cs; }
    }
  
  G4int n = pv->GetVectorLength() - 1;   
  G4double e1 = pv->Energy(0);
  G4double e2 = pv->Energy(n);
  
  if(GammaEnergy <= e1)      { cs = GammaEnergy/(e1*e1)*pv->Value(e1); }
  else if(GammaEnergy <= e2) { cs = pv->Value(GammaEnergy)/GammaEnergy; }
  else if(GammaEnergy > e2)  { cs = pv->Value(e2)/GammaEnergy; }
  
  return cs;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedComptonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					      const G4MaterialCutsCouple* couple,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{
  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used (Nuc Phys 20(1960),15).
  // GEANT4 internal units
  //
  // Note : Effects due to binding of atomic electrons are negliged.

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4LivermorePolarizedComptonModel" << G4endl;

  G4double gammaEnergy0 = aDynamicGamma->GetKineticEnergy();
 
  // do nothing below the threshold
  // should never get here because the XS is zero below the limit
  if (gammaEnergy0 < LowEnergyLimit())     
    return ; 


  G4ThreeVector gammaPolarization0 = aDynamicGamma->GetPolarization();

  // Protection: a polarisation parallel to the
  // direction causes problems;
  // in that case find a random polarization

  G4ThreeVector gammaDirection0 = aDynamicGamma->GetMomentumDirection();

  // Make sure that the polarization vector is perpendicular to the
  // gamma direction. If not

  if(!(gammaPolarization0.isOrthogonal(gammaDirection0, 1e-6))||(gammaPolarization0.mag()==0))
    { // only for testing now
      gammaPolarization0 = GetRandomPolarization(gammaDirection0);
    }
  else
    {
      if ( gammaPolarization0.howOrthogonal(gammaDirection0) != 0)
	{
	  gammaPolarization0 = GetPerpendicularPolarization(gammaDirection0, gammaPolarization0);
	}
    }

  // End of Protection

  G4double E0_m = gammaEnergy0 / electron_mass_c2 ;

  // Select randomly one element in the current material
  //G4int Z = crossSectionHandler->SelectRandomAtom(couple,gammaEnergy0);
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,gammaEnergy0);
  G4int Z = (G4int)elm->GetZ();

  // Sample the energy and the polarization of the scattered photon

  G4double epsilon, epsilonSq, onecost, sinThetaSqr, greject ;

  G4double epsilon0Local = 1./(1. + 2*E0_m);
  G4double epsilon0Sq = epsilon0Local*epsilon0Local;
  G4double alpha1   = - std::log(epsilon0Local);
  G4double alpha2 = 0.5*(1.- epsilon0Sq);

  G4double wlGamma = h_Planck*c_light/gammaEnergy0;
  G4double gammaEnergy1;
  G4ThreeVector gammaDirection1;

  do {
    if ( alpha1/(alpha1+alpha2) > G4UniformRand() )
      {
	epsilon   = G4Exp(-alpha1*G4UniformRand());  
	epsilonSq = epsilon*epsilon; 
      }
    else 
      {
	epsilonSq = epsilon0Sq + (1.- epsilon0Sq)*G4UniformRand();
	epsilon   = std::sqrt(epsilonSq);
      }

    onecost = (1.- epsilon)/(epsilon*E0_m);
    sinThetaSqr   = onecost*(2.-onecost);

    // Protection
    if (sinThetaSqr > 1.)
      {
	G4cout
	  << " -- Warning -- G4LivermorePolarizedComptonModel::SampleSecondaries "
	  << "sin(theta)**2 = "
	  << sinThetaSqr
	  << "; set to 1"
	  << G4endl;
	sinThetaSqr = 1.;
      }
    if (sinThetaSqr < 0.)
      {
	G4cout
	  << " -- Warning -- G4LivermorePolarizedComptonModel::SampleSecondaries "
	  << "sin(theta)**2 = "
	  << sinThetaSqr
	  << "; set to 0"
	  << G4endl;
	sinThetaSqr = 0.;
      }
    // End protection

    G4double x =  std::sqrt(onecost/2.) / (wlGamma/cm);;
    G4double scatteringFunction = scatterFunctionData->FindValue(x,Z-1);
    greject = (1. - epsilon*sinThetaSqr/(1.+ epsilonSq))*scatteringFunction;

  } while(greject < G4UniformRand()*Z);


  // ****************************************************
  //		Phi determination
  // ****************************************************

  G4double phi = SetPhi(epsilon,sinThetaSqr);

  //
  // scattered gamma angles. ( Z - axis along the parent gamma)
  //

  G4double cosTheta = 1. - onecost;

  // Protection

  if (cosTheta > 1.)
    {
      G4cout
	<< " -- Warning -- G4LivermorePolarizedComptonModel::SampleSecondaries "
	<< "cosTheta = "
	<< cosTheta
	<< "; set to 1"
	<< G4endl;
      cosTheta = 1.;
    }
  if (cosTheta < -1.)
    {
      G4cout 
	<< " -- Warning -- G4LivermorePolarizedComptonModel::SampleSecondaries "
	<< "cosTheta = " 
	<< cosTheta
	<< "; set to -1"
	<< G4endl;
      cosTheta = -1.;
    }
  // End protection      
  
  
  G4double sinTheta = std::sqrt (sinThetaSqr);
  
  // Protection
  if (sinTheta > 1.)
    {
      G4cout 
	<< " -- Warning -- G4LivermorePolarizedComptonModel::SampleSecondaries "
	<< "sinTheta = " 
	<< sinTheta
	<< "; set to 1"
	<< G4endl;
      sinTheta = 1.;
    }
  if (sinTheta < -1.)
    {
      G4cout 
	<< " -- Warning -- G4LivermorePolarizedComptonModel::SampleSecondaries "
	<< "sinTheta = " 
	<< sinTheta
	<< "; set to -1" 
	<< G4endl;
      sinTheta = -1.;
    }
  // End protection
  
      
  G4double dirx = sinTheta*std::cos(phi);
  G4double diry = sinTheta*std::sin(phi);
  G4double dirz = cosTheta ;
  

  // oneCosT , eom

  // Doppler broadening -  Method based on:
  // Y. Namito, S. Ban and H. Hirayama, 
  // "Implementation of the Doppler Broadening of a Compton-Scattered Photon Into the EGS4 Code" 
  // NIM A 349, pp. 489-494, 1994
  
  // Maximum number of sampling iterations

  static G4int maxDopplerIterations = 1000;
  G4double bindingE = 0.;
  G4double photonEoriginal = epsilon * gammaEnergy0;
  G4double photonE = -1.;
  G4int iteration = 0;
  G4double eMax = gammaEnergy0;

  G4int shellIdx = 0;

  if (verboseLevel > 3) {
    G4cout << "Started loop to sample broading" << G4endl;
  }

  do
    {
      iteration++;
      // Select shell based on shell occupancy
      shellIdx = shellData->SelectRandomShell(Z);
      bindingE = shellData->BindingEnergy(Z,shellIdx);
      
      if (verboseLevel > 3) {
	G4cout << "Shell ID= " << shellIdx 
	       << " Ebind(keV)= " << bindingE/keV << G4endl;
      }
      eMax = gammaEnergy0 - bindingE;
      
      // Randomly sample bound electron momentum (memento: the data set is in Atomic Units)
      G4double pSample = profileData->RandomSelectMomentum(Z,shellIdx);

      if (verboseLevel > 3) {
       G4cout << "pSample= " << pSample << G4endl;
     }
      // Rescale from atomic units
      G4double pDoppler = pSample * fine_structure_const;
      G4double pDoppler2 = pDoppler * pDoppler;
      G4double var2 = 1. + onecost * E0_m;
      G4double var3 = var2*var2 - pDoppler2;
      G4double var4 = var2 - pDoppler2 * cosTheta;
      G4double var = var4*var4 - var3 + pDoppler2 * var3;
      if (var > 0.)
	{
	  G4double varSqrt = std::sqrt(var);        
	  G4double scale = gammaEnergy0 / var3;  
          // Random select either root
 	  if (G4UniformRand() < 0.5) photonE = (var4 - varSqrt) * scale;               
	  else photonE = (var4 + varSqrt) * scale;
	} 
      else
	{
	  photonE = -1.;
	}
   } while ( iteration <= maxDopplerIterations && 
	     (photonE < 0. || photonE > eMax || photonE < eMax*G4UniformRand()) );
  //while (iteration <= maxDopplerIterations && photonE > eMax); ???


  // End of recalculation of photon energy with Doppler broadening
  // Revert to original if maximum number of iterations threshold has been reached
  if (iteration >= maxDopplerIterations)
    {
      photonE = photonEoriginal;
      bindingE = 0.;
    }

  gammaEnergy1 = photonE;
 
  //
  // update G4VParticleChange for the scattered photon 
  //



  //  gammaEnergy1 = epsilon*gammaEnergy0;


  // New polarization

  G4ThreeVector gammaPolarization1 = SetNewPolarization(epsilon,
							sinThetaSqr,
							phi,
							cosTheta);

  // Set new direction
  G4ThreeVector tmpDirection1( dirx,diry,dirz );
  gammaDirection1 = tmpDirection1;

  // Change reference frame.

  SystemOfRefChange(gammaDirection0,gammaDirection1,
		    gammaPolarization0,gammaPolarization1);

  if (gammaEnergy1 > 0.)
    {
      fParticleChange->SetProposedKineticEnergy( gammaEnergy1 ) ;
      fParticleChange->ProposeMomentumDirection( gammaDirection1 );
      fParticleChange->ProposePolarization( gammaPolarization1 );
    }
  else
    {
      gammaEnergy1 = 0.;
      fParticleChange->SetProposedKineticEnergy(0.) ;
      fParticleChange->ProposeTrackStatus(fStopAndKill);
    }

  //
  // kinematic of the scattered electron
  //

  G4double ElecKineEnergy = gammaEnergy0 - gammaEnergy1 -bindingE;

  // SI -protection against negative final energy: no e- is created
  // like in G4LivermoreComptonModel.cc
  if(ElecKineEnergy < 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(gammaEnergy0 - gammaEnergy1);
    return;
  }
 
  // SI - Removed range test
  
  G4double ElecMomentum = std::sqrt(ElecKineEnergy*(ElecKineEnergy+2.*electron_mass_c2));

  G4ThreeVector ElecDirection((gammaEnergy0 * gammaDirection0 -
				   gammaEnergy1 * gammaDirection1) * (1./ElecMomentum));

  G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),ElecDirection.unit(),ElecKineEnergy) ;
  fvect->push_back(dp);

  // sample deexcitation
  //
  
  if (verboseLevel > 3) {
    G4cout << "Started atomic de-excitation " << fAtomDeexcitation << G4endl;
  }
  
  if(fAtomDeexcitation && iteration < maxDopplerIterations) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      size_t nbefore = fvect->size();
      G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIdx);
      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      size_t nafter = fvect->size();
      if(nafter > nbefore) {
	for (size_t i=nbefore; i<nafter; ++i) {
	  //Check if there is enough residual energy 
	  if (bindingE >= ((*fvect)[i])->GetKineticEnergy())
            {
              //Ok, this is a valid secondary: keep it
	      bindingE -= ((*fvect)[i])->GetKineticEnergy();
            }
	  else
            {
	      //Invalid secondary: not enough energy to create it!
	      //Keep its energy in the local deposit
	      delete (*fvect)[i]; 
              (*fvect)[i]=0;
            }
	} 
      }
    }
  }
  //This should never happen
  if(bindingE < 0.0) 
    G4Exception("G4LivermoreComptonModel::SampleSecondaries()",
		"em2050",FatalException,"Negative local energy deposit");   
  
  fParticleChange->ProposeLocalEnergyDeposit(bindingE);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedComptonModel::SetPhi(G4double energyRate,
					     G4double sinSqrTh)
{
  G4double rand1;
  G4double rand2;
  G4double phiProbability;
  G4double phi;
  G4double a, b;

  do
    {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      phiProbability=0.;
      phi = twopi*rand1;
      
      a = 2*sinSqrTh;
      b = energyRate + 1/energyRate;
      
      phiProbability = 1 - (a/b)*(std::cos(phi)*std::cos(phi));

      
 
    }
  while ( rand2 > phiProbability );
  return phi;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LivermorePolarizedComptonModel::SetPerpendicularVector(G4ThreeVector& a)
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

G4ThreeVector G4LivermorePolarizedComptonModel::GetRandomPolarization(G4ThreeVector& direction0)
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

G4ThreeVector G4LivermorePolarizedComptonModel::GetPerpendicularPolarization
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

G4ThreeVector G4LivermorePolarizedComptonModel::SetNewPolarization(G4double epsilon,
							      G4double sinSqrTh, 
							      G4double phi,
							      G4double costheta) 
{
  G4double rand1;
  G4double rand2;
  G4double cosPhi = std::cos(phi);
  G4double sinPhi = std::sin(phi);
  G4double sinTheta = std::sqrt(sinSqrTh);
  G4double cosSqrPhi = cosPhi*cosPhi;
  //  G4double cossqrth = 1.-sinSqrTh;
  //  G4double sinsqrphi = sinPhi*sinPhi;
  G4double normalisation = std::sqrt(1. - cosSqrPhi*sinSqrTh);
 

  // Determination of Theta 
  
  // ---- MGP ---- Commented out the following 3 lines to avoid compilation 
  // warnings (unused variables)
  // G4double thetaProbability;
  G4double theta;
  // G4double a, b;
  // G4double cosTheta;

  /*

  depaola method
  
  do
  {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      thetaProbability=0.;
      theta = twopi*rand1;
      a = 4*normalisation*normalisation;
      b = (epsilon + 1/epsilon) - 2;
      thetaProbability = (b + a*std::cos(theta)*std::cos(theta))/(a+b);
      cosTheta = std::cos(theta);
    }
  while ( rand2 > thetaProbability );
  
  G4double cosBeta = cosTheta;

  */


  // Dan Xu method (IEEE TNS, 52, 1160 (2005))

  rand1 = G4UniformRand();
  rand2 = G4UniformRand();

  if (rand1<(epsilon+1.0/epsilon-2)/(2.0*(epsilon+1.0/epsilon)-4.0*sinSqrTh*cosSqrPhi))
    {
      if (rand2<0.5)
	theta = pi/2.0;
      else
	theta = 3.0*pi/2.0;
    }
  else
    {
      if (rand2<0.5)
	theta = 0;
      else
	theta = pi;
    }
  G4double cosBeta = std::cos(theta);
  G4double sinBeta = std::sqrt(1-cosBeta*cosBeta);
  
  G4ThreeVector gammaPolarization1;

  G4double xParallel = normalisation*cosBeta;
  G4double yParallel = -(sinSqrTh*cosPhi*sinPhi)*cosBeta/normalisation;
  G4double zParallel = -(costheta*sinTheta*cosPhi)*cosBeta/normalisation;
  G4double xPerpendicular = 0.;
  G4double yPerpendicular = (costheta)*sinBeta/normalisation;
  G4double zPerpendicular = -(sinTheta*sinPhi)*sinBeta/normalisation;

  G4double xTotal = (xParallel + xPerpendicular);
  G4double yTotal = (yParallel + yPerpendicular);
  G4double zTotal = (zParallel + zPerpendicular);
  
  gammaPolarization1.setX(xTotal);
  gammaPolarization1.setY(yTotal);
  gammaPolarization1.setZ(zTotal);
  
  return gammaPolarization1;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedComptonModel::SystemOfRefChange(G4ThreeVector& direction0,
						    G4ThreeVector& direction1,
						    G4ThreeVector& polarization0,
						    G4ThreeVector& polarization1)
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
  
  direction1 = (direction_x*Axis_X0 + direction_y*Axis_Y0 + direction_z*Axis_Z0).unit();
  G4double polarization_x = polarization1.getX();
  G4double polarization_y = polarization1.getY();
  G4double polarization_z = polarization1.getZ();

  polarization1 = (polarization_x*Axis_X0 + polarization_y*Axis_Y0 + polarization_z*Axis_Z0).unit();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AutoLock.hh"
namespace { G4Mutex LivermorePolarizedComptonModelMutex = G4MUTEX_INITIALIZER; }

void 
G4LivermorePolarizedComptonModel::InitialiseForElement(const G4ParticleDefinition*, 
					      G4int Z)
{
  G4AutoLock l(&LivermorePolarizedComptonModelMutex);
  //  G4cout << "G4LivermoreComptonModel::InitialiseForElement Z= " 
  //   << Z << G4endl;
  if(!data[Z]) { ReadData(Z); }
  l.unlock();
}
