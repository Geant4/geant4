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

#include "G4LivermorePolarizedPhotoElectricModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedPhotoElectricModel::G4LivermorePolarizedPhotoElectricModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),isInitialised(false),meanFreePathTable(0),crossSectionHandler(0), shellCrossSectionHandler(0)
{
  lowEnergyLimit = 250 * eV; // SI - Could be 10 eV ?
  highEnergyLimit = 100 * GeV;
  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);
  
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  SetDeexcitationFlag(true);
  ActivateAuger(false);

  G4cout << "Livermore Polarized PhotoElectric is constructed " << G4endl
         << "Energy range: "
         << lowEnergyLimit / keV << " keV - "
         << highEnergyLimit / GeV << " GeV"
         << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedPhotoElectricModel::~G4LivermorePolarizedPhotoElectricModel()
{  
  //  if (meanFreePathTable)   delete meanFreePathTable;
  if (crossSectionHandler) delete crossSectionHandler;
  if (shellCrossSectionHandler) delete shellCrossSectionHandler;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedPhotoElectricModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& cuts)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4LivermorePolarizedPhotoElectricModel::Initialise()" << G4endl;

  if (crossSectionHandler)
  {
    crossSectionHandler->Clear();
    delete crossSectionHandler;
  }

  if (shellCrossSectionHandler)
    {
      shellCrossSectionHandler->Clear();
      delete shellCrossSectionHandler;
    }


  /*
  // Energy limits
  
  if (LowEnergyLimit() < lowEnergyLimit)
  {
  G4cout << "G4LivermorePolarizedPhotoElectricModel: low energy limit increased from " << 
  LowEnergyLimit()/eV << " eV to " << lowEnergyLimit << " eV" << G4endl;
  SetLowEnergyLimit(lowEnergyLimit);
  }
  
  if (HighEnergyLimit() > highEnergyLimit)
  {
  G4cout << "G4LivermorePolarizedPhotoElectricModel: high energy limit decreased from " << 
  HighEnergyLimit()/GeV << " GeV to " << highEnergyLimit << " GeV" << G4endl;
  SetHighEnergyLimit(highEnergyLimit);
  }
  */
  
  // Reading of data files - all materials are read
  
  crossSectionHandler = new G4CrossSectionHandler;
  crossSectionHandler->Clear();
  G4String crossSectionFile = "phot/pe-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  meanFreePathTable = 0;
  //  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();

  shellCrossSectionHandler = new G4CrossSectionHandler();
  shellCrossSectionHandler->Clear();
  G4String shellCrossSectionFile = "phot/pe-ss-cs-";
  shellCrossSectionHandler->LoadShellData(shellCrossSectionFile);


  //
  if (verboseLevel > 2) 
    G4cout << "Loaded cross section files for Livermore Polarized PhotoElectric model" << G4endl;
  
  InitialiseElementSelectors(particle,cuts);

  G4cout << "Livermore Polarized PhotoElectric model is initialized " << G4endl
         << "Energy range: "
         << LowEnergyLimit() / keV << " keV - "
         << HighEnergyLimit() / GeV << " GeV"
         << G4endl;

  //
    
  if(isInitialised) return;

  /*  if(pParticleChange)
    fParticleChange = reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
  else
    fParticleChange = new G4ParticleChangeForGamma();
  */
    
  fParticleChange = GetParticleChangeForGamma();
  
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedPhotoElectricModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4LivermorePolarizedPhotoElectricModel" << G4endl;

  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedPhotoElectricModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
							       const G4MaterialCutsCouple* couple,
							       const G4DynamicParticle* aDynamicGamma,
							       G4double,
							       G4double)
{

  // Fluorescence generated according to:
  // J. Stepanek ,"A program to determine the radiation spectra due to a single atomic
  // subshell ionisation by a particle or due to deexcitation or decay of radionuclides",
  // Comp. Phys. Comm. 1206 pp 1-1-9 (1997)

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4LivermorePolarizedPhotoElectricModel" << G4endl;

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();
  G4ThreeVector gammaPolarization0 = aDynamicGamma->GetPolarization();  
  G4ThreeVector photonDirection = aDynamicGamma->GetMomentumDirection();
  
  // kill incident photon
  
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  
  // low-energy gamma is absorpted by this process

  if (photonEnergy <= lowEnergyLimit)
    {
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
      return;
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

  // Select randomly one element in the current material

  //  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy);

  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple->GetMaterial(),particle,photonEnergy);
  G4int Z = (G4int)elm->GetZ();

  // Select the ionised shell in the current atom according to shell cross sections

  size_t shellIndex = shellCrossSectionHandler->SelectRandomShell(Z,photonEnergy);

  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  const G4AtomicShell* shell = transitionManager->Shell(Z,shellIndex);
  G4double bindingEnergy = shell->BindingEnergy();
  G4int shellId = shell->ShellId();
  
  // Primary outgoing  electron
  
  G4double eKineticEnergy = photonEnergy - bindingEnergy;


  if (eKineticEnergy > 0.)
    { 

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
    }
  else
    {
      bindingEnergy = photonEnergy;
    }
  

  // deexcitation 
  if(DeexcitationFlag() && Z > 5) {
    const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();
    size_t index = couple->GetIndex();
    G4double cutg = (*(theCoupleTable->GetEnergyCutsVector(0)))[index];
    //cutg = std::min(cutForLowEnergySecondaryPhotons,cutg);
    G4double cute = (*(theCoupleTable->GetEnergyCutsVector(1)))[index];
    //cute = std::min(cutForLowEnergySecondaryPhotons,cute);
    
    //   G4DynamicParticle* aPhoton;
    
    // Generation of fluorescence
    // Data in EADL are available only for Z > 5
    // Protection to avoid generating photons in the unphysical case of
    // shell binding energy > photon energy
    if (bindingEnergy > cutg || bindingEnergy > cute)
      {
	G4DynamicParticle* aPhoton;
	deexcitationManager.SetCutForSecondaryPhotons(cutg);
	deexcitationManager.SetCutForAugerElectrons(cute);
	std::vector<G4DynamicParticle*>* photonVector = 
	  deexcitationManager.GenerateParticles(Z,shellId);
	size_t nTotPhotons = photonVector->size();
	for (size_t k=0; k<nTotPhotons; k++)
	  {
	    aPhoton = (*photonVector)[k];
	    if (aPhoton)
	      {
		G4double itsEnergy = aPhoton->GetKineticEnergy();
		if (itsEnergy <= bindingEnergy)
                {
                  // Local energy deposit is given as the sum of the
                  // energies of incident photons minus the energies
                  // of the outcoming fluorescence photons
                  bindingEnergy -= itsEnergy;
		  fvect->push_back(aPhoton);
                }
		else
		  {
		    delete aPhoton;
		    (*photonVector)[k] = 0;
		  }
	      }
	  }
	delete photonVector;
      }
  }
  // excitation energy left
  fParticleChange->ProposeLocalEnergyDeposit(bindingEnergy);
}

/*
  energyDeposit += bindingEnergy;
  // Final state
  
  for (G4int l = 0; l<nElectrons; l++ )
  {
  aPhoton = electronVector[l];
  if(aPhoton) {
  fvect->push_back(aPhoton);
  }
  }
  for ( size_t ll = 0; ll < nTotPhotons; ll++)
  {
  aPhoton = (*photonVector)[ll];
  if(aPhoton) {
  fvect->push_back(aPhoton);
  }
    }
    
    delete photonVector;
    
    if (energyDeposit < 0)
    {
    G4cout << "WARNING - "
    << "G4LowEnergyPhotoElectric::PostStepDoIt - Negative energy deposit"
    << G4endl;
    energyDeposit = 0;
    }
    
    // kill incident photon
  fParticleChange->ProposeMomentumDirection( 0., 0., 0. );
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  fParticleChange->ProposeLocalEnergyDeposit(energyDeposit);

}
*/


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedPhotoElectricModel::ActivateAuger(G4bool augerbool)
{
  if (!DeexcitationFlag() && augerbool)
    {
      G4cout << "WARNING - G4LivermorePolarizedPhotoElectricModel" << G4endl;
      G4cout << "The use of the Atomic Deexcitation Manager is set to false " << G4endl;
      G4cout << "Therefore, Auger electrons will be not generated anyway" << G4endl;
    }
  deexcitationManager.ActivateAugerElectronProduction(augerbool);
  if (verboseLevel > 1)
    G4cout << "Auger production set to " << augerbool << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4LivermorePolarizedPhotoElectricModel::SetCosTheta(G4double energyE)

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



G4double G4LivermorePolarizedPhotoElectricModel::SetPhi(G4double Ph_energy,
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

G4ThreeVector G4LivermorePolarizedPhotoElectricModel::SetPerpendicularVector(G4ThreeVector& a)
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

G4ThreeVector G4LivermorePolarizedPhotoElectricModel::GetRandomPolarization(G4ThreeVector& direction0)
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

G4ThreeVector G4LivermorePolarizedPhotoElectricModel::GetPerpendicularPolarization
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


void G4LivermorePolarizedPhotoElectricModel::SystemOfRefChange
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
   G4double G4LivermorePolarizedPhotoElectricModel::GetMeanFreePath(const G4Track& track,
   G4double,
   G4ForceCondition*)
   {
   
   const G4DynamicParticle* photon = track.GetDynamicParticle();
   G4double energy = photon->GetKineticEnergy();
   G4Material* material = track.GetMaterial();
   //  size_t materialIndex = material->GetIndex();
   
   G4double meanFreePath = DBL_MAX;
   
   //  if (energy > highEnergyLimit)
   //    meanFreePath = meanFreePathTable->FindValue(highEnergyLimit,materialIndex);
   //  else if (energy < lowEnergyLimit) meanFreePath = DBL_MAX;
   //  else meanFreePath = meanFreePathTable->FindValue(energy,materialIndex);
   
   G4double cross = shellCrossSectionHandler->ValueForMaterial(material,energy);
   if(cross > 0.0) meanFreePath = 1.0/cross;
   
   return meanFreePath;
   
   
   }
*/


