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
// $Id: G4LivermorePolarizedComptonModel.cc 82874 2014-07-15 15:25:29Z gcosmo $
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedComptonModel::G4LivermorePolarizedComptonModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),isInitialised(false),
   meanFreePathTable(0),scatterFunctionData(0),crossSectionHandler(0)
{ 
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if( verboseLevel>0 )  
    G4cout << "Livermore Polarized Compton is constructed " << G4endl;
        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedComptonModel::~G4LivermorePolarizedComptonModel()
{  
  if (meanFreePathTable)   delete meanFreePathTable;
  if (crossSectionHandler) delete crossSectionHandler;
  if (scatterFunctionData) delete scatterFunctionData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedComptonModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& cuts)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4LivermorePolarizedComptonModel::Initialise()" << G4endl;

  if (crossSectionHandler)
  {
    crossSectionHandler->Clear();
    delete crossSectionHandler;
  }

  // Reading of data files - all materials are read
  
  crossSectionHandler = new G4CrossSectionHandler;
  crossSectionHandler->Clear();
  G4String crossSectionFile = "comp/ce-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  meanFreePathTable = 0;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();

  G4VDataSetAlgorithm* scatterInterpolation = new G4LogLogInterpolation;
  G4String scatterFile = "comp/ce-sf-";
  scatterFunctionData = new G4CompositeEMDataSet(scatterInterpolation, 1., 1.);
  scatterFunctionData->LoadData(scatterFile);

  // For Doppler broadening
  shellData.SetOccupancyData();
  G4String file = "/doppler/shell-doppler";
  shellData.LoadData(file);

  if (verboseLevel > 2) 
    G4cout << "Loaded cross section files for Livermore Polarized Compton model" << G4endl;

  InitialiseElementSelectors(particle,cuts);

  if(  verboseLevel>0 ) { 
    G4cout << "Livermore Polarized Compton model is initialized " << G4endl
         << "Energy range: "
         << LowEnergyLimit() / eV << " eV - "
         << HighEnergyLimit() / GeV << " GeV"
         << G4endl;
  }
  
  //
    
  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
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

  if (GammaEnergy < LowEnergyLimit()) 
    return 0.0;

  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);
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
	epsilon   = std::exp(-alpha1*G4UniformRand());  
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

  G4int maxDopplerIterations = 1000;
  G4double bindingE = 0.;
  G4double photonEoriginal = epsilon * gammaEnergy0;
  G4double photonE = -1.;
  G4int iteration = 0;
  G4double eMax = gammaEnergy0;

  do
    {
      iteration++;
      // Select shell based on shell occupancy
      G4int shell = shellData.SelectRandomShell(Z);
      bindingE = shellData.BindingEnergy(Z,shell);
      
      eMax = gammaEnergy0 - bindingE;
     
      // Randomly sample bound electron momentum (memento: the data set is in Atomic Units)
      G4double pSample = profileData.RandomSelectMomentum(Z,shell);
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

  fParticleChange->ProposeLocalEnergyDeposit(bindingE);
  
  G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),ElecDirection.unit(),ElecKineEnergy) ;
  fvect->push_back(dp);

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


