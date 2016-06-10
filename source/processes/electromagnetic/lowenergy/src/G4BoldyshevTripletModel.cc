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
// $Id: G4BoldyshevTripletModel.cc 74822 2013-10-22 14:42:13Z gcosmo $
// GEANT4 tag $Name:  $
//
//
// Author: Gerardo Depaola & Francesco Longo
//
// History:
// --------
//   23-06-2010 First implementation as model 


#include "G4BoldyshevTripletModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BoldyshevTripletModel::G4BoldyshevTripletModel(const G4ParticleDefinition*,
								 const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),smallEnergy(4.*MeV),isInitialised(false),
   crossSectionHandler(0),meanFreePathTable(0)
{
  lowEnergyLimit = 4.0*electron_mass_c2;
  highEnergyLimit = 100 * GeV;
  SetHighEnergyLimit(highEnergyLimit);
  	 
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(verboseLevel > 0) {
    G4cout << "Triplet Gamma conversion is constructed " << G4endl
	   << "Energy range: "
	   << lowEnergyLimit / MeV << " MeV - "
	   << highEnergyLimit / GeV << " GeV"
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BoldyshevTripletModel::~G4BoldyshevTripletModel()
{  
  if (crossSectionHandler) delete crossSectionHandler;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4BoldyshevTripletModel::Initialise(const G4ParticleDefinition*,
					    const G4DataVector&)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4BoldyshevTripletModel::Initialise()" << G4endl;

  if (crossSectionHandler)
  {
    crossSectionHandler->Clear();
    delete crossSectionHandler;
  }

  // Read data tables for all materials
  
  crossSectionHandler = new G4CrossSectionHandler();
  crossSectionHandler->Initialise(0,lowEnergyLimit,100.*GeV,400);
  G4String crossSectionFile = "tripdata/pp-trip-cs-"; // here only pair in electron field cs should be used
  crossSectionHandler->LoadData(crossSectionFile);

  //
  
  if (verboseLevel > 0) {
    G4cout << "Loaded cross section files for Livermore GammaConversion" << G4endl;
    G4cout << "To obtain the total cross section this should be used only " << G4endl 
	   << "in connection with G4NuclearGammaConversion " << G4endl;
  }

  if (verboseLevel > 0) { 
    G4cout << "Livermore Electron Gamma Conversion model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / MeV << " MeV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4BoldyshevTripletModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
							    G4double GammaEnergy,
							    G4double Z, G4double,
							    G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4BoldyshevTripletModel" 
	   << G4endl;
  }
  if (GammaEnergy < lowEnergyLimit || GammaEnergy > highEnergyLimit) return 0;

  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BoldyshevTripletModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					      const G4MaterialCutsCouple* ,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{

// The energies of the secondary particles are sampled using
// a modified Wheeler-Lamb model (see PhysRevD 7 (1973), 26) 

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4BoldyshevTripletModel" << G4endl;

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();
  G4ParticleMomentum photonDirection = aDynamicGamma->GetMomentumDirection();

  G4double epsilon ;
  G4double p0 = electron_mass_c2; 
  
  G4double positronTotEnergy, electronTotEnergy, thetaEle, thetaPos;
  G4double ener_re=0., theta_re, phi_re, phi;

  // Calculo de theta - elecron de recoil
  
  G4double energyThreshold = sqrt(2.)*electron_mass_c2; // -> momentumThreshold_N = 1
  energyThreshold = 1.1*electron_mass_c2;
  // G4cout << energyThreshold << G4endl; 
 
  G4double momentumThreshold_c = sqrt(energyThreshold * energyThreshold - electron_mass_c2*electron_mass_c2); // momentun in MeV/c unit
  G4double momentumThreshold_N = momentumThreshold_c/electron_mass_c2; // momentun in mc unit
  
  // Calculation of recoil electron production 
  
  G4double SigmaTot = (28./9.) * std::log ( 2.* photonEnergy / electron_mass_c2 ) - 218. / 27. ; 
  G4double X_0 = 2. * ( sqrt(momentumThreshold_N*momentumThreshold_N + 1) -1 );
  G4double SigmaQ = (82./27. - (14./9.) * log (X_0) + 4./15.*X_0 - 0.0348 * X_0 * X_0); 
  G4double recoilProb = G4UniformRand();
  //G4cout << "SIGMA TOT " << SigmaTot <<  " " << "SigmaQ " << SigmaQ << " " << SigmaQ/SigmaTot << " " << recoilProb << G4endl;
  
  if (recoilProb >= SigmaQ/SigmaTot) // create electron recoil 
    { 
      
      G4double cosThetaMax = (  ( energyThreshold - electron_mass_c2 ) / (momentumThreshold_c) + electron_mass_c2*
	    			( energyThreshold + electron_mass_c2 ) / (photonEnergy*momentumThreshold_c) );
      
      if (cosThetaMax > 1) G4cout << "ERRORE " << G4endl;
      
      G4double r1;
      G4double r2;
      G4double are, bre, loga, f1_re, greject, cost;
      
      do {
	r1 = G4UniformRand();
	r2 = G4UniformRand();
	//	cost = (pow(4./enern,0.5*r1)) ;
	cost = pow(cosThetaMax,r1);
	theta_re = acos(cost);
	are = 1./(14.*cost*cost);
	bre = (1.-5.*cost*cost)/(2.*cost);
	loga = log((1.+ cost)/(1.- cost));
	f1_re = 1. - bre*loga;
	
	if ( theta_re >= 4.47*CLHEP::pi/180.)
	  {
	    greject = are*f1_re;
	  } else {
	  greject = 1. ;
	}
      } while(greject < r2);
      
      // Calculo de phi - elecron de recoil
      
      G4double r3, r4, rt;
      
      do {
	
	r3 = G4UniformRand();
	r4 = G4UniformRand();
	phi_re = twopi*r3 ;
	G4double sint2 = 1. - cost*cost ;
	G4double fp = 1. - sint2*loga/(2.*cost) ;
	rt = (1.-cos(2.*phi_re)*fp/f1_re)/(2.*pi) ;
	
      } while(rt < r4);
      
      // Calculo de la energia - elecron de recoil - relacion momento maximo <-> angulo
      
      G4double S = electron_mass_c2*(2.* photonEnergy + electron_mass_c2);
      G4double D2 = 4.*S * electron_mass_c2*electron_mass_c2
	+ (S - electron_mass_c2*electron_mass_c2)
	*(S - electron_mass_c2*electron_mass_c2)*sin(theta_re)*sin(theta_re);
      ener_re = electron_mass_c2 * (S + electron_mass_c2*electron_mass_c2)/sqrt(D2);
      
      // New Recoil energy calculation 

      G4double momentum_recoil = 2* (electron_mass_c2) * (std::cos(theta_re)/(std::sin(phi_re)*std::sin(phi_re)));
      G4double ener_recoil = sqrt( momentum_recoil*momentum_recoil + electron_mass_c2*electron_mass_c2);
      ener_re = ener_recoil;

      //      G4cout << "electron de retroceso " << ener_re << " " << theta_re << " " << phi_re << G4endl;
      
      // Recoil electron creation
      G4double dxEle_re=sin(theta_re)*std::cos(phi_re),dyEle_re=sin(theta_re)*std::sin(phi_re), dzEle_re=cos(theta_re);
      
      G4double electronRKineEnergy = std::max(0.,ener_re - electron_mass_c2) ;
      
      G4ThreeVector electronRDirection (dxEle_re, dyEle_re, dzEle_re);
      electronRDirection.rotateUz(photonDirection);
      
      G4DynamicParticle* particle3 = new G4DynamicParticle (G4Electron::Electron(),
							    electronRDirection,
							    electronRKineEnergy);
      fvect->push_back(particle3);	    
      
    }
  else
    {
      // deposito la energia  ener_re - electron_mass_c2
      // G4cout << "electron de retroceso " << ener_re << G4endl;
      fParticleChange->ProposeLocalEnergyDeposit(ener_re - electron_mass_c2);
    }
  
  // Depaola (2004) suggested distribution for e+e- energy
  
  //  G4double t = 0.5*asinh(momentumThreshold_N);
  G4double t = 0.5*log(momentumThreshold_N + sqrt(momentumThreshold_N*momentumThreshold_N+1));
  
  G4cout << 0.5*asinh(momentumThreshold_N) << "  " << t << G4endl;

  G4double J1 = 0.5*(t*cosh(t)/sinh(t) - log(2.*sinh(t)));
  G4double J2 = (-2./3.)*log(2.*sinh(t)) + t*cosh(t)/sinh(t) + (sinh(t)-t*pow(cosh(t),3))/(3.*pow(sinh(t),3));
  G4double b = 2.*(J1-J2)/J1;
  
  G4double n = 1 - b/6.;
  G4double re=0.;
  re = G4UniformRand();
  G4double a = 0.;
  G4double b1 =  16. - 3.*b - 36.*b*re*n + 36.*b*pow(re,2.)*pow(n,2.) + 
    6.*pow(b,2.)*re*n;
  a = pow((b1/b),0.5);
  G4double c1 = (-6. + 12.*re*n + b + 2*a)*pow(b,2.);
  epsilon = (pow(c1,1./3.))/(2.*b) + (b-4.)/(2.*pow(c1,1./3.))+0.5;
  
  G4double photonEnergy1 = photonEnergy - ener_re ; // resto al foton la energia del electron de retro.
  positronTotEnergy = epsilon*photonEnergy1;
  electronTotEnergy = photonEnergy1 - positronTotEnergy; // temporarly
  
  G4double momento_e = sqrt(electronTotEnergy*electronTotEnergy - 
			    electron_mass_c2*electron_mass_c2) ;
  G4double momento_p = sqrt(positronTotEnergy*positronTotEnergy - 
			    electron_mass_c2*electron_mass_c2) ;
  
  thetaEle = acos((sqrt(p0*p0/(momento_e*momento_e) +1.)- p0/momento_e)) ;
  thetaPos = acos((sqrt(p0*p0/(momento_p*momento_p) +1.)- p0/momento_p)) ;
  phi  = twopi * G4UniformRand();
  
  G4double dxEle= std::sin(thetaEle)*std::cos(phi),dyEle= std::sin(thetaEle)*std::sin(phi),dzEle=std::cos(thetaEle);
  G4double dxPos=-std::sin(thetaPos)*std::cos(phi),dyPos=-std::sin(thetaPos)*std::sin(phi),dzPos=std::cos(thetaPos);
  
  
  // Kinematics of the created pair:
  // the electron and positron are assumed to have a symetric angular 
  // distribution with respect to the Z axis along the parent photon
  
  G4double electronKineEnergy = std::max(0.,electronTotEnergy - electron_mass_c2) ;
  
  // SI - The range test has been removed wrt original G4LowEnergyGammaconversion class
  
  G4ThreeVector electronDirection (dxEle, dyEle, dzEle);
  electronDirection.rotateUz(photonDirection);
  
  G4DynamicParticle* particle1 = new G4DynamicParticle (G4Electron::Electron(),
							electronDirection,
							electronKineEnergy);
  
  // The e+ is always created (even with kinetic energy = 0) for further annihilation
  G4double positronKineEnergy = std::max(0.,positronTotEnergy - electron_mass_c2) ;
  
  // SI - The range test has been removed wrt original G4LowEnergyGammaconversion class
  
  G4ThreeVector positronDirection (dxPos, dyPos, dzPos);
  positronDirection.rotateUz(photonDirection);   
  
  // Create G4DynamicParticle object for the particle2 
  G4DynamicParticle* particle2 = new G4DynamicParticle(G4Positron::Positron(),
						       positronDirection, positronKineEnergy);
  // Fill output vector
  
  
  fvect->push_back(particle1);
  fvect->push_back(particle2);
  

  
  
  // kill incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);   
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




