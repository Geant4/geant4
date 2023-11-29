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
// File name:   G4eplusTo3GammaOKVIModel
//
// Authors:  Andrei Alkin, Vladimir Ivanchenko, Omrame Kadri
//
// Creation date: 29.03.2018
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eplusTo3GammaOKVIModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4TrackStatus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eplusTo3GammaOKVIModel::G4eplusTo3GammaOKVIModel(const G4ParticleDefinition*,
                                                   const G4String& nam)
  : G4VEmModel(nam), fDelta(0.001)
{
  theGamma = G4Gamma::Gamma();
  fParticleChange = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusTo3GammaOKVIModel::~G4eplusTo3GammaOKVIModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusTo3GammaOKVIModel::Initialise(const G4ParticleDefinition*,
					  const G4DataVector&)
{
  // here particle change is set for the triplet model
  if(fParticleChange) { return; }
  fParticleChange = GetParticleChangeForGamma();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// (A.A.) F_{ijk} calculation method
G4double G4eplusTo3GammaOKVIModel::ComputeF(G4double fr1, G4double fr2, 
                                            G4double fr3, G4double kinEnergy) 
{
  G4double ekin   = std::max(eV,kinEnergy);
  G4double tau    = ekin/electron_mass_c2;
  G4double gam    = tau + 1.0;
  G4double gamma2 = gam*gam;
  G4double bg2    = tau * (tau+2.0);
  G4double bg     = sqrt(bg2);

  G4double rho = (gamma2+4.*gam+1.)*G4Log(gam+bg)/(gamma2-1.) 
    - (gam+3.)/(sqrt(gam*gam - 1.)) + 1.;   
  G4double border;

  if(ekin < 500*MeV) { 
    border = 1. - (electron_mass_c2)/(2*(ekin + electron_mass_c2)); 
  } else { 
    border = 1. - (100*electron_mass_c2)/(2*(ekin + electron_mass_c2)); 
  }

  border = std::min(border, 0.9999);

  if (fr1>border)  { fr1 = border; }
  if (fr2>border)  { fr2 = border; }
  if (fr3>border)  { fr3 = border; }

  G4double  fr1s = fr1*fr1; // "s" for "squared" 
  G4double  fr2s = fr2*fr2;
  G4double  fr3s = fr3*fr3; 

  G4double  aa = (1.-fr1)*(1.-fr2);
  G4double  ab = fr3s + (fr1-fr2)*(fr1-fr2);
  G4double  add= ((1.-fr1)*(1.-fr1) + (1.-fr2)*(1.-fr2))/(fr3s*aa);

  G4double  fres = -rho*(1./fr1s + 1./fr2s) 
    + (ab/(2.*(fr1*fr2*aa)))*(G4Log(2.*gam*aa/(fr1*fr2))) 
    + (ab/(2.*fr1*fr2*(1-fr3)))*G4Log(2.*gam*(1.-fr3)/(fr1*fr2)) - add;

  return fres;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// (A.A.) F_{ijk} calculation method
G4double G4eplusTo3GammaOKVIModel::ComputeF0(G4double fr1, G4double fr2, 
                                             G4double fr3) 
{
  G4double tau    = 0.0;
  G4double gam    = tau + 1.0;
  G4double gamma2 = gam*gam;
  G4double bg2    = tau * (tau+2.0);
  G4double bg     = sqrt(bg2);

  G4double rho = (gamma2+4.*gam+1.)*G4Log(gam+bg)/(gamma2-1.) 
    - (gam+3.)/(sqrt(gam*gam - 1.)) + 1.;   
  G4double border = 0.5;

  if (fr1>border)  { fr1 = border; }
  if (fr2>border)  { fr2 = border; }
  if (fr3>border)  { fr3 = border; }

  G4double  fr1s = fr1*fr1; // "s" for "squared" 
  G4double  fr2s = fr2*fr2;
  G4double  fr3s = fr3*fr3; 

  G4double  aa = (1.-fr1)*(1.-fr2);
  G4double  ab = fr3s + (fr1-fr2)*(fr1-fr2);
  G4double  add= ((1.-fr1)*(1.-fr1) + (1.-fr2)*(1.-fr2))/(fr3s*aa);

  G4double  fres = -rho*(1./fr1s + 1./fr2s) 
    + (ab/(2.*(fr1*fr2*aa)))*(G4Log(2.*gam*aa/(fr1*fr2))) 
    + (ab/(2.*fr1*fr2*(1-fr3)))*G4Log(2.*gam*(1.-fr3)/(fr1*fr2)) - add;

  return fres;
}

//(A.A.) diff x-sections for maximum search and rejection
G4double G4eplusTo3GammaOKVIModel::ComputeFS(G4double fr1, 
         G4double fr2, G4double fr3, G4double kinEnergy) 
{
  G4double ekin  = std::max(eV,kinEnergy);   
  G4double tau   = ekin/electron_mass_c2;
  G4double gam   = tau + 1.0;  

  G4double fsum = fr1*fr1*(ComputeF(fr1,fr2,fr3,ekin) + 
			   ComputeF(fr3,fr1,fr2,ekin) + 
			   ComputeF(fr2,fr3,fr1,ekin));

  G4double dcross = fsum/((3*fr1*fr1*(gam+1.)));

  return dcross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4eplusTo3GammaOKVIModel::ComputeCrossSectionPerElectron(G4double kinEnergy)
{
  // Calculates the cross section per electron of annihilation into 3 photons
  // from the Heilter formula.

  G4double ekin   = std::max(eV,kinEnergy);   
  G4double tau    = ekin/electron_mass_c2;
  G4double gam    = tau + 1.0;
  G4double gamma2 = gam*gam;
  G4double bg2    = tau * (tau+2.0);
  G4double bg     = sqrt(bg2);
  
  G4double rho = (gamma2+4*gam+1.)*G4Log(gam+bg)/(gamma2-1.) 
    - (gam+3.)/(sqrt(gam*gam - 1.));

  G4double cross = alpha_rcl2*(4.2 - (2.*G4Log(fDelta)+1.)*rho*rho)/(gam+1.);
  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eplusTo3GammaOKVIModel::ComputeCrossSectionPerAtom(
                                    const G4ParticleDefinition*,
                                    G4double kineticEnergy, G4double Z,
				    G4double, G4double, G4double)
{
  // Calculates the cross section per atom of annihilation into two photons

  

  G4double cross = Z*ComputeCrossSectionPerElectron(kineticEnergy);
  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eplusTo3GammaOKVIModel::CrossSectionPerVolume(
					const G4Material* material,
					const G4ParticleDefinition*,
					      G4double kineticEnergy,
					      G4double, G4double)
{
  // Calculates the cross section per volume of annihilation into two photons
  
  G4double eDensity = material->GetElectronDensity();
  G4double cross = eDensity*ComputeCrossSectionPerElectron(kineticEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Polarisation of gamma according to M.H.L.Pryce and J.C.Ward, 
// Nature 4065 (1947) 435.

void 
G4eplusTo3GammaOKVIModel::SampleSecondaries(vector<G4DynamicParticle*>* vdp,
					    const G4MaterialCutsCouple*,
					    const G4DynamicParticle* dp,
					    G4double, G4double)
{

  G4double posiKinEnergy = dp->GetKineticEnergy();
  G4DynamicParticle *aGamma1, *aGamma2;
  G4DynamicParticle* aGamma3 = nullptr;
  G4double border;

  if(posiKinEnergy < 500*MeV) { 
    border = 1. - (electron_mass_c2)/(2*(posiKinEnergy + electron_mass_c2)); 
  } else {
    border = 1. - (100*electron_mass_c2)/(2*(posiKinEnergy + electron_mass_c2)); 
  }
  border = std::min(border, 0.9999);

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
   
  // Case at rest
  if(posiKinEnergy == 0.0) {
    G4double cost = 2.*rndmEngine->flat()-1.;
    G4double sint = sqrt((1. - cost)*(1. + cost));
    G4double phi  = twopi * rndmEngine->flat();
    G4ThreeVector dir(sint*cos(phi), sint*sin(phi), cost);
    phi = twopi * rndmEngine->flat();
    G4double cosphi = cos(phi);
    G4double sinphi = sin(phi);
    G4ThreeVector pol(cosphi, sinphi, 0.0);
    pol.rotateUz(dir);
    aGamma1 = new G4DynamicParticle(theGamma, dir, electron_mass_c2);
    aGamma1->SetPolarization(pol.x(),pol.y(),pol.z());
    aGamma2 = new G4DynamicParticle(theGamma,-dir, electron_mass_c2);
    pol.set(-sinphi, cosphi, 0.0);
    pol.rotateUz(dir);
    aGamma2->SetPolarization(pol.x(),pol.y(),pol.z());

  } else {

    G4ThreeVector posiDirection = dp->GetMomentumDirection();

    // (A.A.) LIMITS FOR 1st GAMMA
    G4double xmin = 0.01;
    G4double xmax = 0.667;   // CHANGE to 3/2 
  
    G4double d1, d0, x1, x2, dmax, x2min;

    // (A.A.) sampling of x1 x2 x3 (whole cycle of rejection)
    do {
      x1 = 1/((1/xmin) - ((1/xmin)-(1/xmax))*rndmEngine->flat()); 
      dmax = ComputeFS(posiKinEnergy, x1,1.-x1,border);
      x2min = 1.-x1;
      x2 = 1 - rndmEngine->flat()*(1-x2min);
      d1 = dmax*rndmEngine->flat();
      d0 = ComputeFS(posiKinEnergy,x1,x2,2-x1-x2);
    }  
    while(d0 < d1);

    G4double x3 = 2 - x1 - x2;

    //
    // angles between Gammas  
    //

    G4double psi13 = 2*asin(sqrt(std::abs((x1+x3-1)/(x1*x3))));
    G4double psi12 = 2*asin(sqrt(std::abs((x1+x2-1)/(x1*x2))));

    //          sin^t

    //G4double phi  = twopi * rndmEngine->flat();
    //G4double psi = acos(x3);                 // Angle of the plane

    //
    // kinematic of the created pair
    //

    G4double TotalAvailableEnergy = posiKinEnergy + 2.0*electron_mass_c2;

    G4double phot1Energy = 0.5*x1*TotalAvailableEnergy;      
    G4double phot2Energy = 0.5*x2*TotalAvailableEnergy;
    G4double phot3Energy = 0.5*x3*TotalAvailableEnergy;

    
    // DIRECTIONS
    
    // The azimuthal angles of ql and q3 with respect to some plane 
    // through the beam axis are generated at random. 

    G4ThreeVector phot1Direction(0, 0, 1);
    G4ThreeVector phot2Direction(0, sin(psi12), cos(psi12));
    G4ThreeVector phot3Direction(0, sin(psi13), cos(psi13));

    phot1Direction.rotateUz(posiDirection);                    
    phot2Direction.rotateUz(posiDirection);
    phot3Direction.rotateUz(posiDirection);

    aGamma1 = new G4DynamicParticle (theGamma,phot1Direction, phot1Energy);
    aGamma2 = new G4DynamicParticle (theGamma,phot2Direction, phot2Energy);
    aGamma3 = new G4DynamicParticle (theGamma,phot3Direction, phot3Energy);

    
                                                       //POLARIZATION - ???
   /*


    phi = twopi * rndmEngine->flat();
    G4double cosphi = cos(phi);
    G4double sinphi = sin(phi);
    G4ThreeVector pol(cosphi, sinphi, 0.0);
    pol.rotateUz(phot1Direction);
    aGamma1->SetPolarization(pol.x(),pol.y(),pol.z());

    G4double phot2Energy =(1.-epsil)*TotalAvailableEnergy;
    G4double posiP= sqrt(posiKinEnergy*(posiKinEnergy+2.*electron_mass_c2));
    G4ThreeVector dir = posiDirection*posiP - phot1Direction*phot1Energy;
    G4ThreeVector phot2Direction = dir.unit();

    // create G4DynamicParticle object for the particle2
    aGamma2 = new G4DynamicParticle (theGamma,phot2Direction, phot2Energy);

    //!!! likely problematic direction to be checked
    pol.set(-sinphi, cosphi, 0.0);
    pol.rotateUz(phot1Direction);
    cost = pol*phot2Direction;
    pol -= cost*phot2Direction;
    pol = pol.unit();
    aGamma2->SetPolarization(pol.x(),pol.y(),pol.z());
    
    */  

  }
  /*
    G4cout << "Annihilation in fly: e0= " << posiKinEnergy
    << " m= " << electron_mass_c2
    << " e1= " << phot1Energy 
    << " e2= " << phot2Energy << " dir= " <<  dir 
    << " -> " << phot1Direction << " " 
    << phot2Direction << G4endl;
  */
 
  vdp->push_back(aGamma1);
  vdp->push_back(aGamma2);
  if(aGamma3 != nullptr) { vdp->push_back(aGamma3); }

  // kill primary positron
  fParticleChange->SetProposedKineticEnergy(0.0);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
