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
// File name:   G4eeToTwoGammaModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 18-04-05 Compute CrossSectionPerVolume (V.Ivanchenko)
// 06-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 29-06-06 Fix problem for zero energy incident positron (V.Ivanchenko) 
// 20-10-06 Add theGamma as a member (V.Ivanchenko)
// 18-01-20 Introduce thermal model of annihilation at rest (J.Allison)
//
//
// Class Description:
//
// Implementation of e+ annihilation into 2 gamma
//
// The secondaries Gamma energies are sampled using the Heitler cross section.
//
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1: The initial electron is assumed free and at rest if atomic PDF 
//         is not defined
//
// Note 2: The annihilation processes producing one or more than two photons are
//         ignored, as negligible compared to the two photons process.

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToTwoGammaModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackStatus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4EmParameters.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4bool G4eeToTwoGammaModel::fSampleAtomicPDF = false;

G4eeToTwoGammaModel::G4eeToTwoGammaModel(const G4ParticleDefinition*,
                                         const G4String& nam)
  : G4VEmModel(nam),
    pi_rcl2(pi*classic_electr_radius*classic_electr_radius)
{
  theGamma = G4Gamma::Gamma();
  fParticleChange = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToTwoGammaModel::~G4eeToTwoGammaModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToTwoGammaModel::Initialise(const G4ParticleDefinition*,
                                     const G4DataVector&)
{
  if(IsMaster()) {
    G4int verbose = G4EmParameters::Instance()->Verbose();
    // redo initialisation for each new run
    fSampleAtomicPDF = false;
    const auto& materialTable = G4Material::GetMaterialTable();
    for (const auto& material: *materialTable) {
      const G4double meanEnergyPerIonPair = material->GetIonisation()->GetMeanEnergyPerIonPair();
      if (meanEnergyPerIonPair > 0.) {
	fSampleAtomicPDF = true;
        if(verbose > 0) {
          G4cout << "### G4eeToTwoGammaModel: for " << material->GetName() << " mean energy per ion pair is " 
                 << meanEnergyPerIonPair/CLHEP::eV << " eV" << G4endl;
	}
      }
    }
  }
  // If no materials have meanEnergyPerIonPair set. This is probably the usual
  // case, since most applications are not senstive to the slight
  // non-collinearity of gammas in eeToTwoGamma. Do not issue any warning.

  if(fParticleChange) { return; }
  fParticleChange = GetParticleChangeForGamma();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4eeToTwoGammaModel::ComputeCrossSectionPerElectron(G4double kineticEnergy)
{
  // Calculates the cross section per electron of annihilation into two photons
  // from the Heilter formula.

  G4double ekin  = std::max(eV,kineticEnergy);   

  G4double tau   = ekin/electron_mass_c2;
  G4double gam   = tau + 1.0;
  G4double gamma2= gam*gam;
  G4double bg2   = tau * (tau+2.0);
  G4double bg    = sqrt(bg2);

  G4double cross = pi_rcl2*((gamma2+4*gam+1.)*G4Log(gam+bg) - (gam+3.)*bg)
                 / (bg2*(gam+1.));
  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eeToTwoGammaModel::ComputeCrossSectionPerAtom(
                                    const G4ParticleDefinition*,
                                    G4double kineticEnergy, G4double Z,
				    G4double, G4double, G4double)
{
  // Calculates the cross section per atom of annihilation into two photons
  return Z*ComputeCrossSectionPerElectron(kineticEnergy);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eeToTwoGammaModel::CrossSectionPerVolume(
					const G4Material* material,
					const G4ParticleDefinition*,
					      G4double kineticEnergy,
					      G4double, G4double)
{
  // Calculates the cross section per volume of annihilation into two photons
  return material->GetElectronDensity()*ComputeCrossSectionPerElectron(kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Polarisation of gamma according to M.H.L.Pryce and J.C.Ward,
// Nature 4065 (1947) 435.

void G4eeToTwoGammaModel::SampleSecondaries(vector<G4DynamicParticle*>* vdp,
					    const G4MaterialCutsCouple* pCutsCouple,
					    const G4DynamicParticle* dp,
					    G4double,
					    G4double)
{
  G4double posiKinEnergy = dp->GetKineticEnergy();
  G4DynamicParticle *aGamma1, *aGamma2;

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
   
  // Case at rest
  if(posiKinEnergy == 0.0) {

    const G4double eGamma = electron_mass_c2;

    // In rest frame of positronium gammas are back to back
    const G4ThreeVector& dir1 = G4RandomDirection();
    const G4ThreeVector& dir2 = -dir1;
    aGamma1 = new G4DynamicParticle(G4Gamma::Gamma(),dir1,eGamma);
    aGamma2 = new G4DynamicParticle(G4Gamma::Gamma(),dir2,eGamma);

    // In rest frame the gammas are polarised perpendicular to each other - see
    // Pryce and Ward, Nature No 4065 (1947) p.435.
    // Snyder et al, Physical Review 73 (1948) p.440.
    G4ThreeVector pol1 = (G4RandomDirection().cross(dir1)).unit();
    G4ThreeVector pol2 = (pol1.cross(dir2)).unit();

    // But the positronium is moving...
    // A positron in matter slows down and combines with an atomic electron to
    // make a neutral “atom” called positronium, about half the size of a normal
    // atom. I expect that when the energy of the positron is small enough,
    // less than the binding energy of positronium (6.8 eV), it is
    // energetically favourable for an electron from the outer orbitals of a
    // nearby atom or molecule to transfer and bind to the positron, as in an
    // ionic bond, leaving behind a mildly ionised nearby atom/molecule. I
    // would expect the positronium to come away with a kinetic energy of a
    // few eV on average. In its para (spin 0) state it annihilates into two
    // photons, which in the rest frame of the positronium are collinear
    // (back-to-back) due to momentum conservation. Because of the motion of the
    // positronium, photons will be not quite back-to-back in the laboratory.

    // The positroniuim acquires an energy of order its binding energy and
    // doesn't have time to thermalise. Nevertheless, here we approximate its
    // energy distribution by a Maxwell-Boltzman with mean energy <KE>. In terms
    // of a more familiar concept of temperature, and the law of equipartition
    // of energy of translational motion, <KE>=3kT/2. Each component of velocity
    // has a distribution exp(-mv^2/2kT), which is a Gaussian of mean zero
    // and variance kT/m=2<KE>/3m, where m is the positronium mass.

    // We take <KE> = material->GetIonisation()->GetMeanEnergyPerIonPair().

    if(fSampleAtomicPDF) {
      const G4Material* material = pCutsCouple->GetMaterial();
      const G4double meanEnergyPerIonPair = material->GetIonisation()->GetMeanEnergyPerIonPair();
      const G4double& meanKE = meanEnergyPerIonPair;  // Just an alias
      if (meanKE > 0.) {  // Positronium haas motion
	// Mass of positronium
	const G4double mass = 2.*electron_mass_c2;
	// Mean <KE>=3kT/2, as described above
	// const G4double T = 2.*meanKE/(3.*k_Boltzmann);
	// Component velocities: Gaussian, variance kT/m=2<KE>/3m.
	const G4double sigmav = std::sqrt(2.*meanKE/(3.*mass));
	// This is in units where c=1
	const G4double vx = G4RandGauss::shoot(0.,sigmav);
	const G4double vy = G4RandGauss::shoot(0.,sigmav);
	const G4double vz = G4RandGauss::shoot(0.,sigmav);
	const G4ThreeVector v(vx,vy,vz);  // In unit where c=1
	const G4ThreeVector& beta = v;    // so beta=v/c=v

	aGamma1->Set4Momentum(aGamma1->Get4Momentum().boost(beta));
	aGamma2->Set4Momentum(aGamma2->Get4Momentum().boost(beta));

	// Rotate polarisation vectors
	const G4ThreeVector& newDir1 = aGamma1->GetMomentumDirection();
	const G4ThreeVector& newDir2 = aGamma2->GetMomentumDirection();
	const G4ThreeVector& axis1 = dir1.cross(newDir1);  // No need to be unit
	const G4ThreeVector& axis2 = dir2.cross(newDir2);  // No need to be unit
	const G4double& angle1 = std::acos(dir1*newDir1);
	const G4double& angle2 = std::acos(dir2*newDir2);
	if (axis1 != G4ThreeVector()) pol1.rotate(axis1,angle1);
	if (axis2 != G4ThreeVector()) pol2.rotate(axis2,angle2);
      }
    }
    aGamma1->SetPolarization(pol1.x(),pol1.y(),pol1.z());
    aGamma2->SetPolarization(pol2.x(),pol2.y(),pol2.z());

  } else {  // Positron interacts in flight

    G4ThreeVector posiDirection = dp->GetMomentumDirection();

    G4double tau     = posiKinEnergy/electron_mass_c2;
    G4double gam     = tau + 1.0;
    G4double tau2    = tau + 2.0;
    G4double sqgrate = sqrt(tau/tau2)*0.5;
    G4double sqg2m1  = sqrt(tau*tau2);

    // limits of the energy sampling
    G4double epsilmin = 0.5 - sqgrate;
    G4double epsilmax = 0.5 + sqgrate;
    G4double epsilqot = epsilmax/epsilmin;

    //
    // sample the energy rate of the created gammas
    //
    G4double epsil, greject;

    do {
      epsil = epsilmin*G4Exp(G4Log(epsilqot)*rndmEngine->flat());
      greject = 1. - epsil + (2.*gam*epsil-1.)/(epsil*tau2*tau2);
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while( greject < rndmEngine->flat());

    //
    // scattered Gamma angles. ( Z - axis along the parent positron)
    //

    G4double cost = (epsil*tau2-1.)/(epsil*sqg2m1);
    if(std::abs(cost) > 1.0) {
      G4cout << "### G4eeToTwoGammaModel WARNING cost= " << cost
	     << " positron Ekin(MeV)= " << posiKinEnergy
	     << " gamma epsil= " << epsil
	     << G4endl;
      if(cost > 1.0) cost = 1.0;
      else cost = -1.0; 
    }
    G4double sint = sqrt((1.+cost)*(1.-cost));
    G4double phi  = twopi * rndmEngine->flat();

    //
    // kinematic of the created pair
    //

    G4double totalEnergy = posiKinEnergy + 2.0*electron_mass_c2;
    G4double phot1Energy = epsil*totalEnergy;

    G4ThreeVector phot1Direction(sint*cos(phi), sint*sin(phi), cost);
    phot1Direction.rotateUz(posiDirection);
    aGamma1 = new G4DynamicParticle (theGamma,phot1Direction, phot1Energy);
    phi = twopi * rndmEngine->flat();
    G4double cosphi = cos(phi);
    G4double sinphi = sin(phi);
    G4ThreeVector pol(cosphi, sinphi, 0.0);
    pol.rotateUz(phot1Direction);
    aGamma1->SetPolarization(pol.x(),pol.y(),pol.z());

    G4double phot2Energy =(1.-epsil)*totalEnergy;
    G4double posiP= sqrt(posiKinEnergy*(posiKinEnergy+2.*electron_mass_c2));
    G4ThreeVector dir = posiDirection*posiP - phot1Direction*phot1Energy;
    G4ThreeVector phot2Direction = dir.unit();

    // create G4DynamicParticle object for the particle2
    aGamma2 = new G4DynamicParticle (theGamma, phot2Direction, phot2Energy);

    //!!! likely problematic direction to be checked
    pol.set(-sinphi, cosphi, 0.0);
    pol.rotateUz(phot1Direction);
    cost = pol*phot2Direction;
    pol -= cost*phot2Direction;
    pol = pol.unit();
    aGamma2->SetPolarization(pol.x(),pol.y(),pol.z());
    /*
    G4cout << "Annihilation on fly: e0= " << posiKinEnergy
	   << " m= " << electron_mass_c2
	   << " e1= " << phot1Energy 
	   << " e2= " << phot2Energy << " dir= " <<  dir 
	   << " -> " << phot1Direction << " " 
	   << phot2Direction << G4endl;
    */
  }  
 
  vdp->push_back(aGamma1);
  vdp->push_back(aGamma2);

  // kill primary positron
  fParticleChange->SetProposedKineticEnergy(0.0);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
