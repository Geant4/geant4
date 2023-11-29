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
//
// G4 Model: Low Energy Fission
// F.W. Jones, TRIUMF, 03-DEC-96
// 
// This is a prototype of a low-energy fission process.
// Currently it is based on the GHEISHA routine FISSIO,
// and conforms fairly closely to the original Fortran.
// Note: energy is in MeV and momentum is in MeV/c.
//
// use -scheme for elastic scattering: HPW, 20th June 1997
// the code comes mostly from the old Low-energy Fission class
//
// 25-JUN-98 FWJ: replaced missing Initialize for ParticleChange.

#include <iostream>

#include "G4LFission.hh"
#include "globals.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicsModelCatalog.hh"

G4LFission::G4LFission(const G4String& name)
  : G4HadronicInteraction(name), secID(-1)
{
  init();
  SetMinEnergy(0.0*GeV);
  SetMaxEnergy(DBL_MAX);
  G4PhysicsModelCatalog::GetModelID( "model_" + GetModelName() );
}


G4LFission::~G4LFission()
{
  theParticleChange.Clear();
}


void G4LFission::ModelDescription(std::ostream& outFile) const
{
  outFile << "G4LFission is one of the Low Energy Parameterized\n"
          << "(LEP) models used to implement neutron-induced fission of\n"
          << "nuclei.  It is a re-engineered version of the GHEISHA code\n"
          << "of H. Fesefeldt which emits neutrons and gammas but no\n"
          << "nuclear fragments.  The model is applicable to all incident\n"
          << "neutron energies.\n";
}

void G4LFission::init()
{
   G4int i;
   G4double xx = 1. - 0.5;
   G4double xxx = std::sqrt(2.29*xx);
   spneut[0] = G4Exp(-xx/0.965)*(G4Exp(xxx) - G4Exp(-xxx))/2.;
   for (i = 2; i <= 10; i++) {
      xx = i*1. - 0.5;
      xxx = std::sqrt(2.29*xx);
      spneut[i-1] = spneut[i-2] + G4Exp(-xx/0.965)*(G4Exp(xxx) - G4Exp(-xxx))/2.;
   }
   for (i = 1; i <= 10; i++) {
      spneut[i-1] = spneut[i-1]/spneut[9];
      if (verboseLevel > 1) G4cout << "G4LFission::init: i=" << i << 
         " spneut=" << spneut[i-1] << G4endl;
   }
}


G4HadFinalState* G4LFission::ApplyYourself(const G4HadProjectile& aTrack,
                                           G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  const G4HadProjectile* aParticle = &aTrack;

  G4double N = targetNucleus.GetA_asInt();
  G4double Z = targetNucleus.GetZ_asInt();
  theParticleChange.SetStatusChange(stopAndKill);

  G4double P = aParticle->GetTotalMomentum()/MeV;
  G4double Px = aParticle->Get4Momentum().vect().x();
  G4double Py = aParticle->Get4Momentum().vect().y();
  G4double Pz = aParticle->Get4Momentum().vect().z();
  G4double E = aParticle->GetTotalEnergy()/MeV;
  G4double E0 = aParticle->GetDefinition()->GetPDGMass()/MeV;
  G4double Q = aParticle->GetDefinition()->GetPDGCharge();
  if (verboseLevel > 1) {
      G4cout << "G4LFission:ApplyYourself: incident particle:" << G4endl;
      G4cout << "P      " << P << " MeV/c" << G4endl;
      G4cout << "Px     " << Px << " MeV/c" << G4endl;
      G4cout << "Py     " << Py << " MeV/c" << G4endl;
      G4cout << "Pz     " << Pz << " MeV/c" << G4endl;
      G4cout << "E      " << E << " MeV" << G4endl;
      G4cout << "mass   " << E0 << " MeV" << G4endl;
      G4cout << "charge " << Q << G4endl;
  }
  // GHEISHA ADD operation to get total energy, mass, charge:
   if (verboseLevel > 1) {
      G4cout << "G4LFission:ApplyYourself: material:" << G4endl;
      G4cout << "A      " << N << G4endl;
      G4cout << "Z      " << Z << G4endl;
      G4cout << "atomic mass " << 
        Atomas(N, Z) << "MeV" << G4endl;
   }
  E = E + Atomas(N, Z);
  G4double E02 = E*E - P*P;
  E0 = std::sqrt(std::abs(E02));
  if (E02 < 0) E0 = -E0;
  Q = Q + Z;
  if (verboseLevel > 1) {
      G4cout << "G4LFission:ApplyYourself: total:" << G4endl;
      G4cout << "E      " << E << " MeV" << G4endl;
      G4cout << "mass   " << E0 << " MeV" << G4endl;
      G4cout << "charge " << Q << G4endl;
  }
  Px = -Px;
  Py = -Py;
  Pz = -Pz;

  G4double e1 = aParticle->GetKineticEnergy()/MeV;
   if (e1 < 1.) e1 = 1.;

// Average number of neutrons
   G4double avern = 2.569 + 0.559*G4Log(e1);
   G4bool photofission = 0;      // For now
// Take the following value if photofission is not included
   if (!photofission) avern = 2.569 + 0.900*G4Log(e1);

// Average number of gammas
   G4double averg = 9.500 + 0.600*G4Log(e1);

   G4double ran = G4RandGauss::shoot();
// Number of neutrons
   G4int nn = static_cast<G4int>(avern + ran*1.23 + 0.5);
   ran = G4RandGauss::shoot();
// Number of gammas
   G4int ng = static_cast<G4int>(averg + ran*3. + 0.5);
   if (nn < 1) nn = 1;
   if (ng < 1) ng = 1;
   G4double exn = 0.;
   G4double exg = 0.;

// Make secondary neutrons and distribute kinetic energy
   G4DynamicParticle* aNeutron;
   G4int i;
   for (i = 1; i <= nn; i++) {
      ran = G4UniformRand();
      G4int j;
      for (j = 1; j <= 10; j++) {
         if (ran < spneut[j-1]) goto label12;
      }
      j = 10;
    label12:
      ran = G4UniformRand();
      G4double ekin = (j - 1)*1. + ran;
      exn = exn + ekin;
      aNeutron = new G4DynamicParticle(G4Neutron::NeutronDefinition(),
                                       G4ParticleMomentum(1.,0.,0.),
                                       ekin*MeV);
      theParticleChange.AddSecondary(aNeutron, secID);
   }

// Make secondary gammas and distribute kinetic energy
   G4DynamicParticle* aGamma;
   for (i = 1; i <= ng; i++) {
      ran = G4UniformRand();
      G4double ekin = -0.87*G4Log(ran);
      exg = exg + ekin;
      aGamma = new G4DynamicParticle(G4Gamma::GammaDefinition(),
                                     G4ParticleMomentum(1.,0.,0.),
                                     ekin*MeV);
      theParticleChange.AddSecondary(aGamma, secID);
   }

// Distribute momentum vectors and do Lorentz transformation

   G4HadSecondary* theSecondary;

   for (i = 1; i <= nn + ng; i++) {
      G4double ran1 = G4UniformRand();
      G4double ran2 = G4UniformRand();
      G4double cost = -1. + 2.*ran1;
      G4double sint = std::sqrt(std::abs(1. - cost*cost));
      G4double phi = ran2*twopi;
      //      G4cout << ran1 << " " << ran2 << G4endl;
      //      G4cout << cost << " " << sint << " " << phi << G4endl;
      theSecondary = theParticleChange.GetSecondary(i - 1);
      G4double pp = theSecondary->GetParticle()->GetTotalMomentum()/MeV;
      G4double px = pp*sint*std::sin(phi);
      G4double py = pp*sint*std::cos(phi);
      G4double pz = pp*cost;
      //      G4cout << pp << G4endl;
      //      G4cout << px << " " << py << " " << pz << G4endl;
      G4double e = theSecondary->GetParticle()->GetTotalEnergy()/MeV;
      G4double e0 = theSecondary->GetParticle()->GetDefinition()->GetPDGMass()/MeV;

      G4double a = px*Px + py*Py + pz*Pz;
      a = (a/(E + E0) - e)/E0;

      px = px + a*Px;
      py = py + a*Py;
      pz = pz + a*Pz;
      G4double p2 = px*px + py*py + pz*pz;
      pp = std::sqrt(p2);
      e = std::sqrt(e0*e0 + p2);
      G4double ekin = e - theSecondary->GetParticle()->GetDefinition()->GetPDGMass()/MeV;
      theSecondary->GetParticle()->SetMomentumDirection(G4ParticleMomentum(px/pp,
                                                            py/pp,
                                                            pz/pp));
      theSecondary->GetParticle()->SetKineticEnergy(ekin*MeV);
   }
   
  return &theParticleChange;
}

// Computes atomic mass in MeV (translation of GHEISHA routine ATOMAS)
// Not optimized: conforms closely to original Fortran.

G4double G4LFission::Atomas(const G4double A, const G4double Z)
{
  G4double rmel = G4Electron::ElectronDefinition()->GetPDGMass()/MeV;
  G4double rmp  = G4Proton::ProtonDefinition()->GetPDGMass()/MeV;
  G4double rmn  = G4Neutron::NeutronDefinition()->GetPDGMass()/MeV;
  G4double rmd  = G4Deuteron::DeuteronDefinition()->GetPDGMass()/MeV;
  G4double rma  = G4Alpha::AlphaDefinition()->GetPDGMass()/MeV;

  G4int ia = static_cast<G4int>(A + 0.5);
   if (ia < 1) return 0;
   G4int iz = static_cast<G4int>(Z + 0.5);
   if (iz < 0) return 0;
   if (iz > ia) return 0;

   if (ia == 1) {
      if (iz == 0) return rmn;          //neutron
      if (iz == 1) return rmp + rmel;   //Hydrogen
   }
   else if (ia == 2 && iz == 1) {
      return rmd;                       //Deuteron
   }
   else if (ia == 4 && iz == 2) {
      return rma;                       //Alpha
   }

  G4Pow* Pow=G4Pow::GetInstance();
  G4double mass = (A - Z)*rmn + Z*rmp + Z*rmel - 15.67*A
                  + 17.23*Pow->A23(A)
                  + 93.15*(A/2. - Z)*(A/2. - Z)/A
                  + 0.6984523*Z*Z/Pow->A13(A);
  G4int ipp = (ia - iz)%2;
  G4int izz = iz%2;
  if (ipp == izz) mass = mass + (ipp + izz -1)*12.*Pow->powA(A, -0.5);

  return mass;
}

const std::pair<G4double, G4double> G4LFission::GetFatalEnergyCheckLevels() const
{
  // max energy non-conservation is mass of heavy nucleus
  return std::pair<G4double, G4double>(10.0*perCent, 350.0*CLHEP::GeV);
}
