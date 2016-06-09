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
// $Id: G4LCapture.cc,v 1.14 2007-02-24 05:17:29 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Model: Low-energy Neutron Capture
// F.W. Jones, TRIUMF, 03-DEC-96
// 
// This is a prototype of a low-energy neutron capture process.
// Currently it is based on the GHEISHA routine CAPTUR,
// and conforms fairly closely to the original Fortran.
//
// HPW Capture using models now. the code comes from the
// original G4LCapture class.
//
// 25-JUN-98 FWJ: replaced missing Initialize for ParticleChange.
//

#include "globals.hh"
#include "G4LCapture.hh"
#include "Randomize.hh"

G4LCapture::G4LCapture() : 
   G4HadronicInteraction("G4LCapture")
{   
   SetMinEnergy( 0.0*GeV );
   SetMaxEnergy( DBL_MAX );
}

G4LCapture::~G4LCapture()
{
   theParticleChange.Clear();
}

G4HadFinalState*
G4LCapture::ApplyYourself(const G4HadProjectile & aTrack, G4Nucleus& targetNucleus)
{

   theParticleChange.Clear();
   theParticleChange.SetStatusChange(stopAndKill);

   G4double N = targetNucleus.GetN();
   G4double Z = targetNucleus.GetZ();

   const G4LorentzVector theMom = aTrack.Get4Momentum();
   G4double P = theMom.vect().mag()/GeV;
   G4double Px = theMom.vect().x()/GeV;
   G4double Py = theMom.vect().y()/GeV;
   G4double Pz = theMom.vect().z()/GeV;
   G4double E = theMom.e()/GeV;
   G4double E0 = aTrack.GetDefinition()->GetPDGMass()/GeV;
   G4double Q = aTrack.GetDefinition()->GetPDGCharge();
   if (verboseLevel > 1) {
      G4cout << "G4LCapture:ApplyYourself: incident particle:" << G4endl;
      G4cout << "P      " << P << " GeV/c" << G4endl;
      G4cout << "Px     " << Px << " GeV/c" << G4endl;
      G4cout << "Py     " << Py << " GeV/c" << G4endl;
      G4cout << "Pz     " << Pz << " GeV/c" << G4endl;
      G4cout << "E      " << E << " GeV" << G4endl;
      G4cout << "mass   " << E0 << " GeV" << G4endl;
      G4cout << "charge " << Q << G4endl;
   }

   // GHEISHA ADD operation to get total energy, mass, charge:

   if (verboseLevel > 1) {
      G4cout << "G4LCapture:ApplyYourself: material:" << G4endl;
      G4cout << "A      " << N << G4endl;
      G4cout << "Z   " << Z  << G4endl;
      G4cout << "atomic mass " << 
        Atomas(N, Z) << "GeV" << G4endl;
   }
   E = E + Atomas(N, Z);
   G4double E02 = E*E - P*P;
   E0 = std::sqrt(std::abs(E02));
   if (E02 < 0) E0 = -E0;
   Q = Q + Z;
   if (verboseLevel > 1) {
      G4cout << "G4LCapture:ApplyYourself: total:" << G4endl;
      G4cout << "E      " << E << " GeV" << G4endl;
      G4cout << "mass   " << E0 << " GeV" << G4endl;
      G4cout << "charge " << Q << G4endl;
   }
   Px = -Px;
   Py = -Py;
   Pz = -Pz;

   // Make a gamma...

   G4double p;
   if (Z == 1 && N == 1) { // special case for hydrogen
     p = 0.0022;
   } else {
     G4double ran = G4RandGauss::shoot();
     p = 0.0065 + ran*0.0010;
   }

   G4double ran1 = G4UniformRand();
   G4double ran2 = G4UniformRand();
   G4double cost = -1. + 2.*ran1;
   G4double sint = std::sqrt(std::abs(1. - cost*cost));
   G4double phi = ran2*twopi;

   G4double px = p*sint*std::sin(phi);
   G4double py = p*sint*std::cos(phi);
   G4double pz = p*cost;
   G4double e = p;
   // G4double e0 = 0.;  DHW 19 May 2011: variable set but not used

   G4double a = px*Px + py*Py + pz*Pz;
   a = (a/(E + E0) - e)/E0;

   px = px + a*Px;
   py = py + a*Py;
   pz = pz + a*Pz;

   G4DynamicParticle* aGamma;
   aGamma = new G4DynamicParticle(G4Gamma::GammaDefinition(),
                                  G4ThreeVector(px*GeV, py*GeV, pz*GeV));
   theParticleChange.AddSecondary(aGamma);

   // Make another gamma if there is sufficient energy left over...

   G4double xp = 0.008 - p;
   if (xp > 0.) {
     if (Z > 1 || N > 1) { 
       ran1 = G4UniformRand();
       ran2 = G4UniformRand();
       cost = -1. + 2.*ran1;
       sint = std::sqrt(std::abs(1. - cost*cost));
       phi = ran2*twopi;

       px = xp*sint*std::sin(phi);
       py = xp*sint*std::cos(phi);
       pz = xp*cost;
       e = xp;
       //  e0 = 0.;   DHW 19 May 2011: variable set but not used

       a = px*Px + py*Py + pz*Pz;
       a = (a/(E + E0) - e)/E0;

       px = px + a*Px;
       py = py + a*Py;
       pz = pz + a*Pz;

       aGamma = new G4DynamicParticle(G4Gamma::GammaDefinition(),
                                      G4ThreeVector(px*GeV, py*GeV, pz*GeV));
       theParticleChange.AddSecondary(aGamma);
     }
   }
   return &theParticleChange;
}
