// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LFission.cc,v 1.4 2000-08-03 08:49:30 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//

#include "globals.hh"
#include "G4LFission.hh"
#include "Randomize.hh"

G4LFission::G4LFission() : 
   G4HadronicInteraction()
{
   init();
   theParticleChange.SetNumberOfSecondaries(1000);
   
   SetMinEnergy( 0.0*GeV );
   SetMaxEnergy( DBL_MAX );
   
}

G4LFission::~G4LFission()
{
   theParticleChange.Clear();
}
 

void 
G4LFission::init()
{
   G4int i;
   G4double xx = 1. - 0.5;
   G4double xxx = sqrt(2.29*xx);
   spneut[0] = exp(-xx/0.965)*(exp(xxx) - exp(-xxx))/2.;
   for (i = 2; i <= 10; i++) {
      xx = i*1. - 0.5;
      xxx = sqrt(2.29*xx);
      spneut[i-1] = spneut[i-2] + exp(-xx/0.965)*(exp(xxx) - exp(-xxx))/2.;
   }
   for (i = 1; i <= 10; i++) {
      spneut[i-1] = spneut[i-1]/spneut[9];
      if (verboseLevel > 1) G4cout << "G4LFission::init: i=" << i << 
         " spneut=" << spneut[i-1] << G4endl;
   }
}


G4VParticleChange*
G4LFission::ApplyYourself(const G4Track & aTrack,G4Nucleus & targetNucleus)
{
   theParticleChange.Initialize(aTrack);

   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   const G4Material* aMaterial = aTrack.GetMaterial();

   G4double N = targetNucleus.GetN();
   G4double Z = targetNucleus.GetZ();
   //   theParticleChange.SetKillSignal(1);
   theParticleChange.SetStatusChange(fStopAndKill);

   G4double P = aParticle->GetTotalMomentum()/MeV;
   G4double Px = P*(aParticle->GetMomentumDirection().x());
   G4double Py = P*(aParticle->GetMomentumDirection().y());
   G4double Pz = P*(aParticle->GetMomentumDirection().z());
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
   E0 = sqrt(abs(E02));
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
   G4double avern = 2.569 + 0.559*log(e1);
   G4bool photofission = 0;      // For now
// Take the following value if photofission is not included
   if (!photofission) avern = 2.569 + 0.900*log(e1);

// Average number of gammas
   G4double averg = 9.500 + 0.600*log(e1);

   G4double ran = G4RandGauss::shoot();
// Number of neutrons
   G4int nn = avern + ran*1.23 + 0.5;
   ran = G4RandGauss::shoot();
// Number of gammas
   G4int ng = averg + ran*3. + 0.5;
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
      theParticleChange.AddSecondary(aNeutron);
   }

// Make secondary gammas and distribute kinetic energy
   G4DynamicParticle* aGamma;
   for (i = 1; i <= ng; i++) {
      ran = G4UniformRand();
      G4double ekin = -0.87*log(ran);
      exg = exg + ekin;
      aGamma = new G4DynamicParticle(G4Gamma::GammaDefinition(),
                                     G4ParticleMomentum(1.,0.,0.),
                                     ekin*MeV);
      theParticleChange.AddSecondary(aGamma);
   }

   G4double ex = exn + exg;

// Distribute momentum vectors and do Lorentz transformation

   G4Track* theSecondary;

   for (i = 1; i <= nn + ng; i++) {
      G4double ran1 = G4UniformRand();
      G4double ran2 = G4UniformRand();
      G4double cost = -1. + 2.*ran1;
      G4double sint = sqrt(abs(1. - cost*cost));
      G4double phi = ran2*twopi;
      //      G4cout << ran1 << " " << ran2 << G4endl;
      //      G4cout << cost << " " << sint << " " << phi << G4endl;
      theSecondary = theParticleChange.GetSecondary(i - 1);
      G4double pp = theSecondary->GetDynamicParticle()->GetTotalMomentum()/MeV;
      G4double px = pp*sint*sin(phi);
      G4double py = pp*sint*cos(phi);
      G4double pz = pp*cost;
      //      G4cout << pp << G4endl;
      //      G4cout << px << " " << py << " " << pz << G4endl;
      G4double e = theSecondary->GetTotalEnergy()/MeV;
      G4double e0 = theSecondary->GetDefinition()->GetPDGMass()/MeV;

      G4double a = px*Px + py*Py + pz*Pz;
      a = (a/(E + E0) - e)/E0;

      px = px + a*Px;
      py = py + a*Py;
      pz = pz + a*Pz;
      G4double p2 = px*px + py*py + pz*pz;
      pp = sqrt(p2);
      e = sqrt(e0*e0 + p2);
      G4double ekin = e - theSecondary->GetDefinition()->GetPDGMass()/MeV;
      theSecondary->SetMomentumDirection(G4ParticleMomentum(px/pp,
                                                            py/pp,
                                                            pz/pp));
      theSecondary->SetKineticEnergy(ekin*MeV);
   }
   
   return &theParticleChange;
}

// Computes atomic mass in MeV (translation of GHEISHA routine ATOMAS)
// Not optimized: conforms closely to original Fortran.

G4double
G4LFission::Atomas(const G4double A, const G4double Z)
{
   G4double rmel = G4Electron::ElectronDefinition()->GetPDGMass()/MeV;
   G4double rmp  = G4Proton::ProtonDefinition()->GetPDGMass()/MeV;
   G4double rmn  = G4Neutron::NeutronDefinition()->GetPDGMass()/MeV;
   G4double rmd  = G4Deuteron::DeuteronDefinition()->GetPDGMass()/MeV;
   G4double rma  = G4Alpha::AlphaDefinition()->GetPDGMass()/MeV;

   G4int ia = A + 0.5;
   if (ia < 1) return 0;
   G4int iz = Z + 0.5;
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

   G4double mass = (A - Z)*rmn + Z*rmp + Z*rmel
                   - 15.67*A
                   + 17.23*pow(A, 2./3.)
                   + 93.15*(A/2. - Z)*(A/2. - Z)/A
                   + 0.6984523*Z*Z/pow(A, 1./3.);
   G4int ipp = (ia - iz)%2;
   G4int izz = iz%2;
   if (ipp == izz) mass = mass + (ipp + izz -1)*12.*pow(A, -0.5);

   return mass;
}
