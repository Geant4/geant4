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
// G4RKFieldIntegrator
#include "G4RKFieldIntegrator.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NucleiProperties.hh"
#include "G4FermiMomentum.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4Nucleon.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

// Class G4RKFieldIntegrator
//*************************************************************************************************************************************

// only theActive are propagated, nothing else
// only theSpectators define the field, nothing else

void G4RKFieldIntegrator::Transport(G4KineticTrackVector &theActive, const G4KineticTrackVector &theSpectators, G4double theTimeStep)
{
   (void)theActive;
   (void)theSpectators;
   (void)theTimeStep;
}


G4double G4RKFieldIntegrator::CalculateTotalEnergy(const G4KineticTrackVector& Barions)
{
   const G4double Alpha  =  0.25/fermi/fermi;
   const G4double t1     = -7264.04*fermi*fermi*fermi;
   const G4double tGamma =  87.65*fermi*fermi*fermi*fermi*fermi*fermi;
//   const G4double Gamma  =  1.676;
   const G4double Vo     = -0.498*fermi;
   const G4double GammaY =  1.4*fermi;

   G4double Etot = 0;
   G4int nBarion = (G4int)Barions.size();
   for(G4int c1 = 0; c1 < nBarion; ++c1)
      {
      G4KineticTrack* p1 = Barions.operator[](c1);
   // Ekin
      Etot += p1->Get4Momentum().e();
      for(G4int c2 = c1 + 1; c2 < nBarion; ++c2)
         {
         G4KineticTrack* p2 = Barions.operator[](c2);
         G4double r12 = (p1->GetPosition() - p2->GetPosition()).mag()*fermi;

         //  Esk2
         Etot += t1*G4Pow::GetInstance()->A23(Alpha/pi)*G4Exp(-Alpha*r12*r12);

         // Eyuk
         Etot += Vo*0.5/r12*G4Exp(1/(4*Alpha*GammaY*GammaY))*
            (G4Exp(-r12/GammaY)*(1 - Erf(0.5/GammaY/std::sqrt(Alpha) - std::sqrt(Alpha)*r12)) -
             G4Exp( r12/GammaY)*(1 - Erf(0.5/GammaY/std::sqrt(Alpha) + std::sqrt(Alpha)*r12)));

         // Ecoul
         Etot += 1.44*p1->GetDefinition()->GetPDGCharge()*p2->GetDefinition()->GetPDGCharge()/r12*Erf(std::sqrt(Alpha)*r12);

         // Epaul
         Etot = 0;

         for(G4int c3 = c2 + 1; c3 < nBarion; c3++)
            {
            G4KineticTrack* p3 = Barions.operator[](c3);
            G4double r13 = (p1->GetPosition() - p3->GetPosition()).mag()*fermi;

            // Esk3
            Etot  = tGamma*G4Pow::GetInstance()->powA(4*Alpha*Alpha/3/pi/pi, 1.5)*G4Exp(-Alpha*(r12*r12 + r13*r13));
            }
         }
      }
   return Etot;
}

//************************************************************************************************
// originated from the Numerical recipes error function
G4double G4RKFieldIntegrator::Erf(G4double X)
{
   const G4double Z1 = 1;
   const G4double HF = Z1/2;
   const G4double C1 = 0.56418958;

   const G4double P10 = +3.6767877;
   const G4double Q10 = +3.2584593;
   const G4double P11 = -9.7970465E-2;

//   static G4ThreadLocal G4double P2[5] = { 7.3738883, 6.8650185,  3.0317993, 0.56316962, 4.3187787e-5 };
//   static G4ThreadLocal G4double Q2[5] = { 7.3739609, 15.184908, 12.79553,   5.3542168,  1. };
   const G4double P2[5] = { 7.3738883, 6.8650185,  3.0317993, 0.56316962, 4.3187787e-5 };
   const G4double Q2[5] = { 7.3739609, 15.184908, 12.79553,   5.3542168,  1. };

   const G4double P30 = -1.2436854E-1;
   const G4double Q30 = +4.4091706E-1;
   const G4double P31 = -9.6821036E-2;

   G4double V = std::abs(X);
   G4double H;
   G4double Y;
   G4int c1;

   if(V < HF)
      {
      Y = V*V;
      H = X*(P10 + P11*Y)/(Q10+Y);
      }
   else
      {
      if(V < 4)
         {
	 G4double AP = P2[4];
	 G4double AQ = Q2[4];
	 for(c1 = 3; c1 >= 0; c1--)
            {
            AP = P2[c1] + V*AP;
            AQ = Q2[c1] + V*AQ;
            }
	 H = 1 - G4Exp(-V*V)*AP/AQ;
	 }
      else
        {
        Y = 1./V*V;
        H = 1 - G4Exp(-V*V)*(C1+Y*(P30 + P31*Y)/(Q30 + Y))/V;
        }
     if (X < 0)
        H = -H;
     }
   return H;
}

//************************************************************************************************
//This is a QMD version to calculate excitation energy of a fragment,
//which consists from G4KTV &the Particles
/*
G4double G4RKFieldIntegrator::GetExcitationEnergy(const G4KineticTrackVector &theParticles)
{
   // Excitation energy of a fragment consisting from A nucleons and Z protons
   // is Etot - Z*Mp - (A - Z)*Mn - B(A, Z), where B(A,Z) is the binding energy of fragment
   //  and Mp, Mn are proton and neutron mass, respectively.
   G4int NZ = 0;
   G4int NA = 0;
   G4double Etot = CalculateTotalEnergy(theParticles);
   for(G4int cParticle = 0; cParticle < theParticles.length(); cParticle++)
      {
      G4KineticTrack* pKineticTrack = theParticles.at(cParticle);
      G4int Encoding =  std::abs(pKineticTrack->GetDefinition()->GetPDGEncoding());
      if (Encoding == 2212)
          NZ++, NA++;
      if (Encoding == 2112)
          NA++;
      Etot -= pKineticTrack->GetDefinition()->GetPDGMass();
      }
   return Etot - G4NucleiProperties::GetBindingEnergy(NZ, NA);
}
*/

//*************************************************************************************************************************************
//This is a simplified method to get excitation energy of a residual
// nucleus with nHitNucleons.
G4double G4RKFieldIntegrator::GetExcitationEnergy(G4int nHitNucleons, const G4KineticTrackVector &)
{
   const G4double MeanE = 50;
   G4double Sum = 0;
   for(G4int c1 = 0; c1 < nHitNucleons; ++c1)
       {
       Sum += -MeanE*G4Log(G4UniformRand());
       }
   return Sum;
}
//*************************************************************************************************************************************

/*
//This is free propagation of particles for CASCADE mode. Target nucleons should be frozen
void G4RKFieldIntegrator::Integrate(G4KineticTrackVector& theParticles)
   {
   for(G4int cParticle = 0; cParticle < theParticles.length(); ++cParticle)
      {
      G4KineticTrack* pKineticTrack = theParticles.at(cParticle);
      pKineticTrack->SetPosition(pKineticTrack->GetPosition() + theTimeStep*pKineticTrack->Get4Momentum().boostVector());
      }
   }
*/
//*************************************************************************************************************************************

void G4RKFieldIntegrator::Integrate(const G4KineticTrackVector& theBarions, G4double theTimeStep)
{
   for(std::size_t cParticle = 0; cParticle < theBarions.size(); ++cParticle)
      {
      G4KineticTrack* pKineticTrack = theBarions[cParticle];
      pKineticTrack->SetPosition(pKineticTrack->GetPosition() + theTimeStep*pKineticTrack->Get4Momentum().boostVector());
      }
}

//*************************************************************************************************************************************

// constant to calculate theCoulomb barrier
const G4double G4RKFieldIntegrator::coulomb = 1.44 / 1.14 * MeV;

// kaon's potential constant (real part only)
// 0.35 + i0.82 or 0.63 + i0.89 fermi
const G4double G4RKFieldIntegrator::a_kaon = 0.35;

// pion's potential constant (real part only)
//!! for pions it has todiffer from kaons
// 0.35 + i0.82 or 0.63 + i0.89 fermi
const G4double G4RKFieldIntegrator::a_pion = 0.35;

// antiproton's potential constant (real part only)
// 1.53 + i2.50 fermi
const G4double G4RKFieldIntegrator::a_antiproton = 1.53;

// methods for calculating potentials for different types of particles
// aPosition is relative to the nucleus center
G4double G4RKFieldIntegrator::GetNeutronPotential(G4double )
{
   /*
   const G4double Mn  = 939.56563 * MeV; // mass of nuetron

   G4VNuclearDensity *theDencity;
   if(theA < 17) theDencity = new G4NuclearShellModelDensity(theA, theZ);
   else          theDencity = new G4NuclearFermiDensity(theA, theZ);

   // GetDencity() accepts only G4ThreeVector so build it:
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);
   delete theDencity;

   G4FermiMomentum *fm = new G4FermiMomentum();
   fm->Init(theA, theZ);
   G4double fermiMomentum = fm->GetFermiMomentum(density);
   delete fm;

   return sqr(fermiMomentum)/(2 * Mn)
      + G4CreateNucleus::GetBindingEnergy(theZ, theA)/theA;
      //+ G4NucleiProperties::GetBindingEnergy(theZ, theA)/theA;
   */

   return 0.0;
}

G4double G4RKFieldIntegrator::GetProtonPotential(G4double )
{
   /*
   // calculate Coulomb barrier value
   G4double theCoulombBarrier = coulomb * theZ/(1. + G4Pow::GetInstance()->Z13(theA));
   const G4double Mp  = 938.27231 * MeV; // mass of proton

   G4VNuclearDensity *theDencity;
   if(theA < 17) theDencity = new G4NuclearShellModelDensity(theA, theZ);
   else          theDencity = new G4NuclearFermiDensity(theA, theZ);

   // GetDencity() accepts only G4ThreeVector so build it:
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);
   delete theDencity;

   G4FermiMomentum *fm = new G4FermiMomentum();
   fm->Init(theA, theZ);
   G4double fermiMomentum = fm->GetFermiMomentum(density);
   delete fm;

   return sqr(fermiMomentum)/ (2 * Mp)
      + G4CreateNucleus::GetBindingEnergy(theZ, theA)/theA;
      //+ G4NucleiProperties::GetBindingEnergy(theZ, theA)/theA
      + theCoulombBarrier;
   */

   return 0.0;
}

G4double G4RKFieldIntegrator::GetAntiprotonPotential(G4double )
{
   /*
   //G4double theM = G4NucleiProperties::GetAtomicMass(theA, theZ);
   G4double theM = theZ * G4Proton::Proton()->GetPDGMass()
      + (theA - theZ) * G4Neutron::Neutron()->GetPDGMass()
      + G4CreateNucleus::GetBindingEnergy(theZ, theA);

   const G4double Mp  = 938.27231 * MeV; // mass of proton
   G4double mu = (theM * Mp)/(theM + Mp);

   // antiproton's potential coefficient
   //   V = coeff_antiproton * nucleus_density
   G4double coeff_antiproton = -2.*pi/mu * (1. + Mp) * a_antiproton;

   G4VNuclearDensity *theDencity;
   if(theA < 17) theDencity = new G4NuclearShellModelDensity(theA, theZ);
   else          theDencity = new G4NuclearFermiDensity(theA, theZ);

   // GetDencity() accepts only G4ThreeVector so build it:
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);
   delete theDencity;

   return coeff_antiproton * density;
   */

   return 0.0;
}

G4double G4RKFieldIntegrator::GetKaonPotential(G4double )
{
   /*
   //G4double theM = G4NucleiProperties::GetAtomicMass(theA, theZ);
   G4double theM = theZ * G4Proton::Proton()->GetPDGMass()
      + (theA - theZ) * G4Neutron::Neutron()->GetPDGMass()
      + G4CreateNucleus::GetBindingEnergy(theZ, theA);

   const G4double Mk  = 496. * MeV;      // mass of "kaon"
   G4double mu = (theM * Mk)/(theM + Mk);

   // kaon's potential coefficient
   //   V = coeff_kaon * nucleus_density
   G4double coeff_kaon = -2.*pi/mu * (1. + Mk/theM) * a_kaon;

   G4VNuclearDensity *theDencity;
   if(theA < 17) theDencity = new G4NuclearShellModelDensity(theA, theZ);
   else          theDencity = new G4NuclearFermiDensity(theA, theZ);

   // GetDencity() accepts only G4ThreeVector so build it:
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);
   delete theDencity;

   return coeff_kaon * density;
   */

   return 0.0;
}

G4double G4RKFieldIntegrator::GetPionPotential(G4double )
{
   /*
   //G4double theM = G4NucleiProperties::GetAtomicMass(theA, theZ);
   G4double theM = theZ * G4Proton::Proton()->GetPDGMass()
      + (theA - theZ) * G4Neutron::Neutron()->GetPDGMass()
      + G4CreateNucleus::GetBindingEnergy(theZ, theA);

   const G4double Mpi = 139. * MeV;      // mass of "pion"
   G4double mu = (theM * Mpi)/(theM + Mpi);

   // pion's potential coefficient
   //   V = coeff_pion * nucleus_density
   G4double coeff_pion = -2.*pi/mu * (1. + Mpi) * a_pion;

   G4VNuclearDensity *theDencity;
   if(theA < 17) theDencity = new G4NuclearShellModelDensity(theA, theZ);
   else          theDencity = new G4NuclearFermiDensity(theA, theZ);

   // GetDencity() accepts only G4ThreeVector so build it:
   G4ThreeVector aPosition(0.0, 0.0, radius);
   G4double density = theDencity->GetDensity(aPosition);
   delete theDencity;

   return coeff_pion * density;
   */

   return 0.0;
}
