// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyPolarizedCompton.cc,v 1.1 2001-05-23 16:39:27 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group

// --------- G4LowEnergyPolarizedCompton class -----
//
//           by G.Depaola & F.Longo (21 may 2001)
//
// ************************************************************
//

// Design of G4PolarizedComptonScattering by V.Lara (1998)
// Corrections by Rui Curado da Silva (2000)
// New Implementation by G.Depaola & F.Longo      
//
// - sampling of phi
// - polarization of scattered photon
//
// --------------------------------------------------------------

#include "G4LowEnergyPolarizedCompton.hh"
#include "G4Electron.hh"
#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// constructor
 
G4LowEnergyPolarizedCompton::G4LowEnergyPolarizedCompton(const G4String& processName)
  : G4LowEnergyCompton (processName)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4LowEnergyPolarizedCompton::PostStepDoIt(const G4Track& aTrack,
							     const G4Step&  aStep)
  //
  // The scattered gamma energy is sampled according to Klein - Nishina formula.
// The random number techniques of Butcher & Messel are used (Nuc Phys 20(1960),15).
// GEANT4 internal units
//
// Note : Effects due to binding of atomic electrons are negliged.
 
{
   aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();

   const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle(); 
   
   G4ThreeVector GammaPolarization0 = aDynamicGamma->GetPolarization();  
   
   if (abs(GammaPolarization0.mag() - 1.e0) > 1.e-14)
     G4LowEnergyCompton::PostStepDoIt(aTrack,aStep);
   
   G4double GammaEnergy0 = aDynamicGamma->GetKineticEnergy();
   G4double E0_m = GammaEnergy0 / electron_mass_c2 ;
   
   G4ParticleMomentum GammaDirection0 = aDynamicGamma->GetMomentumDirection();
 

   //
   // sample the energy rate of the scattered gamma 
   //
   
   G4double epsilon, epsilonsq, onecost, sint2, greject ;

   G4double epsilon0 = 1./(1. + 2*E0_m) , epsilon0sq = epsilon0*epsilon0;
   G4double alpha1   = - log(epsilon0)  , alpha2 = 0.5*(1.- epsilon0sq);
   
   do {
       if ( alpha1/(alpha1+alpha2) > G4UniformRand() )
	 { 
	   epsilon   = exp(-alpha1*G4UniformRand());  
	   epsilonsq = epsilon*epsilon; 
	 }
       else 
	 {
	   epsilonsq = epsilon0sq + (1.- epsilon0sq)*G4UniformRand();
	   epsilon   = sqrt(epsilonsq);
	 };
       onecost = (1.- epsilon)/(epsilon*E0_m);
       sint2   = onecost*(2.-onecost);
       greject = 1. - epsilon*sint2/(1.+ epsilonsq);
   } while (greject < G4UniformRand());
   
   
   
   // ****************************************************
   //		Phi determination
   // ****************************************************
   
   
   G4double Phi = SetPhi(epsilon,sint2);
   
   
   //
   // scattered gamma angles. ( Z - axis along the parent gamma)
   //
   
   G4double cosTeta = 1. - onecost, sinTeta = sqrt (sint2);
   
   G4double dirx = sinTeta*cos(Phi), diry = sinTeta*sin(Phi), dirz = cosTeta ;
   
   
   //
   // update G4VParticleChange for the scattered gamma 
   //
   
   G4double GammaEnergy1 = epsilon*GammaEnergy0;
   
   // New polarization
   
   G4ThreeVector GammaPolarization1 = SetNewPolarization(epsilon,sint2,Phi,
                                                         cosTeta);
   
   // Set new direction 
   G4ThreeVector GammaDirection1 ( dirx,diry,dirz );
   


   // Change reference frame.
   SystemOfRefChange(GammaDirection0,GammaDirection1,
                     GammaPolarization0,GammaPolarization1);
   

   if (GammaEnergy1 > 0.)
     {
       aParticleChange.SetEnergyChange( GammaEnergy1 ) ;
     }
   else
     {    
       aParticleChange.SetEnergyChange(0.) ;
       aParticleChange.SetStatusChange(fStopAndKill);
     }
       
   //
   // kinematic of the scattered electron
   //
   
   G4double ElecKineEnergy = GammaEnergy0 - GammaEnergy1 ;
      
   if((G4EnergyLossTables::GetRange(G4Electron::Electron(),
       ElecKineEnergy,aMaterial)>aStep.GetPostStepPoint()->GetSafety())
      ||
      (ElecKineEnergy > 
	(G4Electron::Electron()->GetCutsInEnergy())[aMaterial->GetIndex()]))              
     {
       G4double ElecMomentum = sqrt(ElecKineEnergy*(ElecKineEnergy+2.*electron_mass_c2));
       G4ThreeVector ElecDirection (
				    (GammaEnergy0*GammaDirection0 - GammaEnergy1*GammaDirection1)*(1./ElecMomentum) );
       
       // create G4DynamicParticle object for the electron.  
        G4DynamicParticle* aElectron= new G4DynamicParticle (G4Electron::Electron(),
                                                        ElecDirection, ElecKineEnergy) ;

        aParticleChange.SetNumberOfSecondaries(1) ;
        aParticleChange.AddSecondary( aElectron ) ;
        aParticleChange.SetLocalEnergyDeposit (0.) ; 
       }
   else
     {
       aParticleChange.SetNumberOfSecondaries(0) ;
         aParticleChange.SetLocalEnergyDeposit (ElecKineEnergy) ;
       }


   //  Reset NbOfInteractionLengthLeft and return aParticleChange

   return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyPolarizedCompton::SetPhi(G4double EnergyRate,
					     G4double sinsqrth)
{
  G4double Rand1;
  G4double Rand2;
  G4double PhiProbability;
  G4double Phi;
  G4double a, b;

  do
    {
      Rand1 = G4UniformRand();
      Rand2 = G4UniformRand();
      PhiProbability=0.;
      Phi = twopi*Rand1;
      
      a = 2*sinsqrth;
      b = EnergyRate + 1/EnergyRate;
      
      PhiProbability = 1 - (a/b)*(cos(Phi)*cos(Phi)); 
    }
  while ( Rand2 > PhiProbability );
  
  return Phi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LowEnergyPolarizedCompton::SetNewPolarization(G4double epsilon, G4double sinsqrth, G4double phi,G4double costheta) 
{
  G4double Rand1;
  G4double Rand2;
  G4double cosphi = cos(phi), sinphi = sin(phi);
  G4double sintheta = sqrt(sinsqrth);
  G4double cossqrphi = cosphi*cosphi;
  G4double cossqrth = 1.-sinsqrth;
  G4double sinsqrphi = sinphi*sinphi;
  G4double Normalisation = sqrt(1-cossqrphi*sinsqrth);
  
  // Determination of Theta 
  
  G4double ThetaProbability;
  G4double Theta;
  G4double a, b;

  do
    {
      Rand1 = G4UniformRand();
      Rand2 = G4UniformRand();
      ThetaProbability=0.;
      Theta = twopi*Rand1;
      
      a = 4;
      b = (epsilon + 1/epsilon) - 2;
      
      ThetaProbability = (b + a*cos(Theta)*cos(Theta))/(a+b); 
    }
  while ( Rand2 > ThetaProbability );
  
  G4double cosTheta = cos(Theta);

  G4double cosbeta = cosTheta/Normalisation;
  G4double sinbeta = sqrt(1-cosbeta*cosbeta);
  G4ThreeVector GammaPolarization1;

  G4double Xparallel = Normalisation*cosbeta;
  G4double Yparallel = -(sinsqrth*cosphi*sinphi)*cosbeta/Normalisation;
  G4double Zparallel = -(costheta*sintheta*cosphi)*cosbeta/Normalisation;

  G4double Xperpendicular = 0.;
  G4double Yperpendicular = (costheta)*sinbeta/Normalisation;
  G4double Zperpendicular = -(sintheta*sinphi)*sinbeta/Normalisation;

  G4double Xtotal = (Xparallel + Xperpendicular);
  G4double Ytotal = (Yparallel + Yperpendicular);
  G4double Ztotal = (Zparallel + Zperpendicular);

  GammaPolarization1.setX(Xtotal);
  GammaPolarization1.setY(Ytotal);
  GammaPolarization1.setZ(Ztotal);

  return GammaPolarization1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



void G4LowEnergyPolarizedCompton::SystemOfRefChange
(G4ThreeVector& Direction0,G4ThreeVector& Direction1,
 G4ThreeVector& Polarization0,G4ThreeVector& Polarization1)

{
  // Angles for go back to the original RS
  G4double cosTeta0 = Direction0.cosTheta(), sinTeta0 = sin(Direction0.theta());
  G4double cosPhi0  = cos(Direction0.phi()), sinPhi0  = sin(Direction0.phi());

  G4double cosPsi, sinPsi;

  if (sinTeta0 != 0. ) {
    cosPsi = -Polarization0.z()/sinTeta0;
    if (cosPhi0 != 0. ) {
      sinPsi = (Polarization0.y() - cosTeta0*sinPhi0*cosPsi)/cosPhi0;
    } else {
      sinPsi = -Polarization0.x()/sinPhi0;
    }
  } else {
    cosPsi = Polarization0.x()/cosTeta0;
    sinPsi = Polarization0.y();
  }
  G4double Psi = atan(sinPsi/cosPsi); 
  

  // Rotation along Z axis
  Direction1.rotateZ(Psi);
  // 
  Direction1.rotateUz(Direction0);
  aParticleChange.SetMomentumChange( Direction1 ) ;  

  // 3 Euler angles rotation for scattered photon polarization

  Polarization1.rotateZ(Psi);
  Polarization1.rotateUz(Direction0);
  aParticleChange.SetPolarizationChange( Polarization1 );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....














