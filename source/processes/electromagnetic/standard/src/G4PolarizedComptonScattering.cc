// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PolarizedComptonScattering.cc,v 1.2 1999-12-15 14:51:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ---------- G4PolarizedComptonScattering physics process --------
//                   by Vicente Lara, March 1998
// **************************************************************
//
// --------------------------------------------------------------

#include "G4PolarizedComptonScattering.hh"
#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// constructor
 
G4PolarizedComptonScattering::G4PolarizedComptonScattering(const G4String& processName)
  : G4ComptonScattering (processName)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4PolarizedComptonScattering::PostStepDoIt(const G4Track& aTrack,
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
   G4ComptonScattering::PostStepDoIt(aTrack,aStep);
       
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
            { epsilon   = exp(-alpha1*G4UniformRand());  // pow(epsilon0,G4UniformRand())
              epsilonsq = epsilon*epsilon; }
       else {
             epsilonsq = epsilon0sq + (1.- epsilon0sq)*G4UniformRand();
             epsilon   = sqrt(epsilonsq);
       };
       onecost = (1.- epsilon)/(epsilon*E0_m);
       sint2   = onecost*(2.-onecost);
       greject = 1. - epsilon*sint2/(1.+ epsilonsq);
   } while (greject < G4UniformRand());
 
   //
   // scattered gamma angles. ( Z - axis along the parent gamma)
   //

   G4double cosTeta = 1. - onecost , sinTeta = sqrt (sint2);
   G4double Phi     = twopi * G4UniformRand() ;
   G4double dirx = sinTeta*cos(Phi) , diry = sinTeta*sin(Phi) , dirz = cosTeta ;

   //
   // update G4VParticleChange for the scattered gamma 
   //
   
   G4double GammaEnergy1 = epsilon*GammaEnergy0;

   // New polarization
   G4ThreeVector GammaPolarization1 = SetNewPolarization(epsilon,sint2,Phi,
                                                         cosTeta,
							 GammaPolarization0);

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
 
G4ThreeVector G4PolarizedComptonScattering::SetNewPolarization(G4double EnergyRate, 
							       G4double sinsqrth,
							       G4double phi, 
							       G4double costheta, 
							       G4ThreeVector& GammaPolarization0) 
{
  G4double cosphi = cos(phi), sinphi = sin(phi);
  G4double ParallelIntensityPolar = EnergyRate + 1./EnergyRate + 2. - 4.*sinsqrth*cosphi*cosphi;
  G4double PerpendiIntensityPolar = EnergyRate + 1./EnergyRate - 2.;
  G4double PolarizationDegree = sqrt(sinsqrth*sinphi*sinphi + costheta*costheta);
  G4double sintheta = sqrt(sinsqrth);
  
  G4ThreeVector GammaPolarization1;
  // depolarization probability (1-P)
  if ( G4UniformRand() > 0.5*(PerpendiIntensityPolar/ParallelIntensityPolar) )
    {
     // Parallel to initial polarization
     GammaPolarization1.setX(PolarizationDegree);
     GammaPolarization1.setY(-sinsqrth*sinphi*cosphi/PolarizationDegree);
     GammaPolarization1.setZ(-sintheta*costheta*cosphi/PolarizationDegree);
    }
  else 
    {
     // Perpendicular to initial polarization
     GammaPolarization1.setX(0.);
     GammaPolarization1.setY(costheta/PolarizationDegree);
     GammaPolarization1.setZ(-sintheta*sinphi/PolarizationDegree);
    };
  
  return GammaPolarization1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedComptonScattering::SystemOfRefChange(G4ThreeVector& Direction0,
                                                     G4ThreeVector& Direction1,
						  G4ThreeVector& Polarization0,
						  G4ThreeVector& Polarization1)
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
   

  // Rotation along Z axe
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
