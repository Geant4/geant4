// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ionIonisation.cc,v 1.1 1999-01-07 16:11:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4ionIonisation physics process -----------
//                by Laszlo Urban, 08 Dec 1998 
// **************************************************************
// It is the first implementation of the ionisation for IONS    
// --------------------------------------------------------------
 

#include "G4ionIonisation.hh"
#include "G4UnitsTable.hh"

// constructor and destructor
 
G4ionIonisation::G4ionIonisation(const G4String& processName)
   : G4hEnergyLoss(processName)
{ PrintInfoDefinition() ; }

     
G4ionIonisation::~G4ionIonisation() 
{ }
 
 G4double G4ionIonisation::GetMeanFreePath(
                                           const G4Track& trackData,
                                           G4double previousStepSize,
                                           G4ForceCondition* condition)
{
   const G4DynamicParticle* aParticle ;
   G4Material* aMaterial ;
   G4double MeanFreePath;

   *condition = NotForced ;

   aParticle = trackData.GetDynamicParticle() ;
   aMaterial = trackData.GetMaterial() ;

   G4double KineticEnergy = aParticle->GetKineticEnergy() ;
   G4double ChargeSquare=(aParticle->GetDefinition()->GetPDGCharge())*
                         (aParticle->GetDefinition()->GetPDGCharge());

  // compute the (macroscopic) cross section first
 
    const G4ElementVector* theElementVector=
                           aMaterial->GetElementVector() ;
    const G4double* theAtomicNumDensityVector =
                           aMaterial->GetAtomicNumDensityVector();
    const G4int NumberOfElements=
                           aMaterial->GetNumberOfElements() ;
    G4int index = aMaterial->GetIndex() ;
 
    DeltaCutInKineticEnergy = G4Electron::Electron()->GetCutsInEnergy();
    DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[index] ;

    G4double sigma = 0. ;
    for (G4int iel=0; iel<NumberOfElements; iel++ )
    {
       sigma +=  theAtomicNumDensityVector[iel]*
                  ComputeMicroscopicCrossSection(aParticle,
                         KineticEnergy,
                        (*theElementVector)(iel)->GetZ() ) ;
    }

    sigma *= twopi_mc2_rcl2 * ChargeSquare ;  

    // mean free path = 1./macroscopic cross section
    MeanFreePath = sigma<=0 ? DBL_MAX : 1./sigma ;     

   return MeanFreePath ;
}

G4double G4ionIonisation::ComputeMicroscopicCrossSection(
                                 const G4DynamicParticle* aParticle,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber)
{
    G4double TotalEnergy,
             betasquare,
             MaxKineticEnergyTransfer,TotalCrossSection,tempvar;

    // get particle data ...................................
    ParticleMass = aParticle->GetDefinition()->GetPDGMass() ;
    TotalEnergy=KineticEnergy + ParticleMass;

    // some kinematics......................
    betasquare = KineticEnergy*(TotalEnergy+ParticleMass)
                 /(TotalEnergy*TotalEnergy);
    tempvar = ParticleMass+electron_mass_c2;
    MaxKineticEnergyTransfer = 2.*electron_mass_c2*KineticEnergy
                     *(TotalEnergy+ParticleMass)
                     /(tempvar*tempvar+2.*electron_mass_c2*KineticEnergy);

    // now you can calculate the total cross section ------------------
    if( MaxKineticEnergyTransfer > DeltaCutInKineticEnergyNow )
    {
       tempvar=DeltaCutInKineticEnergyNow/MaxKineticEnergyTransfer;
       TotalCrossSection = (1.-tempvar*(1.-betasquare*log(tempvar)))
                           /DeltaCutInKineticEnergyNow;

       TotalCrossSection *=  AtomicNumber/betasquare;
    }
    else
       TotalCrossSection= 0. ;
 
    return TotalCrossSection ;
}
 
 
G4VParticleChange* G4ionIonisation::PostStepDoIt(
                                              const G4Track& trackData,   
                                              const G4Step& stepData)         
{
  const G4DynamicParticle* aParticle ;
  G4Material* aMaterial;
  G4double KineticEnergy,TotalEnergy,TotalMomentum,
           betasquare,MaxKineticEnergyTransfer,
           DeltaKineticEnergy,DeltaTotalMomentum,costheta,sintheta,phi,
           dirx,diry,dirz,finalKineticEnergy,finalPx,finalPy,finalPz,
           x,xc,grej,Psquare,Esquare,summass,rate,grejc,finalMomentum ;

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;

  aParticle = trackData.GetDynamicParticle() ;

  KineticEnergy=aParticle->GetKineticEnergy();
  TotalEnergy=KineticEnergy + ParticleMass ;
  Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ;
  Esquare=TotalEnergy*TotalEnergy ;
  summass = ParticleMass + electron_mass_c2 ;    
  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection() ;

  DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[aMaterial->GetIndex()];

  // some kinematics......................

  betasquare=Psquare/Esquare ;
  MaxKineticEnergyTransfer = 2.*electron_mass_c2*Psquare
                      /(summass*summass+2.*electron_mass_c2*KineticEnergy);

  // sampling kinetic energy of the delta ray 

  if( MaxKineticEnergyTransfer <= DeltaCutInKineticEnergyNow )
  {
    // there is no change at all).....
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }
  else
  {
   // normal case ......................................
      xc=DeltaCutInKineticEnergyNow/MaxKineticEnergyTransfer ;
      rate=MaxKineticEnergyTransfer/TotalEnergy ;

   // sampling follows ...
     grejc=1.-betasquare*xc ;

     do {
          x=xc/(1.-(1.-xc)*G4UniformRand());
          grej=(1.-x*betasquare)/grejc ;
        } while( G4UniformRand()>grej );
   }
   
   DeltaKineticEnergy = x * MaxKineticEnergyTransfer ;

   if(DeltaKineticEnergy <= 0.)
     return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

   DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
                                               2. * electron_mass_c2 )) ;
   TotalMomentum = sqrt(Psquare) ;
   costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2)
            /(DeltaTotalMomentum * TotalMomentum) ;

   //  protection against costheta > 1 or < -1   ---------------
   if ( costheta < -1. ) 
          costheta = -1. ;
   if ( costheta > +1. ) 
          costheta = +1. ;

   //  direction of the delta electron  ........
   phi = twopi * G4UniformRand() ; 
   sintheta = sqrt((1.+costheta)*(1.-costheta));
   dirx = sintheta * cos(phi) ;
   diry = sintheta * sin(phi) ;
   dirz = costheta ;

   G4ThreeVector DeltaDirection(dirx,diry,dirz) ;
   DeltaDirection.rotateUz(ParticleDirection) ;

   // create G4DynamicParticle object for delta ray
   G4DynamicParticle *theDeltaRay = new G4DynamicParticle;
   theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
   theDeltaRay->SetMomentumDirection(
                   DeltaDirection.x(),DeltaDirection.y(),DeltaDirection.z()); 
   theDeltaRay->SetDefinition(G4Electron::Electron());

   // fill aParticleChange 
   finalKineticEnergy = KineticEnergy - DeltaKineticEnergy ;
   if (finalKineticEnergy > 0.)
     {
       // changed energy and momentum of the actual particle
       finalMomentum=sqrt(finalKineticEnergy*
                         (finalKineticEnergy+2.*ParticleMass)) ;

       finalPx = (TotalMomentum*ParticleDirection.x()
                 -DeltaTotalMomentum*DeltaDirection.x())/finalMomentum ; 
       finalPy = (TotalMomentum*ParticleDirection.y()
                 -DeltaTotalMomentum*DeltaDirection.y())/finalMomentum ; 
       finalPz = (TotalMomentum*ParticleDirection.z()
                 -DeltaTotalMomentum*DeltaDirection.z())/finalMomentum ; 

       aParticleChange.SetMomentumChange( finalPx,finalPy,finalPz );
     }
   else
     {
       finalKineticEnergy = 0. ;
             aParticleChange.SetStatusChange(fStopAndKill);
     }

   aParticleChange.SetEnergyChange( finalKineticEnergy );
   aParticleChange.SetNumberOfSecondaries(1);   
   aParticleChange.AddSecondary( theDeltaRay );
      
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

void G4ionIonisation::PrintInfoDefinition()
{
  G4String comments = "  Knock-on electron cross sections . ";
           comments += "\n         MeanFreePath is computed at tracking time.\n";
           comments += "         delta ray energy sampled from  differential Xsection.";

  G4cout << endl << GetProcessName() << ":  " << comments << endl;
}

