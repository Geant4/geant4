// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ionIonisation.cc,v 1.4 1999-12-15 14:51:54 gunter Exp $
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
   : G4VContinuousDiscreteProcess(processName),
     ParticleMass(proton_mass_c2),Charge(eplus),
     dEdx(1.*MeV/mm),MinKineticEnergy(1.*keV)
{ PrintInfoDefinition() ; }

     
G4ionIonisation::~G4ionIonisation() 
{ }

G4double G4ionIonisation::GetConstraints(const G4DynamicParticle *aParticle,
                                              G4Material *aMaterial)
{
  // returns the Step limit
  // dRoverRange is the max. allowed relative range loss in one step
  // it calculates dEdx and the range as well....
  const G4double minstep=0.01*mm ;

  G4double KineticEnergy,StepLimit;

  Charge = aParticle->GetDefinition()->GetPDGCharge()/eplus ;

  KineticEnergy = aParticle->GetKineticEnergy();

  G4double massratio=proton_mass_c2/
           aParticle->GetDefinition()->GetPDGMass() ;

  G4double Tscaled= KineticEnergy*massratio ;
  G4double ChargeSquare = Charge*Charge ;

  dEdx=ComputedEdx(aParticle,aMaterial) ;
  StepLimit = 0.2*KineticEnergy/dEdx ;
  if(StepLimit < minstep)
    StepLimit = minstep ;

  return StepLimit ;
}

G4VParticleChange* G4ionIonisation::AlongStepDoIt(
                              const G4Track& trackData,const G4Step& stepData)
 // compute the energy loss after a step
{
  const G4DynamicParticle* aParticle;
  G4Material* aMaterial;
  G4double E,finalT,Step,ChargeSquare,MeanLoss ;

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;
 
  // get the actual (true) Step length from stepData
  Step = stepData.GetStepLength() ;

  aParticle = trackData.GetDynamicParticle() ;
  G4double massratio=proton_mass_c2/
           aParticle->GetDefinition()->GetPDGMass() ;
  ChargeSquare = Charge*Charge ;

  G4int index = aMaterial->GetIndex() ;
  E = aParticle->GetKineticEnergy() ;

  if(E < MinKineticEnergy) MeanLoss = E ;
  else
  {
     MeanLoss = Step*dEdx ;
     MeanLoss /= (massratio*ChargeSquare) ;
  }
  finalT = E - MeanLoss ;

  if(finalT < MinKineticEnergy) finalT = 0. ;

   //  kill the particle if the kinetic energy <= 0
  if (finalT <= 0. )
  {
    finalT = 0.;
    aParticleChange.SetStatusChange(fStopAndKill);
  }

  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(E-finalT) ;

  return &aParticleChange ;
}

 
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
   Charge=(aParticle->GetDefinition()->GetPDGCharge())/eplus;
   G4double ChargeSquare=Charge*Charge ;

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

G4double G4ionIonisation::ComputedEdx(const G4DynamicParticle* aParticle,
                     G4Material* material)
{
  // cuts for  electron ....................
  DeltaCutInKineticEnergy = G4Electron::Electron()->GetCutsInEnergy() ;

  G4double KineticEnergy , ionloss ;
  G4double  RateMass ;
  G4bool isOutRange ;
  const G4double twoln10 = 2.*log(10.) ;
  const G4double Factor = twopi_mc2_rcl2 ;
  const G4double bg2lim = 0.0169 , taulim = 8.4146e-3 ;

  RateMass = electron_mass_c2/proton_mass_c2 ;

    // get material parameters needed for the energy loss calculation

    G4double ElectronDensity,Eexc,Eexc2,Cden,Mden,Aden,X0den,X1den,taul ;
    G4double* ShellCorrectionVector;

    ElectronDensity = material->GetElectronDensity();
    Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
    Eexc2 = Eexc*Eexc ;
    Cden = material->GetIonisation()->GetCdensity();
    Mden = material->GetIonisation()->GetMdensity();
    Aden = material->GetIonisation()->GetAdensity();
    X0den = material->GetIonisation()->GetX0density();
    X1den = material->GetIonisation()->GetX1density();
    taul = material->GetIonisation()->GetTaul() ;
    ShellCorrectionVector = material->GetIonisation()->
                                          GetShellCorrectionVector();

    // get elements in the actual material,
    // they are needed for the low energy part ....

    const G4ElementVector* theElementVector=
                   material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector=
                   material->GetAtomicNumDensityVector() ;
    const G4int NumberOfElements=
                   material->GetNumberOfElements() ;

    // get  electron cut in kin. energy for the material
    DeltaCutInKineticEnergyNow =
         DeltaCutInKineticEnergy[material->GetIndex()] ;

    // some local variables -------------------
    G4double tau,tau0,Tmax,gamma,bg2,beta2,rcut,delta,x,sh ;

    KineticEnergy=aParticle->GetKineticEnergy();


      tau = KineticEnergy/proton_mass_c2 ;

      if ( tau < taul )
      //  low energy part , parametrized energy loss formulae
      {
        ionloss = 0. ;
        //  loop for the elements in the material
        for (G4int iel=0; iel<NumberOfElements; iel++)
        {
          const G4Element* element = (*theElementVector)(iel);

          if ( tau < element->GetIonisation()->GetTau0())
            ionloss += theAtomicNumDensityVector[iel]
                       *( element->GetIonisation()->GetAlow()*sqrt(tau)
                       +element->GetIonisation()->GetBlow()*tau) ;
          else
            ionloss += theAtomicNumDensityVector[iel]
                       *  element->GetIonisation()->GetClow()/sqrt(tau) ;
        }
      }
      else
      // high energy part , Bethe-Bloch formula
      {
        gamma = tau +1. ;
        bg2 = tau*(tau+2.) ;
        beta2 = bg2/(gamma*gamma) ;
        Tmax = 2.*electron_mass_c2*bg2
               /(1.+2.*gamma*RateMass+RateMass*RateMass) ;

        if ( DeltaCutInKineticEnergyNow < Tmax)
          rcut = DeltaCutInKineticEnergyNow/Tmax ;
        else
          rcut = 1.;

        ionloss = log(2.*electron_mass_c2*bg2*Tmax/Eexc2)

                  +log(rcut)-(1.+rcut)*beta2 ;

        // density correction

        x = log(bg2)/twoln10 ;
        if ( x < X0den )
          delta = 0. ;
        else
        {
          delta = twoln10*x - Cden ;
          if ( x < X1den )
            delta += Aden*pow((X1den-x),Mden) ;
        }

        // shell correction

        if ( bg2 > bg2lim ) {
          sh = 0. ;
          x = 1. ;
          for (G4int k=0; k<=2; k++) {
            x *= bg2 ;
            sh += ShellCorrectionVector[k]/x;
          }
        }
        else {
          sh = 0. ;
          x = 1. ;
          for (G4int k=0; k<=2; k++) {
             x *= bg2lim ;
             sh += ShellCorrectionVector[k]/x;
          }
          sh *= log(tau/taul)/log(taulim/taul) ;
        }

        // now you can compute the total ionization loss

        ionloss -= delta + sh ;
        ionloss *= Factor*ElectronDensity/beta2 ;
      }
      if ( ionloss <= 0.)
        ionloss = 0. ;

      dEdx = ionloss ;
   return dEdx ;
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

  G4cout << G4endl << GetProcessName() << ":  " << comments << G4endl;
}

