// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hIonisation.cc,v 1.3 1999-04-13 09:05:42 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hIonisation physics process -----------
//                by Laszlo Urban, 30 May 1997 
// **************************************************************
// It is the first implementation of the NEW IONISATION PROCESS.
// It calculates the ionisation of charged hadrons.
// **************************************************************
// corrected by L.Urban on 24/09/97
// several bugs corrected by L.Urban on 13/01/98
// 07-04-98: remove 'tracking cut' of the ionizing particle, MMa
// 22/10/98: cleanup L.Urban
// 02/02/99: bugs fixed , L.Urban
// --------------------------------------------------------------
 

#include "G4hIonisation.hh"
#include "G4UnitsTable.hh"

// constructor and destructor
 
G4hIonisation::G4hIonisation(const G4String& processName)
   : G4hEnergyLoss(processName),
     theMeanFreePathTable(NULL),
     LowestKineticEnergy(1.00*keV),
     HighestKineticEnergy(100.*TeV),
     TotBin(100),
     theProton (G4Proton::Proton()),
     theAntiProton (G4AntiProton::AntiProton()),
     theElectron ( G4Electron::Electron() )
{ }
     
G4hIonisation::~G4hIonisation() 
{
     if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
     }
}
 
void G4hIonisation::SetPhysicsTableBining(G4double lowE, G4double highE,
                                 G4int nBins)
{
  LowestKineticEnergy = lowE;  HighestKineticEnergy = highE;
  TotBin = nBins ;
}

// methods.............................................

void G4hIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
//  just call BuildLossTable+BuildLambdaTable
{
  ParticleMass = aParticleType.GetPDGMass() ;

  Charge = aParticleType.GetPDGCharge();

  G4double ElectronCutInRange = G4Electron::Electron()->GetCuts(); 

  DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;

  if(Charge>0.)
  {
    if( (ptableElectronCutInRange != ElectronCutInRange)  
                       || (theDEDXpTable == NULL))
    {
      BuildLossTable(aParticleType) ;
      RecorderOfpProcess[CounterOfpProcess] = theLossTable ;
      CounterOfpProcess++;
    }
  }
  else
  {
    if( (pbartableElectronCutInRange != ElectronCutInRange)  
                        || (theDEDXpbarTable == NULL))
    {
      BuildLossTable(aParticleType) ;
      RecorderOfpbarProcess[CounterOfpbarProcess] = theLossTable ;
      CounterOfpbarProcess++;
    }
  }
 
  BuildLambdaTable(aParticleType) ;

  BuildDEDXTable(aParticleType) ;

  if(&aParticleType == G4Proton::Proton())
    PrintInfoDefinition();

}

void G4hIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
  // cuts for  electron ....................
  DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;

  G4double LowEdgeEnergy , ionloss ;
  G4double  RateMass ;
  G4bool isOutRange ;
  static const G4MaterialTable* theMaterialTable=
                                   G4Material::GetMaterialTable();
  const G4double twoln10 = 2.*log(10.) ;
  const G4double Factor = twopi_mc2_rcl2 ;
  const G4double bg2lim = 0.0169 , taulim = 8.4146e-3 ;

  RateMass = electron_mass_c2/proton_mass_c2 ;

  //  create table

  G4int numOfMaterials = theMaterialTable->length();

  if ( theLossTable) {
     theLossTable->clearAndDestroy();
     delete theLossTable;
  }
  theLossTable = new G4PhysicsTable(numOfMaterials);

  //  loop for materials

  for (G4int J=0; J<numOfMaterials; J++)
  {

    // create physics vector and fill it

    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);
 
    // get material parameters needed for the energy loss calculation

    G4double ElectronDensity,Eexc,Eexc2,Cden,Mden,Aden,X0den,X1den,taul ;
    G4double* ShellCorrectionVector;
   
    const G4Material* material= (*theMaterialTable)[J];

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

    DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;

    // some local variables -------------------
    G4double tau,tau0,Tmax,gamma,bg2,beta2,rcut,delta,x,sh ;

    // now comes the loop for the kinetic energy values*****************

    for (G4int i = 0 ; i < TotBin ; i++)
    {
      LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
      tau = LowEdgeEnergy/proton_mass_c2 ;

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

      aVector->PutValue(i,ionloss) ;

    }
    theLossTable->insert(aVector);
  }

}

void G4hIonisation::BuildLambdaTable(const G4ParticleDefinition& aParticleType)
{
  // Build mean free path tables for the delta ray production process
  //     tables are built for MATERIALS 

      G4double LowEdgeEnergy , Value ,sigma ;
      G4bool isOutRange ;
      const G4MaterialTable* theMaterialTable=
                                         G4Material::GetMaterialTable();

      //create table

      G4int numOfMaterials = theMaterialTable->length();

      if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
      }

      theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);

  // get electron and particle cuts in kinetic energy

      DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
 
  // loop for materials 

      for (G4int J=0 ; J < numOfMaterials; J++)
      { 
        //create physics vector then fill it ....

        G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
               LowestKineticEnergy, HighestKineticEnergy, TotBin);

  // compute the (macroscopic) cross section first
 
        const G4Material* material= (*theMaterialTable)[J];
        
        const G4ElementVector* theElementVector=
                         material->GetElementVector() ;
        const G4double* theAtomicNumDensityVector =
                         material->GetAtomicNumDensityVector();
        const G4int NumberOfElements=
                         material->GetNumberOfElements() ;
 
  // get the electron kinetic energy cut for the actual material,
  //  it will be used in ComputeMicroscopicCrossSection
  // ( it is the SAME for ALL the ELEMENTS in THIS MATERIAL )
  //   ------------------------------------------------------

        DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;

        for ( G4int i = 0 ; i < TotBin ; i++ )
        {
           LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;

           sigma = 0. ;
           
           for (G4int iel=0; iel<NumberOfElements; iel++ )
           {
                sigma +=  theAtomicNumDensityVector[iel]*
                        ComputeMicroscopicCrossSection(aParticleType,
                          LowEdgeEnergy,
                          (*theElementVector)(iel)->GetZ() ) ;
           }

  // mean free path = 1./macroscopic cross section

           Value = sigma<=0 ? DBL_MAX : 1./sigma ;     

           aVector->PutValue(i, Value) ;
        }


        theMeanFreePathTable->insert(aVector);
      }
}


G4double G4hIonisation::ComputeMicroscopicCrossSection(
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber)
{
  //******************************************************************
  // cross section formula is OK for spin=0 and 1/2 only !
  // *****************************************************************

  // calculates the microscopic cross section in GEANT4 internal units
  //    ( it is called for elements , AtomicNumber = Z )

    G4double TotalEnergy,
             betasquare,
             MaxKineticEnergyTransfer,TotalCrossSection,tempvar;

    // get particle data ...................................

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

  // +term for spin=1/2 particle
     if(aParticleType.GetPDGSpin() == 0.5)
     {
       TotalCrossSection +=  0.5
                       *(MaxKineticEnergyTransfer-DeltaCutInKineticEnergyNow)
                       /(TotalEnergy*TotalEnergy);
       G4double  Section =  0.5
                       *(MaxKineticEnergyTransfer-DeltaCutInKineticEnergyNow)
                       /(TotalEnergy*TotalEnergy);
     }
       TotalCrossSection = twopi_mc2_rcl2 * AtomicNumber
                           *TotalCrossSection/betasquare;
    }
    else
       TotalCrossSection= 0. ;

    return TotalCrossSection ;
}
 
 
 
G4VParticleChange* G4hIonisation::PostStepDoIt(
                                              const G4Track& trackData,   
                                              const G4Step& stepData)         
{
  // Units are expressed in GEANT4 internal units.

  const G4DynamicParticle* aParticle ;
  G4Material* aMaterial;
  G4double KineticEnergy,TotalEnergy,TotalMomentum,
           betasquare,MaxKineticEnergyTransfer,
           DeltaKineticEnergy,DeltaTotalMomentum,costheta,sintheta,phi,
           dirx,diry,dirz,finalKineticEnergy,finalPx,finalPy,finalPz,
           x,xc,te2,grej,Psquare,Esquare,summass,rate,grejc,finalMomentum ;

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;

  aParticle = trackData.GetDynamicParticle() ;

  ParticleMass=aParticle->GetDefinition()->GetPDGMass();
  KineticEnergy=aParticle->GetKineticEnergy();
  TotalEnergy=KineticEnergy + ParticleMass ;
  Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ;
  Esquare=TotalEnergy*TotalEnergy ;
  summass = ParticleMass + electron_mass_c2 ;    
  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection() ;

  //  get kinetic energy cut for the electron....
  DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[aMaterial->GetIndex()];

  // some kinematics......................

  betasquare=Psquare/Esquare ;
  MaxKineticEnergyTransfer = 2.*electron_mass_c2*Psquare
                      /(summass*summass+2.*electron_mass_c2*KineticEnergy);

  // sampling kinetic energy of the delta ray 

  if( MaxKineticEnergyTransfer <= DeltaCutInKineticEnergyNow )
  {
    // pathological case (it should not happen ,
    // there is no change at all).....

    // return &aParticleChange;
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }
  else
  {
   // normal case ......................................
      xc=DeltaCutInKineticEnergyNow/MaxKineticEnergyTransfer ;
      rate=MaxKineticEnergyTransfer/TotalEnergy ;

     if(aParticle->GetDefinition()->GetPDGSpin() == 1)     
       te2=0.5*rate*rate ;
     else
       te2=0. ;

   // sampling follows ...
     grejc=1.-betasquare*xc+te2*xc*xc ;

     do {
          x=xc/(1.-(1.-xc)*G4UniformRand());
          grej=(1.-x*(betasquare-x*te2))/grejc ;
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
   G4double Edep = 0 ;

   if (finalKineticEnergy > MinKineticEnergy)
     {
      finalPx = TotalMomentum*ParticleDirection.x()
                        - DeltaTotalMomentum*DeltaDirection.x();
      finalPy = TotalMomentum*ParticleDirection.y()
                        - DeltaTotalMomentum*DeltaDirection.y();
      finalPz = TotalMomentum*ParticleDirection.z()
                        - DeltaTotalMomentum*DeltaDirection.z();
      finalMomentum =
                sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz) ;
      finalPx /= finalMomentum ;
      finalPy /= finalMomentum ;
      finalPz /= finalMomentum ;

      aParticleChange.SetMomentumChange( finalPx,finalPy,finalPz );
     }
   else
     {
       finalKineticEnergy = 0. ;
       Edep = finalKineticEnergy ;
       if (aParticle->GetDefinition()->GetParticleName() == "proton")
             aParticleChange.SetStatusChange(fStopAndKill);
       else  aParticleChange.SetStatusChange(fStopButAlive);
     }

   aParticleChange.SetEnergyChange( finalKineticEnergy );
   aParticleChange.SetNumberOfSecondaries(1);   
   aParticleChange.AddSecondary( theDeltaRay );
   aParticleChange.SetLocalEnergyDeposit (Edep);
      
   //ResetNumberOfInteractionLengthLeft();
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

void G4hIonisation::PrintInfoDefinition()
{
  G4String comments = "  Knock-on electron cross sections . ";
           comments += "\n         Good description above the mean excitation energy.\n";
           comments += "         delta ray energy sampled from  differential Xsection.";

  G4cout << endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
                                                  "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins. \n";
}

