// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuIonisation.cc,v 1.6 2000-02-10 08:32:21 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4MuIonisation physics process -------------
//                by Laszlo Urban, September 1997
// ------------------------------------------------------------------
// It is the implementation of the NEW IONISATION PROCESS.
// It calculates the ionisation of muons.
// **************************************************************
// 08-04-98: remove 'tracking cut' of the ionizing particle, MMa
// 26/10/98: new stuff from R.Kokoulin + cleanup , L.Urban
// --------------------------------------------------------------
 

#include "G4MuIonisation.hh"
#include "G4UnitsTable.hh"

#include "G4ios.hh"

// constructor and destructor
 
G4MuIonisation::G4MuIonisation(const G4String& processName)
   : G4MuEnergyLoss(processName),
     theMeanFreePathTable(NULL),
     LowerBoundLambda(1.*keV),
     UpperBoundLambda(10000.*TeV),
     NbinLambda(100)
{  }

G4MuIonisation::~G4MuIonisation() 
{
     if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
     }

}

void G4MuIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
//  just call BuildLossTable+BuildLambdaTable
{
    // get bining from EnergyLoss
    LowestKineticEnergy  = GetLowerBoundEloss() ;
    HighestKineticEnergy = GetUpperBoundEloss() ;
    TotBin               = GetNbinEloss() ;

  BuildLossTable(aParticleType) ;
  G4double Charge = aParticleType.GetPDGCharge();     

  if(Charge>0.)
  {
    RecorderOfmuplusProcess[CounterOfmuplusProcess] = (*this).theLossTable ;
    CounterOfmuplusProcess++;
  }
  else
  {
    RecorderOfmuminusProcess[CounterOfmuminusProcess] = (*this).theLossTable ;
    CounterOfmuminusProcess++;
  }
 
  G4double electronCutInRange = G4Electron::Electron()->GetCuts();  
  if(electronCutInRange != lastelectronCutInRange)
    BuildLambdaTable(aParticleType) ;
 
   G4MuEnergyLoss::BuildDEDXTable(aParticleType) ;

  if(&aParticleType == theMuonPlus)  
     PrintInfoDefinition() ;
}

void G4MuIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
  DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;

  G4double LowEdgeEnergy , ionloss ;
  G4double deltaloss ;
  G4double RateMass ;
  G4bool isOutRange ;
  static const G4MaterialTable* theMaterialTable=
                                   G4Material::GetMaterialTable();
  const G4double twoln10 = 2.*log(10.) ;
  const G4double Factor = twopi_mc2_rcl2 ;
  const G4double bg2lim = 0.0169 , taulim = 8.4146e-3 ;

  ParticleMass = aParticleType.GetPDGMass() ;
  RateMass = electron_mass_c2/ParticleMass ;

  G4int numOfMaterials = theMaterialTable->length();

  if ( theLossTable) {
     theLossTable->clearAndDestroy();
     delete theLossTable;
  }
  theLossTable = new G4PhysicsTable(numOfMaterials);

  for (G4int J=0; J<numOfMaterials; J++)
  {
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);
  
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
    ShellCorrectionVector = material->GetIonisation()
                                    ->GetShellCorrectionVector();

    const G4ElementVector* theElementVector=
                   material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector=
                   material->GetAtomicNumDensityVector() ;
    const G4int NumberOfElements=
                   material->GetNumberOfElements() ;
    DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;

    G4double tau,tau0,Tmax,gamma,bg2,beta2,rcut,delta,x,sh ;
    for (G4int i = 0 ; i < TotBin ; i++)
    {
      LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
      tau = LowEdgeEnergy/ParticleMass ;
      gamma = tau +1. ;
      bg2 = tau*(tau+2.) ;
      beta2 = bg2/(gamma*gamma) ;
      Tmax = 2.*electron_mass_c2*bg2
             /(1.+2.*gamma*RateMass+RateMass*RateMass) ;

      if ( tau < taul )
      //  low energy part , parametrized energy loss formulae
      {
        ionloss = 0. ;
        deltaloss = 0. ;

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
        if ( DeltaCutInKineticEnergyNow < Tmax)
        {
          deltaloss = log(Tmax/DeltaCutInKineticEnergyNow)-
                      beta2*(1.-DeltaCutInKineticEnergyNow/Tmax) ;
          if(aParticleType.GetPDGSpin() == 0.5)
            deltaloss += 0.25*(Tmax-DeltaCutInKineticEnergyNow)*
                              (Tmax-DeltaCutInKineticEnergyNow)/
                        (LowEdgeEnergy*LowEdgeEnergy+proton_mass_c2*proton_mass_c2) ;
            deltaloss *= Factor*ElectronDensity/beta2 ;
        }
        ionloss -= deltaloss ;
      }
      else
      // high energy part , Bethe-Bloch formula 
      {

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

        ionloss -= delta + sh ;
        ionloss /= beta2 ;
         
        // correction of R. Kokoulin  
        G4double E = LowEdgeEnergy+ParticleMass ;
        G4double epmax = RateMass*E*E/(RateMass*E+ParticleMass) ;
        G4double apar = log(2.*epmax/electron_mass_c2) ;
        ionloss += fine_structure_const*(log(2.*E/ParticleMass)-apar/3.)*
                                        apar*apar/twopi ; 

        ionloss *= Factor*ElectronDensity ;
      }
      if ( ionloss <= 0.)
        ionloss = 0. ;
      aVector->PutValue(i,ionloss) ;
    }
    theLossTable->insert(aVector);
  }
}

void G4MuIonisation::BuildLambdaTable(const G4ParticleDefinition& aParticleType)
{
  // Build mean free path tables for the delta ray production process
  G4double LowEdgeEnergy,Tmax , Value ,sigma ;
  G4bool isOutRange ;
  const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();

  G4int numOfMaterials = theMaterialTable->length();

  if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
     }

  theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);

  // get electron and particle cuts in kinetic energy
  DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;

  for (G4int J=0 ; J < numOfMaterials; J++)
  { 
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
               LowerBoundLambda,UpperBoundLambda,NbinLambda);

    const G4Material* material= (*theMaterialTable)[J];
    const G4ElementVector* theElementVector=
                         material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector =
                         material->GetAtomicNumDensityVector();
    const G4int NumberOfElements=
                         material->GetNumberOfElements() ;
    DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;

    for ( G4int i = 0 ; i < NbinLambda ; i++ )
    {
      LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
      sigma = 0. ;

      // check threshold here !
      G4double Tmax = 2.*electron_mass_c2*LowEdgeEnergy*
                     (LowEdgeEnergy+2.*ParticleMass)/
                     (ParticleMass*ParticleMass+2.*electron_mass_c2*
                     (LowEdgeEnergy+ParticleMass)+
                       electron_mass_c2*electron_mass_c2) ;

      if(Tmax > DeltaCutInKineticEnergyNow)
      {
        for (G4int iel=0; iel<NumberOfElements; iel++ )
        {
          sigma +=  theAtomicNumDensityVector[iel]*
                      ComputeMicroscopicCrossSection(aParticleType,
                          LowEdgeEnergy,
                          (*theElementVector)(iel)->GetZ() ) ;
        }
      }

      Value = sigma<=0 ? DBL_MAX : 1./sigma ;     

      aVector->PutValue(i, Value) ;
    }

    theMeanFreePathTable->insert(aVector);
  }
}


G4double G4MuIonisation::ComputeMicroscopicCrossSection(
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber)
{
  const G4double xgi[] = {0.06943,0.33001,0.66999,0.93057} ;
  const G4double wgi[] = {0.17393,0.32607,0.32607,0.17393} ;
  const G4double ak1 = 4.6 ;
  const G4int k2 = 2 ;
  const G4double masspar = 0.5*ParticleMass*ParticleMass/electron_mass_c2 ;
  G4double TotalEnergy=KineticEnergy + ParticleMass;
  G4double KnockonMaxEnergy = TotalEnergy/(1.+masspar/TotalEnergy) ; 

  G4double TotalCrossSection= 0. ; 

  if( KnockonMaxEnergy > DeltaCutInKineticEnergyNow )
  {
    G4double aaa = log(DeltaCutInKineticEnergyNow);
    G4double bbb = log(KnockonMaxEnergy) ;
    G4int kkk = int((bbb-aaa)/ak1)+k2 ;
    G4double hhh = (bbb-aaa)/kkk ;
    G4double step = exp(hhh) ;
    G4double ymax = 1./KnockonMaxEnergy ;
      
    for (G4int k=0; k<kkk; k++)
    {
      G4double ymin = ymax ;
      ymax = ymin*step ;
      G4double hhy = ymax-ymin ;
      for (G4int i=0; i<4; i++)
      {
        G4double y = ymin+hhy*xgi[i];
        G4double ep = 1./y ;
        TotalCrossSection += ep*ep*wgi[i]*hhy*
                             ComputeDMicroscopicCrossSection(
                             aParticleType,KineticEnergy,
                             AtomicNumber,ep) ;
      }
    }
  }
  return TotalCrossSection ;
}
 
G4double G4MuIonisation::ComputeDMicroscopicCrossSection(
                                 const G4ParticleDefinition& ParticleType,
                                 G4double KineticEnergy, G4double AtomicNumber,
                                 G4double KnockonEnergy)
 // Calculates the differential (D) microscopic cross section
 //   using the cross section formula of R.P. Kokoulin (10/98)
{
  const G4double masspar=0.5*ParticleMass*ParticleMass/electron_mass_c2 ;
  const G4double alphaprime = fine_structure_const/twopi ;
  G4double TotalEnergy = KineticEnergy + ParticleMass ;
  G4double KnockonMaxEnergy = TotalEnergy/(1.+masspar/TotalEnergy) ;

  G4double DCrossSection = 0. ;

  if(KnockonEnergy >=  KnockonMaxEnergy)  return DCrossSection ;

  G4double v = KnockonEnergy/TotalEnergy ;
  DCrossSection = twopi_mc2_rcl2*AtomicNumber*
                 (1.-KnockonEnergy/KnockonMaxEnergy+0.5*v*v)/
                 (KnockonEnergy*KnockonEnergy) ;
  G4double a1 = log(1.+2.*KnockonEnergy/electron_mass_c2) ;
  G4double a3 = log(4.*TotalEnergy*(TotalEnergy-KnockonEnergy)/
                    (ParticleMass*ParticleMass)) ;
   DCrossSection *= (1.+alphaprime*a1*(a3-a1)) ; 

  return DCrossSection ;
}
 
G4VParticleChange* G4MuIonisation::PostStepDoIt(
                                              const G4Track& trackData,   
                                              const G4Step& stepData)         
{
  const G4DynamicParticle* aParticle ;
  const G4double alphaprime = fine_structure_const/twopi ;
  G4Material* aMaterial;
  G4double KineticEnergy,TotalEnergy,TotalMomentum,
           betasquare,MaxKineticEnergyTransfer,
           DeltaKineticEnergy,DeltaTotalMomentum,costheta,sintheta,phi,
           dirx,diry,dirz,finalKineticEnergy,finalPx,finalPy,finalPz,
           x,xc,te2,grej,Psquare,Esquare,summass,rate,grejc,finalMomentum ;
  G4double Charge ;

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;
  aParticle = trackData.GetDynamicParticle() ;
  Charge=aParticle->GetDefinition()->GetPDGCharge();
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
    // pathological case (it should not happen ,
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }
  else
  {
   // normal case ......................................
    xc=DeltaCutInKineticEnergyNow/MaxKineticEnergyTransfer ;
    rate=MaxKineticEnergyTransfer/TotalEnergy ;
    te2=0.5*rate*rate ;

   // sampling follows ...
    G4double a0=log(2.*TotalEnergy/ParticleMass) ;
    grejc=(1.-betasquare*xc+te2*xc*xc)*
            (1.+ alphaprime*a0*a0) ;
    do {
        x=xc/(1.-(1.-xc)*G4UniformRand());
        G4double twoep = 2.*x*MaxKineticEnergyTransfer ;
        grej=(1.-x*(betasquare-x*te2))*
             (1.+alphaprime*log(1.+twoep/electron_mass_c2)*
             (a0+log((2.*TotalEnergy-twoep)/ParticleMass)-
              log(1.+twoep/electron_mass_c2)))
              /grejc ;
       
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
   finalKineticEnergy = KineticEnergy - DeltaKineticEnergy ;
   if (finalKineticEnergy > 0. )
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
       aParticleChange.SetStatusChange(fStopButAlive);
     }

   aParticleChange.SetEnergyChange( finalKineticEnergy );
   aParticleChange.SetNumberOfSecondaries(1);
   aParticleChange.AddSecondary( theDeltaRay );
   aParticleChange.SetLocalEnergyDeposit (0.);
      
   return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

void G4MuIonisation::PrintInfoDefinition()
{
  G4String comments = "knock-on electron cross sections .\n ";
           comments += "         Good description above the mean excitation energy.\n";
           comments += "         delta ray energy sampled from  differential Xsection." ;

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n      PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
                                               "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins. \n";
}

