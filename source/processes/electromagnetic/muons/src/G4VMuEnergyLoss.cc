// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VMuEnergyLoss.cc,v 1.5 2000-06-13 16:36:25 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4VMuEnergyLoss physics process -----------
//                by Laszlo Urban, September 1997
// **************************************************************
// It is the implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of muons.
// **************************************************************
//
// corrections by L.Urban on 27/05/98 (other corrs come soon!)
// cleanup  L.Urban on 23/10/98
// corrections due to new e.m. structure L.Urban 10/02/00
// --------------------------------------------------------------
 

#include "G4VMuEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4Poisson.hh"

// Initialisation of static members *******************************************

G4int       G4VMuEnergyLoss::NbOfProcesses           = 3 ;
G4PhysicsTable** G4VMuEnergyLoss::RecorderOfmuplusProcess  =
                                           new G4PhysicsTable*[10] ;
G4PhysicsTable** G4VMuEnergyLoss::RecorderOfmuminusProcess =
                                           new G4PhysicsTable*[10] ;
G4int       G4VMuEnergyLoss::CounterOfmuplusProcess  = 0 ;
G4int       G4VMuEnergyLoss::CounterOfmuminusProcess = 0 ;


G4PhysicsTable* G4VMuEnergyLoss::theDEDXmuplusTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::theRangemuplusTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::theInverseRangemuplusTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::theLabTimemuplusTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::theProperTimemuplusTable = NULL ;

G4PhysicsTable* G4VMuEnergyLoss::theDEDXmuminusTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::theRangemuminusTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::theInverseRangemuminusTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::theLabTimemuminusTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::theProperTimemuminusTable = NULL ;

G4PhysicsTable* G4VMuEnergyLoss::themuplusRangeCoeffATable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::themuplusRangeCoeffBTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::themuplusRangeCoeffCTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::themuminusRangeCoeffATable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::themuminusRangeCoeffBTable = NULL ;
G4PhysicsTable* G4VMuEnergyLoss::themuminusRangeCoeffCTable = NULL ;
 
G4double G4VMuEnergyLoss::LowerBoundEloss = 1.*keV ;
G4double G4VMuEnergyLoss::UpperBoundEloss = 1000000.*TeV ;
G4int    G4VMuEnergyLoss::NbinEloss = 150 ;
G4double G4VMuEnergyLoss::RTable,G4VMuEnergyLoss::LOGRTable;

G4EnergyLossMessenger* G4VMuEnergyLoss::eLossMessenger = NULL ;

// constructor and destructor
 
G4VMuEnergyLoss::G4VMuEnergyLoss(const G4String& processName)
   : G4VEnergyLoss (processName),
     theLossTable(NULL),
     theRangeCoeffATable(NULL),
     theRangeCoeffBTable(NULL),
     theRangeCoeffCTable(NULL),
     lastgammaCutInRange(0.),
     lastelectronCutInRange(0.),
     theElectron ( G4Electron::Electron() ),
     thePositron ( G4Positron::Positron() ),
     theMuonPlus ( G4MuonPlus::MuonPlus() ),
     theMuonMinus ( G4MuonMinus::MuonMinus() )
{ 
}

G4VMuEnergyLoss::~G4VMuEnergyLoss() 
{
     if(theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable; theLossTable = NULL; 
     }

}
 
   void G4VMuEnergyLoss::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{

  //  calculate data members LOGRTable,RTable first
  G4double lrate = log(UpperBoundEloss/LowerBoundEloss);
  LOGRTable=lrate/NbinEloss;
  RTable   =exp(LOGRTable);

  G4bool MakeTable ;
  ParticleMass = aParticleType.GetPDGMass() ; 
  G4double Charge = aParticleType.GetPDGCharge()/eplus ;
  G4double gammaCutInRange = G4Gamma::Gamma()->GetCuts(); 
  G4double electronCutInRange = G4Electron::Electron()->GetCuts(); 

  MakeTable = false ;
  // Create tables only if there are new cut values 
  if((gammaCutInRange == lastgammaCutInRange) &&
     (electronCutInRange == lastelectronCutInRange))
  {
     ;
  }
  else
  {
    if((Charge > 0.)&&(CounterOfmuplusProcess==NbOfProcesses))
       MakeTable = true ;
    if((Charge < 0.)&&(CounterOfmuminusProcess==NbOfProcesses))
       MakeTable = true ;
  }
      
  if( MakeTable )
  {
   // Build energy loss table as a sum of the energy loss due to the
   //           different processes.                                           
    const G4MaterialTable* theMaterialTable=
                                     G4Material::GetMaterialTable();

    G4int numOfMaterials = theMaterialTable->length();

    if( Charge >0.)    
    {
      RecorderOfProcess=RecorderOfmuplusProcess;
      CounterOfProcess=CounterOfmuplusProcess;

      if(CounterOfProcess == NbOfProcesses)
      {
        if(theDEDXmuplusTable)
        { theDEDXmuplusTable->clearAndDestroy();
          delete theDEDXmuplusTable; }

        theDEDXmuplusTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXmuplusTable;
      }
    }
    else
    {
      RecorderOfProcess=RecorderOfmuminusProcess;
      CounterOfProcess=CounterOfmuminusProcess;

      if(CounterOfProcess == NbOfProcesses)
      {
        if(theDEDXmuminusTable)
        { theDEDXmuminusTable->clearAndDestroy();
          delete theDEDXmuminusTable; }
        theDEDXmuminusTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXmuminusTable;
      }
    }

    if(CounterOfProcess == NbOfProcesses)
    {
      //  loop for materials
      G4double LowEdgeEnergy , Value ;
      G4bool isOutRange ;
      G4int J;
      G4PhysicsTable* pointer ;

      for (J=0; J<numOfMaterials; J++)
      {
        // create physics vector and fill it
        G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowerBoundEloss, UpperBoundEloss, NbinEloss);   
        // loop for the kinetic energy
        for (G4int i=0; i<NbinEloss; i++)
        {
          LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;      
          Value = 0. ;
          for (G4int process=0; process < NbOfProcesses; process++)
          {
            pointer= RecorderOfProcess[process];
            Value += (*pointer)[J]->
                               GetValue(LowEdgeEnergy,isOutRange) ;
          }
          aVector->PutValue(i,Value) ; 
        }
        theDEDXTable->insert(aVector) ;
      }
 }
      //  reset counter to zero ..................
      if( Charge >0.)    
        CounterOfmuplusProcess=0 ;
      else
        CounterOfmuminusProcess=0 ;

      ParticleMass = aParticleType.GetPDGMass() ;

      if(Charge > 0.)
      {
        // Build range table
        theRangemuplusTable = BuildRangeTable(
                  theDEDXmuplusTable,theRangemuplusTable,  
                  LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build lab/proper time tables
        theLabTimemuplusTable = BuildLabTimeTable(theDEDXmuplusTable,
                          theLabTimemuplusTable,
                          LowerBoundEloss,UpperBoundEloss,NbinEloss);
        theProperTimemuplusTable = BuildProperTimeTable(theDEDXmuplusTable,
                          theProperTimemuplusTable,
                          LowerBoundEloss,UpperBoundEloss,NbinEloss);

      // Build coeff tables for the energy loss calculation
        themuplusRangeCoeffATable = BuildRangeCoeffATable(theRangemuplusTable,
                              themuplusRangeCoeffATable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);

        themuplusRangeCoeffBTable = BuildRangeCoeffBTable(theRangemuplusTable,
                              themuplusRangeCoeffBTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);

        themuplusRangeCoeffCTable = BuildRangeCoeffCTable(theRangemuplusTable,
                              themuplusRangeCoeffCTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);

        // invert the range table
        theInverseRangemuplusTable = BuildInverseRangeTable(theRangemuplusTable,
                               themuplusRangeCoeffATable,
                               themuplusRangeCoeffBTable,
                               themuplusRangeCoeffCTable,
                               theInverseRangemuplusTable,
                               LowerBoundEloss,UpperBoundEloss,NbinEloss);
      }
      else
      {
        // Build range table
        theRangemuminusTable = BuildRangeTable(
                  theDEDXmuminusTable,theRangemuminusTable,  
                  LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build lab/proper time tables
        theLabTimemuminusTable = BuildLabTimeTable(theDEDXmuminusTable,
                          theLabTimemuminusTable,
                          LowerBoundEloss,UpperBoundEloss,NbinEloss);
        theProperTimemuminusTable = BuildProperTimeTable(theDEDXmuminusTable,
                          theProperTimemuminusTable,
                          LowerBoundEloss,UpperBoundEloss,NbinEloss);

      // Build coeff tables for the energy loss calculation
        themuminusRangeCoeffATable = BuildRangeCoeffATable(theRangemuminusTable,
                              themuminusRangeCoeffATable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);

        themuminusRangeCoeffBTable = BuildRangeCoeffBTable(theRangemuminusTable,
                              themuminusRangeCoeffBTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);

        themuminusRangeCoeffCTable = BuildRangeCoeffCTable(theRangemuminusTable,
                              themuminusRangeCoeffCTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);

        // invert the range table
        theInverseRangemuminusTable = BuildInverseRangeTable(theRangemuminusTable,
                               themuminusRangeCoeffATable,
                               themuminusRangeCoeffBTable,
                               themuminusRangeCoeffCTable,
                               theInverseRangemuminusTable,
                               LowerBoundEloss,UpperBoundEloss,NbinEloss);
        // invert the range table
      }

    }

    // make the energy loss and the range table available
    G4EnergyLossTables::Register(&aParticleType,  
      (Charge > 0)? theDEDXmuplusTable: theDEDXmuminusTable,
      (Charge > 0)? theRangemuplusTable: theRangemuminusTable,
      (Charge > 0)? theInverseRangemuplusTable: theInverseRangemuminusTable,
      (Charge > 0)? theLabTimemuplusTable: theLabTimemuminusTable,
      (Charge > 0)? theProperTimemuplusTable: theProperTimemuminusTable,
      LowerBoundEloss, UpperBoundEloss, 1.,NbinEloss);

    lastgammaCutInRange = gammaCutInRange ;
    lastelectronCutInRange = electronCutInRange ;

}
      

G4double G4VMuEnergyLoss::GetConstraints(const G4DynamicParticle *aParticle,
                                                      G4Material *aMaterial)
{

  // returns the Step limit
  // dRoverRange is the max. allowed relative range loss in one Step
  // it calculates dEdx and the range as well....

  G4double KineticEnergy,StepLimit;
  G4bool isOutRange ;
  G4int index,bin ;

  if(aParticle->GetDefinition()->GetPDGCharge()>0.)
  {
    theDEDXTable = theDEDXmuplusTable ;
    theRangeTable = theRangemuplusTable ;
    theRangeCoeffATable=themuplusRangeCoeffATable ;
    theRangeCoeffBTable=themuplusRangeCoeffBTable ;
    theRangeCoeffCTable=themuplusRangeCoeffCTable ;
  }
  else
  {
    theDEDXTable = theDEDXmuminusTable ;
    theRangeTable = theRangemuminusTable ;
    theRangeCoeffATable=themuminusRangeCoeffATable ;
    theRangeCoeffBTable=themuminusRangeCoeffBTable ;
    theRangeCoeffCTable=themuminusRangeCoeffCTable ;
  }


  G4double Thigh = UpperBoundEloss/RTable ;
  KineticEnergy = aParticle->GetKineticEnergy();

  bin = G4int(log(KineticEnergy/LowerBoundEloss)/LOGRTable) ;
  EnergyBinNumber = bin ;

  index = aMaterial->GetIndex() ;

  if( KineticEnergy < LowerBoundEloss )
    {
      fdEdx = sqrt(KineticEnergy/LowerBoundEloss)*
             (*theDEDXTable)(index)->GetValue(LowerBoundEloss,isOutRange) ;

      fRangeNow = sqrt(KineticEnergy/LowerBoundEloss)*
             (*theRangeTable)(index)->GetValue(LowerBoundEloss,isOutRange) ;

      StepLimit = fRangeNow ;
    }
  else if (KineticEnergy > Thigh )
    {
      fdEdx = (*theDEDXTable)(index)->GetValue(Thigh,isOutRange);
      fRangeNow = (*theRangeTable)(index)->GetValue(Thigh,isOutRange);
      if (fdEdx > 0.) fRangeNow += (KineticEnergy-Thigh)/fdEdx;
      StepLimit = c1lim*fRangeNow;
    }
  else
    {
         fdEdx = (*theDEDXTable)(index)->
                            GetValue(KineticEnergy,isOutRange) ;

         RangeCoeffA = (*(*theRangeCoeffATable)(index))(EnergyBinNumber) ;
         RangeCoeffB = (*(*theRangeCoeffBTable)(index))(EnergyBinNumber) ;
         RangeCoeffC = (*(*theRangeCoeffCTable)(index))(EnergyBinNumber) ;

         fRangeNow = (RangeCoeffA*KineticEnergy+RangeCoeffB)
                    *KineticEnergy+RangeCoeffC ;

  //   compute the (random) Step limit ..............
     if(fRangeNow>finalRange)
       {
         StepLimit = c1lim*fRangeNow+c2lim+c3lim/fRangeNow ;

        //  randomise this value
         if(rndmStepFlag) StepLimit = finalRange+(StepLimit-finalRange)*
                                                           G4UniformRand() ;
            if(StepLimit > fRangeNow) StepLimit = fRangeNow ;
       }
       else
         StepLimit = fRangeNow ;
  }

  return StepLimit ;

}

G4VParticleChange* G4VMuEnergyLoss::AlongStepDoIt( 
                              const G4Track& trackData,const G4Step& stepData) 
{
 // compute the energy loss after a Step

  static const G4double faclow = 1.5 ;

  // get particle and material pointers from trackData
  const G4DynamicParticle* aParticle = trackData.GetDynamicParticle();
  G4double E      = aParticle->GetKineticEnergy() ;
  G4double charge = aParticle->GetDefinition()->GetPDGCharge();
 
  G4Material* aMaterial = trackData.GetMaterial();
  G4int index = aMaterial->GetIndex();
 
  G4double Step = stepData.GetStepLength();
 
  aParticleChange.Initialize(trackData);
 
  // do not track further if kin.energy < 1. eV
   const G4double MinKineticEnergy = 1.*eV;
   const G4double linLossLimit = 0.05 ;
  
  G4double MeanLoss, finalT;
 
  if (E < MinKineticEnergy)    finalT = 0.; 
  else if ( E< faclow*LowerBoundEloss)
  {
    if (Step >= fRangeNow)  finalT = 0.;
    else finalT = E*(1.-Step/fRangeNow) ;
  }
   
  else if (E>=UpperBoundEloss) finalT = E - Step*fdEdx;

  else if (Step >= fRangeNow)  finalT = 0.;
 
  else
  {
    if(Step/fRangeNow < linLossLimit) finalT = E-Step*fdEdx ;
    else
    {
       if (charge<0.) finalT = G4EnergyLossTables::GetPreciseEnergyFromRange(
                               theMuonMinus,fRangeNow-Step,aMaterial);
       else           finalT = G4EnergyLossTables::GetPreciseEnergyFromRange(
                               theMuonPlus,fRangeNow-Step,aMaterial);
     }
  }

  if(finalT < MinKineticEnergy) finalT = 0. ;

  MeanLoss = E-finalT ; 

  //now the loss with fluctuation
  if ((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowerBoundEloss))

  {
    finalT = E-GetLossWithFluct(aParticle,aMaterial,MeanLoss);
    if (finalT < 0.) finalT = 0. ;
  }

  // kill the particle if the kinetic energy <= 0
  if (finalT <= 0. )
    {
      finalT = 0.;
      aParticleChange.SetStatusChange(fStopButAlive);
    }

  aParticleChange.SetNumberOfSecondaries(0);
  aParticleChange.SetEnergyChange(finalT);
  aParticleChange.SetLocalEnergyDeposit(E-finalT);

  return &aParticleChange;
}



