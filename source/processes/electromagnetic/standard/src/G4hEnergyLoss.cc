// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hEnergyLoss.cc,v 1.19 2000-04-10 09:55:05 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hEnergyLoss physics process -----------
//                by Laszlo Urban, 30 May 1997 
//
// **************************************************************
// It is the first implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of charged hadrons.
// **************************************************************
//
// 7/10/98: bug fixes + some cleanup , L.Urban 
// 22/10/98 : cleanup , L.Urban
// 07/12/98 : works for ions as well+ bug corrected, L.Urban
// 02/02/99 : several bugs fixed, L.Urban
// 10/02/00  modifications , new e.m. structure, L.Urban
// --------------------------------------------------------------

#include "G4hEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4Poisson.hh"

// Initialisation of static members ******************************************
G4int            G4hEnergyLoss::NbOfProcesses  = 1 ;

G4int            G4hEnergyLoss::CounterOfProcess = 0 ;
G4PhysicsTable** G4hEnergyLoss::RecorderOfProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4hEnergyLoss::CounterOfpProcess = 0 ;
G4PhysicsTable** G4hEnergyLoss::RecorderOfpProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4hEnergyLoss::CounterOfpbarProcess = 0 ;
G4PhysicsTable** G4hEnergyLoss::RecorderOfpbarProcess =
                                           new G4PhysicsTable*[10] ;

G4PhysicsTable* G4hEnergyLoss::theDEDXpTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::theDEDXpbarTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::theRangepTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::theRangepbarTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::theInverseRangepTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::theInverseRangepbarTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::theLabTimepTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::theLabTimepbarTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::theProperTimepTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::theProperTimepbarTable = NULL ;

G4PhysicsTable* G4hEnergyLoss::thepRangeCoeffATable = NULL ;
G4PhysicsTable* G4hEnergyLoss::thepRangeCoeffBTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::thepRangeCoeffCTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::thepbarRangeCoeffATable = NULL ;
G4PhysicsTable* G4hEnergyLoss::thepbarRangeCoeffBTable = NULL ;
G4PhysicsTable* G4hEnergyLoss::thepbarRangeCoeffCTable = NULL ;

G4PhysicsTable* G4hEnergyLoss::theDEDXTable = NULL ;

const G4Proton* G4hEnergyLoss::theProton=G4Proton::Proton() ;
const G4AntiProton* G4hEnergyLoss::theAntiProton=G4AntiProton::AntiProton() ;

G4double G4hEnergyLoss::ptableElectronCutInRange = 0.0*mm ;
G4double G4hEnergyLoss::pbartableElectronCutInRange = 0.0*mm ;

G4double         G4hEnergyLoss::Charge ;   

G4double G4hEnergyLoss::LowerBoundEloss = 1.*keV ;
G4double G4hEnergyLoss::UpperBoundEloss = 100.*TeV ;	
G4int    G4hEnergyLoss::NbinEloss = 100 ; 
G4double G4hEnergyLoss::RTable,G4hEnergyLoss::LOGRTable;
// just to keep hLowEnergyIonisation working
// ****************************************
G4double G4hEnergyLoss::LowestKineticEnergy ;
G4double G4hEnergyLoss::HighestKineticEnergy ;
G4int    G4hEnergyLoss::TotBin ; 
// ****************************************

// constructor and destructor
 
G4hEnergyLoss::G4hEnergyLoss(const G4String& processName)
   : G4VEnergyLoss (processName),
     theLossTable (NULL),
     MinKineticEnergy(1.*eV), 
     linLossLimit(0.05)
{
  // just to keep hLowEnergyIonisation working
  // ****************************************
  LowestKineticEnergy  = LowerBoundEloss ;
  HighestKineticEnergy = UpperBoundEloss ;
  TotBin               = NbinEloss       ;
  // ****************************************
}

G4hEnergyLoss::~G4hEnergyLoss() 
{
     if(theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable;
     }
}
 
void G4hEnergyLoss::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{

  //  calculate data members LOGRTable,RTable first
  G4double lrate = log(UpperBoundEloss/LowerBoundEloss);
  LOGRTable=lrate/NbinEloss;
  RTable   =exp(LOGRTable);

  // create table if there is no table or there is a new cut value
  G4bool MakeTable = false ;
     
  G4double ElectronCutInRange = G4Electron::Electron()->GetCuts();

  // create/fill proton or antiproton tables depending on the charge 
  Charge = aParticleType.GetPDGCharge()/eplus;
  ParticleMass = aParticleType.GetPDGMass() ;

  if (Charge>0.) {theDEDXTable= theDEDXpTable;}
  else           {theDEDXTable= theDEDXpbarTable;}

  if(
     ((Charge>0.) && ((theDEDXTable==NULL) || 
     (ElectronCutInRange != ptableElectronCutInRange)))
     ||  
     ((Charge<0.) && ((theDEDXTable==NULL) || 
     (ElectronCutInRange != pbartableElectronCutInRange)))
    )
      MakeTable = true ;
  
  const G4MaterialTable* theMaterialTable=
                                   G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();

  if( MakeTable )
  {

  // Build energy loss table as a sum of the energy loss due to the
  //              different processes.                                           
    if( Charge >0.)    
    {
      RecorderOfProcess=RecorderOfpProcess;
      CounterOfProcess=CounterOfpProcess;

      if(CounterOfProcess == NbOfProcesses)
      {
        if(theDEDXpTable)
        { theDEDXpTable->clearAndDestroy();
          delete theDEDXpTable; }
        theDEDXpTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXpTable;
        ptableElectronCutInRange = ElectronCutInRange ;
      }
    }
    else
    {
      RecorderOfProcess=RecorderOfpbarProcess;
      CounterOfProcess=CounterOfpbarProcess;

      if(CounterOfProcess == NbOfProcesses)
      {
        if(theDEDXpbarTable)
        { theDEDXpbarTable->clearAndDestroy();
          delete theDEDXpbarTable; }
        theDEDXpbarTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXpbarTable;
        pbartableElectronCutInRange = ElectronCutInRange ;
      }
    }

    if(CounterOfProcess == NbOfProcesses)
    {
      //  loop for materials
      G4double LowEdgeEnergy , Value ;
      G4bool isOutRange ;
      G4PhysicsTable* pointer ;

      for (G4int J=0; J<numOfMaterials; J++)
      { 
        // create physics vector and fill it
        G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowerBoundEloss, UpperBoundEloss, NbinEloss);   

        // loop for the kinetic energy
        for (G4int i=0; i<NbinEloss; i++)
        {
          LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;      
          Value = 0. ;
    
          // loop for the contributing processes
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
      
      //  reset counter to zero ..................
      if( Charge >0.)    
        CounterOfpProcess=0 ;
      else
        CounterOfpbarProcess=0 ;

     // ParticleMass = aParticleType.GetPDGMass() ;
      ParticleMass = proton_mass_c2 ;

      if(Charge > 0.)
      {
       // Build range table
       theRangepTable = BuildRangeTable(theDEDXpTable,
                        theRangepTable,
                        LowerBoundEloss,UpperBoundEloss,NbinEloss);
       // Build lab/proper time tables
       theLabTimepTable = BuildLabTimeTable(theDEDXpTable,
                         theLabTimepTable,
                         LowerBoundEloss,UpperBoundEloss,NbinEloss);
       theProperTimepTable = BuildProperTimeTable(theDEDXpTable,
                            theProperTimepTable,
                            LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build coeff tables for the energy loss calculation
       thepRangeCoeffATable = BuildRangeCoeffATable(theRangepTable,
                             thepRangeCoeffATable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepRangeCoeffBTable = BuildRangeCoeffBTable(theRangepTable,
                             thepRangeCoeffBTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepRangeCoeffCTable = BuildRangeCoeffCTable(theRangepTable,
                             thepRangeCoeffCTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // invert the range table
       theInverseRangepTable = BuildInverseRangeTable(theRangepTable,
                              thepRangeCoeffATable,
                              thepRangeCoeffBTable,
                              thepRangeCoeffCTable,
                              theInverseRangepTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);
  
      }
      else
      {
       // Build range table
       theRangepbarTable = BuildRangeTable(theDEDXpbarTable,
                        theRangepbarTable,
                        LowerBoundEloss,UpperBoundEloss,NbinEloss);
       // Build lab/proper time tables
       theLabTimepbarTable = BuildLabTimeTable(theDEDXpbarTable,
                         theLabTimepbarTable,
                         LowerBoundEloss,UpperBoundEloss,NbinEloss);
       theProperTimepbarTable = BuildProperTimeTable(theDEDXpbarTable,
                            theProperTimepbarTable,
                            LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build coeff tables for the energy loss calculation
       thepbarRangeCoeffATable = BuildRangeCoeffATable(theRangepbarTable,
                             thepbarRangeCoeffATable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepbarRangeCoeffBTable = BuildRangeCoeffBTable(theRangepbarTable,
                             thepbarRangeCoeffBTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepbarRangeCoeffCTable = BuildRangeCoeffCTable(theRangepbarTable,
                             thepbarRangeCoeffCTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // invert the range table
       theInverseRangepbarTable = BuildInverseRangeTable(theRangepbarTable,
                              thepbarRangeCoeffATable,
                              thepbarRangeCoeffBTable,
                              thepbarRangeCoeffCTable,
                              theInverseRangepbarTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);
  
      }

    }
  }
  // make the energy loss and the range table available

  G4EnergyLossTables::Register(&aParticleType,  
    (Charge>0)?
      theDEDXpTable: theDEDXpbarTable,
    (Charge>0)?
      theRangepTable: theRangepbarTable,
    (Charge>0)?
      theInverseRangepTable: theInverseRangepbarTable,
    (Charge>0)?
      theLabTimepTable: theLabTimepbarTable,
    (Charge>0)?
      theProperTimepTable: theProperTimepbarTable,
    LowerBoundEloss, UpperBoundEloss,
    proton_mass_c2/aParticleType.GetPDGMass(),NbinEloss);

}
      

G4double G4hEnergyLoss::GetConstraints(const G4DynamicParticle *aParticle,
                                              G4Material *aMaterial)
{
  // returns the Step limit
  // dRoverRange is the max. allowed relative range loss in one step
  // it calculates dEdx and the range as well....

  G4double KineticEnergy,StepLimit;
  G4bool isOut ;

  Charge = aParticle->GetDefinition()->GetPDGCharge()/eplus ;

  KineticEnergy = aParticle->GetKineticEnergy();

  G4double massratio=proton_mass_c2/
           aParticle->GetDefinition()->GetPDGMass() ;

  G4double Tscaled= KineticEnergy*massratio ; 
  G4double ChargeSquare = Charge*Charge ;

     if(Charge>0.)
     {
       fRangeNow = G4EnergyLossTables::GetRange( theProton,
                                            Tscaled,aMaterial) ;
        fdEdx     = G4EnergyLossTables::GetDEDX( theProton,
                                            Tscaled,aMaterial) ;
     }
     else
     {
       fRangeNow = G4EnergyLossTables::GetRange( theAntiProton,
                                             Tscaled,aMaterial) ;
       fdEdx     = G4EnergyLossTables::GetDEDX( theAntiProton,
                                             Tscaled,aMaterial) ;
     }
     fdEdx     *= ChargeSquare ;
     fRangeNow /= (ChargeSquare*massratio) ;

  // compute the (random) Step limit ..............
  if(fRangeNow > finalRange)
  {
    StepLimit = (c1lim*fRangeNow+c2lim+c3lim/fRangeNow) ;

    //  randomise this value
    if(rndmStepFlag) StepLimit = 
                finalRange+(StepLimit-finalRange)*G4UniformRand() ;
    if(StepLimit > fRangeNow) StepLimit = fRangeNow ;
  }
  else StepLimit = fRangeNow ;


  return StepLimit ;
}

G4VParticleChange* G4hEnergyLoss::AlongStepDoIt( 
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
  ChargeSquare = Charge*Charge ;

  G4int index = aMaterial->GetIndex() ;
  E = aParticle->GetKineticEnergy() ;

  if(E < MinKineticEnergy) MeanLoss = E ;
  else
  {
    if(Step >= fRangeNow ) MeanLoss = E ;

    else if(( E > UpperBoundEloss)||( E <= LowerBoundEloss))
              MeanLoss = Step*fdEdx ; 
     
    else
    {
      if(Step>linLossLimit*fRangeNow)
      {
        G4double massratio=proton_mass_c2/
                 aParticle->GetDefinition()->GetPDGMass() ;

        G4double rscaled= fRangeNow*massratio*ChargeSquare ;
        G4double sscaled=   Step   *massratio*ChargeSquare ;

        if(Charge>0.)
        {
          MeanLoss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theProton,
                                         rscaled        ,aMaterial) -
                     G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theProton,
                                         rscaled-sscaled,aMaterial) ;
        }
        else
        {
          MeanLoss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theAntiProton,
                                         rscaled        ,aMaterial) -
                     G4EnergyLossTables::GetPreciseEnergyFromRange(
                                         theAntiProton,
                                         rscaled-sscaled,aMaterial) ;
        }
        MeanLoss /= (massratio*ChargeSquare) ;
      }
      else MeanLoss = Step*fdEdx ;
    }
  } 
  finalT = E - MeanLoss ;

  if(finalT < MinKineticEnergy) finalT = 0. ;

  //  now the loss with fluctuation
  if((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowerBoundEloss))
  {
    MeanLoss /= ChargeSquare ;
    finalT = E-GetLossWithFluct(aParticle,aMaterial,MeanLoss)*ChargeSquare ;
    if (finalT < 0.) finalT = 0.  ;
  }

  //  kill the particle if the kinetic energy <= 0  
  if (finalT <= 0. )
  {
    finalT = 0.;
    if(aParticle->GetDefinition()->GetParticleName() == "proton")
      aParticleChange.SetStatusChange(fStopAndKill);
    else  
      aParticleChange.SetStatusChange(fStopButAlive); 
  } 

  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(E-finalT) ;

  return &aParticleChange ;
}


