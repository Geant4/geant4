// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4eLowEnergyLoss.cc,v 1.3 2000-04-11 10:08:28 lefebure Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//  
// -----------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4eLowEnergyLoss physics process -----------
//                by Laszlo Urban, 20 March 1997 
// **************************************************************
// It is the first implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of e+/e-.
// --------------------------------------------------------------
//
// 08-05-97: small changes by L.Urban
// 27-05-98: several bugs and inconsistencies are corrected,
//           new table (the inverse of the range table) added ,
//           AlongStepDoit uses now this new table. L.Urban
// 08-09-98: cleanup
// 24-09-98: rndmStepFlag false by default (no randomization of the step)
// 14-10-98: messenger file added.
// 16-10-98: public method SetStepFunction() 
// 20-01-99: important correction in AlongStepDoIt , L.Urban
// 10/02/00  modifications , new e.m. structure, L.Urban
// 11/04/00: Bug fix in dE/dx fluctuation simulation, Veronique Lefebure
// --------------------------------------------------------------
 
#include "G4eLowEnergyLoss.hh"
#include "G4EnergyLossMessenger.hh"
#include "G4Poisson.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Initialisation of static data members
// -------------------------------------
// Contributing processes : ion.loss + soft brems->NbOfProcesses is initialized
// to 2 . YOU DO NOT HAVE TO CHANGE this variable for a 'normal' run.
//
// You have to change NbOfProcesses if you invent a new process contributing
// to the continuous energy loss.
// The NbOfProcesses data member can be changed using the (public static)
// functions Get/Set/Plus/MinusNbOfProcesses (see G4eLowEnergyLoss.hh)

G4int            G4eLowEnergyLoss::NbOfProcesses = 1;

G4int            G4eLowEnergyLoss::CounterOfElectronProcess = 0;
G4int            G4eLowEnergyLoss::CounterOfPositronProcess = 0;
G4PhysicsTable** G4eLowEnergyLoss::RecorderOfElectronProcess =
                                           new G4PhysicsTable*[10];
G4PhysicsTable** G4eLowEnergyLoss::RecorderOfPositronProcess =
                                           new G4PhysicsTable*[10];
                                           

G4PhysicsTable*  G4eLowEnergyLoss::theDEDXElectronTable         = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theDEDXPositronTable         = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theRangeElectronTable        = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theRangePositronTable        = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theInverseRangeElectronTable = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theInverseRangePositronTable = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theLabTimeElectronTable      = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theLabTimePositronTable      = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theProperTimeElectronTable   = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theProperTimePositronTable   = NULL;

G4PhysicsTable*  G4eLowEnergyLoss::theeRangeCoeffATable         = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theeRangeCoeffBTable         = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::theeRangeCoeffCTable         = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::thepRangeCoeffATable         = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::thepRangeCoeffBTable         = NULL;
G4PhysicsTable*  G4eLowEnergyLoss::thepRangeCoeffCTable         = NULL;

G4double         G4eLowEnergyLoss::LowerBoundEloss = 250.*eV ;
G4double         G4eLowEnergyLoss::UpperBoundEloss = 100.*GeV ;
G4int            G4eLowEnergyLoss::NbinEloss = 100 ;
G4double         G4eLowEnergyLoss::RTable ;
G4double         G4eLowEnergyLoss::LOGRTable ;


G4EnergyLossMessenger* G4eLowEnergyLoss::eLossMessenger         = NULL;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// constructor and destructor
 
G4eLowEnergyLoss::G4eLowEnergyLoss(const G4String& processName)
   : G4VeLowEnergyLoss (processName),
     theLossTable(NULL),
     theDEDXTable(NULL),
     Charge(-1.),lastCharge(0.),
     MinKineticEnergy(1.*eV),
     //linLossLimit(0.02)
     linLossLimit(0.05)
{
 //create (only once) EnergyLoss messenger 
 if(!eLossMessenger) eLossMessenger = new G4EnergyLossMessenger();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eLowEnergyLoss::~G4eLowEnergyLoss() 
{
     if (theLossTable) 
       {
         theLossTable->clearAndDestroy();
         delete theLossTable;
       }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

void G4eLowEnergyLoss::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{
  ParticleMass = aParticleType.GetPDGMass(); 
  Charge = aParticleType.GetPDGCharge()/eplus;

  //  calculate data members LOGRTable,RTable first

  G4double lrate = log(UpperBoundEloss/LowerBoundEloss);
  LOGRTable=lrate/NbinEloss;
  RTable   =exp(LOGRTable);
  // Build energy loss table as a sum of the energy loss due to the
  // different processes.                                           
  //

  const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();
  
  // create table for the total energy loss

  if (&aParticleType==G4Electron::Electron())
    {
      RecorderOfProcess=RecorderOfElectronProcess;
      CounterOfProcess=CounterOfElectronProcess;
      if (CounterOfProcess == NbOfProcesses)
        {
         if (theDEDXElectronTable)
           { 
             theDEDXElectronTable->clearAndDestroy();
             delete theDEDXElectronTable; 
           }
         theDEDXElectronTable = new G4PhysicsTable(numOfMaterials);
         theDEDXTable = theDEDXElectronTable;
        }
    }
  if (&aParticleType==G4Positron::Positron())
    {
     RecorderOfProcess=RecorderOfPositronProcess;
     CounterOfProcess=CounterOfPositronProcess;
     if (CounterOfProcess == NbOfProcesses)
       {
        if (theDEDXPositronTable)
          { 
            theDEDXPositronTable->clearAndDestroy();
            delete theDEDXPositronTable; 
          }
        theDEDXPositronTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXPositronTable;
       }
    }

  if (CounterOfProcess == NbOfProcesses)
    {
     // fill the tables
     // loop for materials
     G4double LowEdgeEnergy , Value;
     G4bool isOutRange;
     G4PhysicsTable* pointer;

     for (G4int J=0; J<numOfMaterials; J++)
        {
         // create physics vector and fill it

         G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowerBoundEloss, UpperBoundEloss, NbinEloss);   

         // loop for the kinetic energy
   
         for (G4int i=0; i<NbinEloss; i++) 
            {
              LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;      
              //here comes the sum of the different tables created by the  
              //processes (ionisation,bremsstrahlung,etc...)              
              Value = 0.;    
              for (G4int process=0; process < NbOfProcesses; process++)
                 {
                   pointer= RecorderOfProcess[process];
                   Value += (*pointer)[J]->GetValue(LowEdgeEnergy,isOutRange);
                 }

              aVector->PutValue(i,Value) ; 
            }

         theDEDXTable->insert(aVector) ;

        }

 
     //reset counter to zero
     if (&aParticleType==G4Electron::Electron()) CounterOfElectronProcess=0;
     if (&aParticleType==G4Positron::Positron()) CounterOfPositronProcess=0;

     ParticleMass = aParticleType.GetPDGMass(); 

     if (&aParticleType==G4Electron::Electron())
     {
       // Build range table
       theRangeElectronTable = BuildRangeTable(theDEDXElectronTable,
                                               theRangeElectronTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build lab/proper time tables
       theLabTimeElectronTable = BuildLabTimeTable(theDEDXElectronTable,
                         theLabTimeElectronTable,
                         LowerBoundEloss,UpperBoundEloss,NbinEloss);
       theProperTimeElectronTable = BuildProperTimeTable(theDEDXElectronTable,
                            theProperTimeElectronTable,
                            LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build coeff tables for the energy loss calculation
       theeRangeCoeffATable = BuildRangeCoeffATable(theRangeElectronTable,
                             theeRangeCoeffATable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       theeRangeCoeffBTable = BuildRangeCoeffBTable(theRangeElectronTable,
                             theeRangeCoeffBTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       theeRangeCoeffCTable = BuildRangeCoeffCTable(theRangeElectronTable,
                             theeRangeCoeffCTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // invert the range table
       theInverseRangeElectronTable = BuildInverseRangeTable(theRangeElectronTable,
                              theeRangeCoeffATable,
                              theeRangeCoeffBTable,
                              theeRangeCoeffCTable,
                              theInverseRangeElectronTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);
     }
     if (&aParticleType==G4Positron::Positron())
     {
       // Build range table
       theRangePositronTable = BuildRangeTable(theDEDXPositronTable,
                                               theRangePositronTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);


       // Build lab/proper time tables
       theLabTimePositronTable = BuildLabTimeTable(theDEDXPositronTable,
                         theLabTimePositronTable,
                         LowerBoundEloss,UpperBoundEloss,NbinEloss);
       theProperTimePositronTable = BuildProperTimeTable(theDEDXPositronTable,
                            theProperTimePositronTable,
                            LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // Build coeff tables for the energy loss calculation
       thepRangeCoeffATable = BuildRangeCoeffATable(theRangePositronTable,
                             thepRangeCoeffATable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepRangeCoeffBTable = BuildRangeCoeffBTable(theRangePositronTable,
                             thepRangeCoeffBTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       thepRangeCoeffCTable = BuildRangeCoeffCTable(theRangePositronTable,
                             thepRangeCoeffCTable,
                             LowerBoundEloss,UpperBoundEloss,NbinEloss);

       // invert the range table
       theInverseRangePositronTable = BuildInverseRangeTable(theRangePositronTable,
                              thepRangeCoeffATable,
                              thepRangeCoeffBTable,
                              thepRangeCoeffCTable,
                              theInverseRangePositronTable,
                              LowerBoundEloss,UpperBoundEloss,NbinEloss);
     }

     // make the energy loss and the range table available
     G4EnergyLossTables::Register(&aParticleType,  
       (&aParticleType==G4Electron::Electron())?
       theDEDXElectronTable: theDEDXPositronTable,
       (&aParticleType==G4Electron::Electron())?
       theRangeElectronTable: theRangePositronTable,
       (&aParticleType==G4Electron::Electron())?
       theInverseRangeElectronTable: theInverseRangePositronTable,
       (&aParticleType==G4Electron::Electron())?
       theLabTimeElectronTable: theLabTimePositronTable,
       (&aParticleType==G4Electron::Electron())?
       theProperTimeElectronTable: theProperTimePositronTable,
       LowerBoundEloss, UpperBoundEloss, 1.,NbinEloss);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
      
G4VParticleChange* G4eLowEnergyLoss::AlongStepDoIt( const G4Track& trackData,
                                                 const G4Step&  stepData)
{                              
 // compute the energy loss after a Step

  static const G4double faclow = 1.5 ;

  // get particle and material pointers from trackData 
  const G4DynamicParticle* aParticle = trackData.GetDynamicParticle();
  G4double E      = aParticle->GetKineticEnergy() ;
  
  G4Material* aMaterial = trackData.GetMaterial();
  G4int index = aMaterial->GetIndex();
  
  G4double Step = stepData.GetStepLength();

  fParticleChange.Initialize(trackData);  
  
  G4double MeanLoss, finalT; 

  if (E < MinKineticEnergy)   finalT = 0.; 
  
  else if ( E< faclow*LowerBoundEloss)  
  {
    if (Step >= fRangeNow)  finalT = 0.; 
   //  else finalT = E*(1.-Step/fRangeNow) ;
    else finalT = E*(1.-sqrt(Step/fRangeNow)) ;
  }
    
  else if (E>=UpperBoundEloss) finalT = E - Step*fdEdx;
     
  else if (Step >= fRangeNow)  finalT = 0.; 
  
  else
  {
    if(Step/fRangeNow < linLossLimit) finalT = E-Step*fdEdx ;
    else
    {
      if (Charge<0.) finalT = G4EnergyLossTables::GetPreciseEnergyFromRange
                             (G4Electron::Electron(),fRangeNow-Step,aMaterial);
      else           finalT = G4EnergyLossTables::GetPreciseEnergyFromRange
                             (G4Positron::Positron(),fRangeNow-Step,aMaterial);
     }
  }

  if(finalT < MinKineticEnergy) finalT = 0. ;

  MeanLoss = E-finalT ;  
  
  //now the loss with fluctuation
  if ((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowerBoundEloss))
  {
    finalT = E-GetLossWithFluct(aParticle,aMaterial,MeanLoss);
    if (finalT < 0.) finalT = 0.;
  }


  // kill the particle if the kinetic energy <= 0  
  if (finalT <= 0. )
  {
    finalT = 0.;
    if (Charge < 0.) fParticleChange.SetStatusChange(fStopAndKill);
    else             fParticleChange.SetStatusChange(fStopButAlive); 
  } 

  fParticleChange.SetEnergyChange(finalT);
  fParticleChange.SetLocalEnergyDeposit(E-finalT);

  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


