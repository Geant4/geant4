// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIenergyLoss.cc,v 1.4 2000-02-10 09:06:29 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4PAIenergyLoss physics process -----------
//                by V. Grichine, 30 Nov 1997 
// **************************************************************
// It is the first implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of charged hadrons.
// **************************************************************
//
// corrected by V. Grichine on 24/11/97
// corrected by L. Urban    on 27/05/98   ( other corrections come soon!)
// 10/02/00  modifications , new e.m. structure, L.Urban
//
 

#include "G4PAIenergyLoss.hh"
#include "G4PAIonisation.hh"
#include "G4EnergyLossTables.hh"

////////////////////////////////////////////////////////////////////////////
//
// Initialisation of static members 

G4int            G4PAIenergyLoss::NbOfProcesses      = 1 ;
G4PhysicsTable** G4PAIenergyLoss::RecorderOfpProcess =
                                           new G4PhysicsTable*[10] ;
G4int       G4PAIenergyLoss::CounterOfpProcess = 0 ;

G4PhysicsTable* G4PAIenergyLoss::theDEDXpTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::theRangepTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::theInverseRangepTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::theLabTimepTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::theProperTimepTable = NULL ;

G4PhysicsTable* G4PAIenergyLoss::thepRangeCoeffATable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::thepRangeCoeffBTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::thepRangeCoeffCTable = NULL ;

G4PhysicsTable** G4PAIenergyLoss::RecorderOfpbarProcess =
                                           new G4PhysicsTable*[10] ;
G4int       G4PAIenergyLoss::CounterOfpbarProcess = 0 ;

G4PhysicsTable* G4PAIenergyLoss::theDEDXpbarTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::theRangepbarTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::theInverseRangepbarTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::theLabTimepbarTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::theProperTimepbarTable = NULL ;

G4PhysicsTable* G4PAIenergyLoss::thepbarRangeCoeffATable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::thepbarRangeCoeffBTable = NULL ;
G4PhysicsTable* G4PAIenergyLoss::thepbarRangeCoeffCTable = NULL ;

// G4PhysicsTable* G4PAIenergyLoss::fPAItransferBank = NULL ;

G4double G4PAIenergyLoss::CutInRange = 0;

G4double G4PAIenergyLoss::LowerBoundEloss= 1.00*keV ;
G4double G4PAIenergyLoss::UpperBoundEloss= 100.*TeV ;
G4int G4PAIenergyLoss::NbinEloss =100 ;
G4double G4PAIenergyLoss::RTable,G4PAIenergyLoss::LOGRTable;


// constructor and destructor
 
G4PAIenergyLoss::G4PAIenergyLoss(const G4String& processName)
   : G4VEnergyLoss (processName),
     dToverTini(0.20),   // max.relative range loss in one Step = 20%
     theElectron ( G4Electron::Electron() ),
     theProton ( G4Proton::Proton() ),
     theAntiProton ( G4AntiProton::AntiProton() )
{
     theLossTable = NULL ;
     theDEDXTable = NULL ;

//  calculate data members LOGRTable,RTable first
  G4double lrate ;
  lrate = log(UpperBoundEloss/LowerBoundEloss) ;
  LOGRTable=lrate/NbinEloss;
  RTable   =exp(LOGRTable);

}

G4PAIenergyLoss::~G4PAIenergyLoss() 
{
     if(theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable;
     }

}
 
 
/////////////////////////////////////////////////////////////////////////
//
//

void G4PAIenergyLoss::BuildDEDXTable(const G4ParticleDefinition& aParticleType)
{

  G4bool MakeTable = false ;
     
  G4double newCutInRange = aParticleType.GetLengthCuts();

// Create tables only if there is a new cut value !

  // create/fill proton or antiproton tables depending on the charge of the particle
  G4double Charge = aParticleType.GetPDGCharge();

  if (Charge>0.) 
  {
    theDEDXTable= theDEDXpTable;
  } 
  else 
  {
    theDEDXTable= theDEDXpbarTable;
  }

  if ((CutInRange != newCutInRange) || (theDEDXTable==NULL)) 
  {
    MakeTable = true ;
    CutInRange = newCutInRange ;
  }
  
  if( MakeTable )
  {

// Build energy loss table as a sum of the energy loss due to the
//              different processes.                                           
//
//  different processes.                                           

    const G4MaterialTable* theMaterialTable=
                                     G4Material::GetMaterialTable();

//  create table for the total energy loss

    G4int numOfMaterials = theMaterialTable->length();

  G4PhysicsTable** RecorderOfProcess;
  int CounterOfProcess;

 if( Charge >0.)    
 {
    RecorderOfProcess=RecorderOfpProcess;
    CounterOfProcess=CounterOfpProcess;

    if(CounterOfProcess == NbOfProcesses)
    {
  // create tables
      if(theDEDXpTable)
      { theDEDXpTable->clearAndDestroy();
        delete theDEDXpTable; }
      theDEDXpTable = new G4PhysicsTable(numOfMaterials);
      theDEDXTable = theDEDXpTable;
    }
  }
  else
  {
    RecorderOfProcess=RecorderOfpbarProcess;
    CounterOfProcess=CounterOfpbarProcess;

    if(CounterOfProcess == NbOfProcesses)
    {
  // create tables
      if(theDEDXpbarTable)
      { theDEDXpbarTable->clearAndDestroy();
        delete theDEDXpbarTable; }
      theDEDXpbarTable = new G4PhysicsTable(numOfMaterials);
      theDEDXTable = theDEDXpbarTable;
    }
  }


  if(CounterOfProcess == NbOfProcesses)
  {
 // fill the tables

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
//     here comes the sum of the different tables created by the  
//     processes (ionisation,etc...)              

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
 
  //  reset counter to zero 

  if( Charge >0.)    CounterOfpProcess=0 ;
  else  CounterOfpbarProcess=0 ;

      ParticleMass = aParticleType.GetPDGMass() ;

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


//////////////////////////////////////////////////////////////////////
   

























