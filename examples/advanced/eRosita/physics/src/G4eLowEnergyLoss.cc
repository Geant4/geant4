//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//  
// -----------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4eLowEnergyLoss physics process -----------
//                by Laszlo Urban, 20 March 1997 
// **************************************************************
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
// 19-09-00  change of fluctuation sampling V.Ivanchenko
// 20/09/00  update fluctuations V.Ivanchenko
// 18/10/01  add fluorescence AlongStepDoIt V.Ivanchenko
// 18/10/01  Revision to improve code quality and consistency with design, MGP
// 19/10/01  update according to new design, V.Ivanchenko
// 24/10/01  MGP - Protection against negative energy loss in AlongStepDoIt
// 26/10/01  VI Clean up access to deexcitation
// 23/11/01  VI Move static member-functions from header to source
// 28/05/02  VI Remove flag fStopAndKill
// 03/06/02  MGP - Restore fStopAndKill
// 28/10/02  VI Optimal binning for dE/dx
// 21/01/03  VI cut per region
// 01/06/04  VI check if stopped particle has AtRest processes
//
// --------------------------------------------------------------

#include "G4eLowEnergyLoss.hh"
#include "G4SystemOfUnits.hh"
#include "G4EnergyLossMessenger.hh"
#include "G4Poisson.hh"
#include "G4ProductionCutsTable.hh"

//

// Initialisation of static data members
// -------------------------------------
// Contributing processes : ion.loss + soft brems->NbOfProcesses is initialized
// to 2 . YOU DO NOT HAVE TO CHANGE this variable for a 'normal' run.
//
// You have to change NbOfProcesses if you invent a new process contributing
// to the continuous energy loss.
// The NbOfProcesses data member can be changed using the (public static)
// functions Get/Set/Plus/MinusNbOfProcesses (see G4eLowEnergyLoss.hh)

G4int            G4eLowEnergyLoss::NbOfProcesses = 2;

G4int            G4eLowEnergyLoss::CounterOfElectronProcess = 0;
G4int            G4eLowEnergyLoss::CounterOfPositronProcess = 0;
G4PhysicsTable** G4eLowEnergyLoss::RecorderOfElectronProcess =
                                           new G4PhysicsTable*[10];
G4PhysicsTable** G4eLowEnergyLoss::RecorderOfPositronProcess =
                                           new G4PhysicsTable*[10];


G4PhysicsTable*  G4eLowEnergyLoss::theDEDXElectronTable         = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theDEDXPositronTable         = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theRangeElectronTable        = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theRangePositronTable        = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theInverseRangeElectronTable = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theInverseRangePositronTable = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theLabTimeElectronTable      = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theLabTimePositronTable      = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theProperTimeElectronTable   = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theProperTimePositronTable   = 0;

G4PhysicsTable*  G4eLowEnergyLoss::theeRangeCoeffATable         = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theeRangeCoeffBTable         = 0;
G4PhysicsTable*  G4eLowEnergyLoss::theeRangeCoeffCTable         = 0;
G4PhysicsTable*  G4eLowEnergyLoss::thepRangeCoeffATable         = 0;
G4PhysicsTable*  G4eLowEnergyLoss::thepRangeCoeffBTable         = 0;
G4PhysicsTable*  G4eLowEnergyLoss::thepRangeCoeffCTable         = 0;

G4double         G4eLowEnergyLoss::LowerBoundEloss = 10.*eV ;
G4double         G4eLowEnergyLoss::UpperBoundEloss = 100.*GeV ;
G4int            G4eLowEnergyLoss::NbinEloss = 360 ;
G4double         G4eLowEnergyLoss::RTable ;
G4double         G4eLowEnergyLoss::LOGRTable ;


G4EnergyLossMessenger* G4eLowEnergyLoss::eLossMessenger         = 0;

//    
 
// constructor and destructor
 
G4eLowEnergyLoss::G4eLowEnergyLoss(const G4String& processName)
   : G4RDVeLowEnergyLoss (processName),
     theLossTable(0),
     MinKineticEnergy(1.*eV),
     Charge(-1.),
     lastCharge(0.),
     theDEDXTable(0),
     CounterOfProcess(0),
     RecorderOfProcess(0),
     fdEdx(0),
     fRangeNow(0),
     linLossLimit(0.05),
     theFluo(false)
{
 
 //create (only once) EnergyLoss messenger 
 if(!eLossMessenger) eLossMessenger = new G4EnergyLossMessenger();
}

//

G4eLowEnergyLoss::~G4eLowEnergyLoss() 
{
     if (theLossTable) 
       {
         theLossTable->clearAndDestroy();
         delete theLossTable;
       }
}

void G4eLowEnergyLoss::SetNbOfProcesses(G4int nb) 
{
    NbOfProcesses=nb;
}

void G4eLowEnergyLoss::PlusNbOfProcesses()        
{
    NbOfProcesses++;
}

void G4eLowEnergyLoss::MinusNbOfProcesses() 
{
    NbOfProcesses--;
}                                      

G4int G4eLowEnergyLoss::GetNbOfProcesses() 
{
    return NbOfProcesses;
}
    
void G4eLowEnergyLoss::SetLowerBoundEloss(G4double val) 
{
    LowerBoundEloss=val;
}
    
void G4eLowEnergyLoss::SetUpperBoundEloss(G4double val) 
{
    UpperBoundEloss=val;
} 

void G4eLowEnergyLoss::SetNbinEloss(G4int nb)
{
    NbinEloss=nb;
}
 
G4double G4eLowEnergyLoss::GetLowerBoundEloss()
{
    return LowerBoundEloss;
} 
    
G4double G4eLowEnergyLoss::GetUpperBoundEloss() 
{
    return UpperBoundEloss;
} 

G4int G4eLowEnergyLoss::GetNbinEloss() 
{
    return NbinEloss;
} 
//     

void G4eLowEnergyLoss::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{
  ParticleMass = aParticleType.GetPDGMass(); 
  Charge = aParticleType.GetPDGCharge()/eplus;

  //  calculate data members LOGRTable,RTable first

  G4double lrate = std::log(UpperBoundEloss/LowerBoundEloss);
  LOGRTable=lrate/NbinEloss;
  RTable   =std::exp(LOGRTable);
  // Build energy loss table as a sum of the energy loss due to the
  // different processes.
  //

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

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
         theDEDXElectronTable = new G4PhysicsTable(numOfCouples);
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
        theDEDXPositronTable = new G4PhysicsTable(numOfCouples);
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

     for (size_t J=0; J<numOfCouples; J++)
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

     if(verboseLevel > 1) {
       G4cout << (*theDEDXElectronTable) << G4endl;
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

//

G4VParticleChange* G4eLowEnergyLoss::AlongStepDoIt( const G4Track& trackData,
                                                 const G4Step&  stepData)
{
 // compute the energy loss after a Step

  static const G4double faclow = 1.5 ;

  // get particle and material pointers from trackData
  const G4DynamicParticle* aParticle = trackData.GetDynamicParticle();
  G4double E      = aParticle->GetKineticEnergy();

  const G4MaterialCutsCouple* couple = trackData.GetMaterialCutsCouple();

  G4double Step = stepData.GetStepLength();

  aParticleChange.Initialize(trackData);
  //fParticleChange.Initialize(trackData);

  G4double MeanLoss, finalT;

  if (E < MinKineticEnergy)   finalT = 0.;

  else if ( E< faclow*LowerBoundEloss)
  {
    if (Step >= fRangeNow)  finalT = 0.;
   //  else finalT = E*(1.-Step/fRangeNow) ;
    else finalT = E*(1.-std::sqrt(Step/fRangeNow)) ;
  }

  else if (E>=UpperBoundEloss) finalT = E - Step*fdEdx;

  else if (Step >= fRangeNow)  finalT = 0.;

  else
  {
    if(Step/fRangeNow < linLossLimit) finalT = E-Step*fdEdx ;
    else
    {
      if (Charge<0.) finalT = G4EnergyLossTables::GetPreciseEnergyFromRange
                             (G4Electron::Electron(),fRangeNow-Step,couple);
      else           finalT = G4EnergyLossTables::GetPreciseEnergyFromRange
                             (G4Positron::Positron(),fRangeNow-Step,couple);
     }
  }

  if(finalT < MinKineticEnergy) finalT = 0. ;

  MeanLoss = E-finalT ;

  //now the loss with fluctuation
  if ((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowerBoundEloss))
  {
    finalT = E-GetLossWithFluct(aParticle,couple,MeanLoss,Step);
    if (finalT < 0.) finalT = 0.;
  }

  // kill the particle if the kinetic energy <= 0
  if (finalT <= 0. )
  {
    finalT = 0.;
    if(Charge > 0.0) aParticleChange.ProposeTrackStatus(fStopButAlive);
    else             aParticleChange.ProposeTrackStatus(fStopAndKill);
  }

  G4double edep = E - finalT;

  aParticleChange.ProposeEnergy(finalT);

  // Deexcitation of ionised atoms
  std::vector<G4DynamicParticle*>* deexcitationProducts = 0;
  if (theFluo) deexcitationProducts = DeexciteAtom(couple,E,edep);

  size_t nSecondaries = 0;
  if (deexcitationProducts != 0) nSecondaries = deexcitationProducts->size();
  aParticleChange.SetNumberOfSecondaries(nSecondaries);

  if (nSecondaries > 0) {

    const G4StepPoint* preStep = stepData.GetPreStepPoint();
    const G4StepPoint* postStep = stepData.GetPostStepPoint();
    G4ThreeVector r = preStep->GetPosition();
    G4ThreeVector deltaR = postStep->GetPosition();
    deltaR -= r;
    G4double t = preStep->GetGlobalTime();
    G4double deltaT = postStep->GetGlobalTime();
    deltaT -= t;
    G4double time, q;
    G4ThreeVector position;

    for (size_t i=0; i<nSecondaries; i++) {

      G4DynamicParticle* part = (*deexcitationProducts)[i];
      if (part != 0) {
        G4double eSecondary = part->GetKineticEnergy();
        edep -= eSecondary;
	if (edep > 0.)
	  {
	    q = G4UniformRand();
	    time = deltaT*q + t;
	    position  = deltaR*q;
	    position += r;
	    G4Track* newTrack = new G4Track(part, time, position);
	    aParticleChange.AddSecondary(newTrack);
	  }
	else
	  {
	    edep += eSecondary;
	    delete part;
	    part = 0;
	  }
      }
    }
  }
  delete deexcitationProducts;

  aParticleChange.ProposeLocalEnergyDeposit(edep);

  return &aParticleChange;
}

//

