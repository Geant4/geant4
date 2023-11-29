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
//      ---------- G4hEnergyLoss physics process -----------
//                by Laszlo Urban, 30 May 1997
//
// **************************************************************
// It is the first implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of charged hadrons.
// **************************************************************
//
// 7/10/98    bug fixes + some cleanup , L.Urban
// 22/10/98   cleanup , L.Urban
// 07/12/98   works for ions as well+ bug corrected, L.Urban
// 02/02/99   several bugs fixed, L.Urban
// 31/03/00   rename to lowenergy as G4hRDEnergyLoss.cc V.Ivanchenko
// 05/11/00   new method to calculate particle ranges
// 10/05/01   V.Ivanchenko Clean up againist Linux compilation with -Wall
// 23/11/01   V.Ivanchenko Move static member-functions from header to source
// 28/10/02   V.Ivanchenko Optimal binning for dE/dx
// 21/01/03   V.Ivanchenko Cut per region
// 23/01/03   V.Ivanchenko Fix in table build

// 31 Jul 2008 MGP     Short term supply of energy loss of hadrons through clone of 
//                     former G4hLowEnergyLoss (with some initial cleaning)
//                     To be replaced by reworked class to deal with condensed/discrete 
//                     issues properly

// 25 Nov 2010 MGP     Added some protections for FPE (mostly division by 0)
//                     The whole energy loss domain would profit from a design iteration

// 20 Jun 2011 MGP     Corrected some compilation warnings. The whole class will be heavily refactored anyway.

// --------------------------------------------------------------

#include "G4hRDEnergyLoss.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EnergyLossTables.hh"
#include "G4Poisson.hh"
#include "G4ProductionCutsTable.hh"


// Initialisation of static members ******************************************
// contributing processes : ion.loss ->NumberOfProcesses is initialized
//   to 1 . YOU DO NOT HAVE TO CHANGE this variable for a 'normal' run.
// You have to change NumberOfProcesses
// if you invent a new process contributing to the cont. energy loss,
//   NumberOfProcesses should be 2 in this case,
//  or for debugging purposes.
//  The NumberOfProcesses data member can be changed using the (public static)
//  functions Get/Set/Plus/MinusNumberOfProcesses (see G4hRDEnergyLoss.hh)


// The following vectors have a fixed dimension this means that if they are
// filled in with more than 100 elements will corrupt the memory.
G4ThreadLocal G4int G4hRDEnergyLoss::NumberOfProcesses = 1 ;

G4ThreadLocal G4int            G4hRDEnergyLoss::CounterOfProcess = 0 ;
G4ThreadLocal G4PhysicsTable** G4hRDEnergyLoss::RecorderOfProcess  = 0 ;

G4ThreadLocal G4int            G4hRDEnergyLoss::CounterOfpProcess = 0 ;
G4ThreadLocal G4PhysicsTable** G4hRDEnergyLoss::RecorderOfpProcess  = 0 ;

G4ThreadLocal G4int            G4hRDEnergyLoss::CounterOfpbarProcess = 0 ;
G4ThreadLocal G4PhysicsTable** G4hRDEnergyLoss::RecorderOfpbarProcess  = 0 ;

G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theDEDXpTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theDEDXpbarTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theRangepTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theRangepbarTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theInverseRangepTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theInverseRangepbarTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theLabTimepTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theLabTimepbarTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theProperTimepTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theProperTimepbarTable = 0 ;

G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::thepRangeCoeffATable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::thepRangeCoeffBTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::thepRangeCoeffCTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::thepbarRangeCoeffATable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::thepbarRangeCoeffBTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::thepbarRangeCoeffCTable = 0 ;

G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theDEDXTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theRangeTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theInverseRangeTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theLabTimeTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theProperTimeTable = 0 ;

G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theRangeCoeffATable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theRangeCoeffBTable = 0 ;
G4ThreadLocal G4PhysicsTable* G4hRDEnergyLoss::theRangeCoeffCTable = 0 ;

//const G4Proton* G4hRDEnergyLoss::theProton=G4Proton::Proton() ;
//const G4AntiProton* G4hRDEnergyLoss::theAntiProton=G4AntiProton::AntiProton() ;

G4ThreadLocal G4double G4hRDEnergyLoss::ParticleMass;
G4ThreadLocal G4double G4hRDEnergyLoss::ptableElectronCutInRange = 0.0 ;
G4ThreadLocal G4double G4hRDEnergyLoss::pbartableElectronCutInRange = 0.0 ;

G4ThreadLocal G4double G4hRDEnergyLoss::Mass,
                       G4hRDEnergyLoss::taulow, 
                       G4hRDEnergyLoss::tauhigh, 
                       G4hRDEnergyLoss::ltaulow, 
                       G4hRDEnergyLoss::ltauhigh; 

G4ThreadLocal G4double G4hRDEnergyLoss::dRoverRange = 0.20 ;
G4ThreadLocal G4double G4hRDEnergyLoss::finalRange = 0.2; // 200.*micrometer

G4ThreadLocal G4double G4hRDEnergyLoss::c1lim = 0.20 ;
G4ThreadLocal G4double G4hRDEnergyLoss::c2lim = 0.32 ;
 // 2.*(1.-(0.20))*(200.*micrometer)
G4ThreadLocal G4double G4hRDEnergyLoss::c3lim = -0.032 ;
 // -(1.-(0.20))*(200.*micrometer)*(200.*micrometer)

G4ThreadLocal G4double G4hRDEnergyLoss::Charge ;   

G4ThreadLocal G4bool   G4hRDEnergyLoss::rndmStepFlag   = false ;
G4ThreadLocal G4bool   G4hRDEnergyLoss::EnlossFlucFlag = true ;

G4ThreadLocal G4double G4hRDEnergyLoss::LowestKineticEnergy = 1e-05 ; // 10.*eV;
G4ThreadLocal G4double G4hRDEnergyLoss::HighestKineticEnergy= 1.e5; // 100.*GeV;
G4ThreadLocal G4int    G4hRDEnergyLoss::TotBin = 360;
G4ThreadLocal G4double G4hRDEnergyLoss::RTable, G4hRDEnergyLoss::LOGRTable;

//--------------------------------------------------------------------------------

// constructor and destructor
 
G4hRDEnergyLoss::G4hRDEnergyLoss(const G4String& processName)
  : G4VContinuousDiscreteProcess (processName),
    MaxExcitationNumber (1.e6),
    probLimFluct (0.01),
    nmaxDirectFluct (100),
    nmaxCont1(4),
    nmaxCont2(16),
    theLossTable(0),
    linLossLimit(0.05),
    MinKineticEnergy(0.0) 
{
  if (!RecorderOfpbarProcess) RecorderOfpbarProcess = new G4PhysicsTable*[100];
  if (!RecorderOfpProcess) RecorderOfpProcess = new G4PhysicsTable*[100];
  if (!RecorderOfProcess) RecorderOfProcess = new G4PhysicsTable*[100];
}

//--------------------------------------------------------------------------------

G4hRDEnergyLoss::~G4hRDEnergyLoss() 
{
  if(theLossTable)
  {
    theLossTable->clearAndDestroy();
    delete theLossTable;
  }
}

//--------------------------------------------------------------------------------

G4int G4hRDEnergyLoss::GetNumberOfProcesses()    
{
  return NumberOfProcesses; 
} 

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::SetNumberOfProcesses(G4int number)
{
  NumberOfProcesses=number; 
} 

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::PlusNumberOfProcesses()
{
  NumberOfProcesses++; 
} 

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::MinusNumberOfProcesses()
{
  NumberOfProcesses--; 
} 

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::SetdRoverRange(G4double value) 
{
  dRoverRange = value;
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::SetRndmStep (G4bool   value) 
{
  rndmStepFlag = value;
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::SetEnlossFluc (G4bool value) 
{
  EnlossFlucFlag = value;
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::SetStepFunction (G4double c1, G4double c2)
{
  dRoverRange = c1; 
  finalRange = c2;
  c1lim=dRoverRange;
  c2lim=2.*(1-dRoverRange)*finalRange;
  c3lim=-(1.-dRoverRange)*finalRange*finalRange;
}

//--------------------------------------------------------------------------------
 
void G4hRDEnergyLoss::BuildDEDXTable(const G4ParticleDefinition& aParticleType)
{
  //  calculate data members TotBin,LOGRTable,RTable first

  if (!RecorderOfpbarProcess) RecorderOfpbarProcess = new G4PhysicsTable*[100];
  if (!RecorderOfpProcess) RecorderOfpProcess = new G4PhysicsTable*[100];
  if (!RecorderOfProcess) RecorderOfProcess = new G4PhysicsTable*[100];

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();

  // create/fill proton or antiproton tables depending on the charge
  Charge = aParticleType.GetPDGCharge()/eplus;
  ParticleMass = aParticleType.GetPDGMass() ;

  if (Charge>0.) {theDEDXTable= theDEDXpTable;}
  else           {theDEDXTable= theDEDXpbarTable;}

  if( ((Charge>0.) && (theDEDXTable==0)) ||
      ((Charge<0.) && (theDEDXTable==0)) 
      )
    {

      // Build energy loss table as a sum of the energy loss due to the
      //              different processes.
      if( Charge >0.)
	{
	  RecorderOfProcess=RecorderOfpProcess;
	  CounterOfProcess=CounterOfpProcess;

	  if(CounterOfProcess == NumberOfProcesses)
	    {
	      if(theDEDXpTable)
		{ theDEDXpTable->clearAndDestroy();
		  delete theDEDXpTable; }
	      theDEDXpTable = new G4PhysicsTable(numOfCouples);
	      theDEDXTable = theDEDXpTable;
	    }
	}
      else
	{
	  RecorderOfProcess=RecorderOfpbarProcess;
	  CounterOfProcess=CounterOfpbarProcess;

	  if(CounterOfProcess == NumberOfProcesses)
	    {
	      if(theDEDXpbarTable)
		{ theDEDXpbarTable->clearAndDestroy();
		  delete theDEDXpbarTable; }
	      theDEDXpbarTable = new G4PhysicsTable(numOfCouples);
	      theDEDXTable = theDEDXpbarTable;
	    }
	}

      if(CounterOfProcess == NumberOfProcesses)
	{
	  //  loop for materials
	  G4double LowEdgeEnergy , Value ;
	  G4bool isOutRange ;
	  G4PhysicsTable* pointer ;

	  for (std::size_t J=0; J<numOfCouples; ++J)
	    {
	      // create physics vector and fill it
	      G4PhysicsLogVector* aVector =
                new G4PhysicsLogVector(LowestKineticEnergy,
                                       HighestKineticEnergy, TotBin);

	      // loop for the kinetic energy
	      for (G4int i=0; i<TotBin; i++)
		{
		  LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
		  Value = 0. ;

		  // loop for the contributing processes
		  for (G4int process=0; process < NumberOfProcesses; process++)
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

	  // Build range table
	  BuildRangeTable( aParticleType);

	  // Build lab/proper time tables
	  BuildTimeTables( aParticleType) ;

	  // Build coeff tables for the energy loss calculation
	  BuildRangeCoeffATable( aParticleType);
	  BuildRangeCoeffBTable( aParticleType);
	  BuildRangeCoeffCTable( aParticleType);

	  // invert the range table

	  BuildInverseRangeTable(aParticleType);
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
			       LowestKineticEnergy, HighestKineticEnergy,
			       proton_mass_c2/aParticleType.GetPDGMass(),
			       TotBin);

}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::BuildRangeTable(const G4ParticleDefinition& aParticleType)
{
  // Build range table from the energy loss table

  Mass = aParticleType.GetPDGMass();

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  if( Charge >0.)
    {
      if(theRangepTable)
	{ theRangepTable->clearAndDestroy();
	  delete theRangepTable; }
      theRangepTable = new G4PhysicsTable(numOfCouples);
      theRangeTable = theRangepTable ;
    }
  else
    {
      if(theRangepbarTable)
	{ theRangepbarTable->clearAndDestroy();
	  delete theRangepbarTable; }
      theRangepbarTable = new G4PhysicsTable(numOfCouples);
      theRangeTable = theRangepbarTable ;
    }

  // loop for materials

  for (G4int J=0;  J<numOfCouples; ++J)
    {
      G4PhysicsLogVector* aVector;
      aVector = new G4PhysicsLogVector(LowestKineticEnergy,
				       HighestKineticEnergy,TotBin);

      BuildRangeVector(J, aVector);
      theRangeTable->insert(aVector);
    }
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::BuildTimeTables(const G4ParticleDefinition& aParticleType)
{
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  if(&aParticleType == G4Proton::Proton())
    {
      if(theLabTimepTable)
	{ theLabTimepTable->clearAndDestroy();
	  delete theLabTimepTable; }
      theLabTimepTable = new G4PhysicsTable(numOfCouples);
      theLabTimeTable = theLabTimepTable ;

      if(theProperTimepTable)
	{ theProperTimepTable->clearAndDestroy();
	  delete theProperTimepTable; }
      theProperTimepTable = new G4PhysicsTable(numOfCouples);
      theProperTimeTable = theProperTimepTable ;
    }

  if(&aParticleType == G4AntiProton::AntiProton())
    {
      if(theLabTimepbarTable)
	{ theLabTimepbarTable->clearAndDestroy();
	  delete theLabTimepbarTable; }
      theLabTimepbarTable = new G4PhysicsTable(numOfCouples);
      theLabTimeTable = theLabTimepbarTable ;

      if(theProperTimepbarTable)
	{ theProperTimepbarTable->clearAndDestroy();
	  delete theProperTimepbarTable; }
      theProperTimepbarTable = new G4PhysicsTable(numOfCouples);
      theProperTimeTable = theProperTimepbarTable ;
    }

  for (G4int J=0;  J<numOfCouples; ++J)
    {
      G4PhysicsLogVector* aVector;
      G4PhysicsLogVector* bVector;

      aVector = new G4PhysicsLogVector(LowestKineticEnergy,
				       HighestKineticEnergy,TotBin);

      BuildLabTimeVector(J, aVector);
      theLabTimeTable->insert(aVector);

      bVector = new G4PhysicsLogVector(LowestKineticEnergy,
				       HighestKineticEnergy,TotBin);

      BuildProperTimeVector(J, bVector);
      theProperTimeTable->insert(bVector);
    }
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::BuildRangeVector(G4int materialIndex,
				       G4PhysicsLogVector* rangeVector)
{
  //  create range vector for a material

  G4bool isOut;
  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];
  G4double energy1 = rangeVector->GetLowEdgeEnergy(0);
  G4double dedx    = physicsVector->GetValue(energy1,isOut);
  G4double range   = 0.5*energy1/dedx;
  rangeVector->PutValue(0,range);
  G4int n = 100;
  G4double del = 1.0/(G4double)n ;

  for (G4int j=1; j<TotBin; j++) {

    G4double energy2 = rangeVector->GetLowEdgeEnergy(j);
    G4double de = (energy2 - energy1) * del ;
    G4double dedx1 = dedx ;

    for (G4int i=1; i<n; i++) {
      G4double energy = energy1 + i*de ;
      G4double dedx2  = physicsVector->GetValue(energy,isOut);
      range  += 0.5*de*(1.0/dedx1 + 1.0/dedx2);
      dedx1   = dedx2;
    }
    rangeVector->PutValue(j,range);
    dedx = dedx1 ;
    energy1 = energy2 ;
  }
}    

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::BuildLabTimeVector(G4int materialIndex,
					 G4PhysicsLogVector* timeVector)
{
  //  create lab time vector for a material

  G4int nbin=100;
  G4bool isOut;
  G4double tlim=5.*keV,parlowen=0.4,ppar=0.5-parlowen ;
  //G4double losslim,clim,taulim,timelim,ltaulim,ltaumax,
  G4double losslim,clim,taulim,timelim,
    LowEdgeEnergy,tau,Value ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];

  // low energy part first...
  losslim = physicsVector->GetValue(tlim,isOut);
  taulim=tlim/ParticleMass ;
  clim=std::sqrt(ParticleMass*tlim/2.)/(c_light*losslim*ppar) ;
  //ltaulim = std::log(taulim);
  //ltaumax = std::log(HighestKineticEnergy/ParticleMass) ;

  G4int i=-1;
  G4double oldValue = 0. ;
  G4double tauold ;
  do  
    {
      i += 1 ;
      LowEdgeEnergy = timeVector->GetLowEdgeEnergy(i);
      tau = LowEdgeEnergy/ParticleMass ;
      if ( tau <= taulim )
	{
	  Value = clim*std::exp(ppar*std::log(tau/taulim)) ;
	}
      else
	{
	  timelim=clim ;
	  ltaulow = std::log(taulim);
	  ltauhigh = std::log(tau);
	  Value = timelim+LabTimeIntLog(physicsVector,nbin);
	}
      timeVector->PutValue(i,Value);
      oldValue = Value ;
      tauold = tau ;
    } while (tau<=taulim) ;

  i += 1 ;
  for (G4int j=i; j<TotBin; j++)
    {
      LowEdgeEnergy = timeVector->GetLowEdgeEnergy(j);
      tau = LowEdgeEnergy/ParticleMass ;
      ltaulow = std::log(tauold);
      ltauhigh = std::log(tau);
      Value = oldValue+LabTimeIntLog(physicsVector,nbin);
      timeVector->PutValue(j,Value);
      oldValue = Value ;
      tauold = tau ;
    }
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::BuildProperTimeVector(G4int materialIndex,
					    G4PhysicsLogVector* timeVector)
{
  //  create proper time vector for a material

  G4int nbin=100;
  G4bool isOut;
  G4double tlim=5.*keV,parlowen=0.4,ppar=0.5-parlowen ;
  //G4double losslim,clim,taulim,timelim,ltaulim,ltaumax,
  G4double losslim,clim,taulim,timelim,
    LowEdgeEnergy,tau,Value ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];
 
  // low energy part first...
  losslim = physicsVector->GetValue(tlim,isOut);
  taulim=tlim/ParticleMass ;
  clim=std::sqrt(ParticleMass*tlim/2.)/(c_light*losslim*ppar) ;
  //ltaulim = std::log(taulim);
  //ltaumax = std::log(HighestKineticEnergy/ParticleMass) ;

  G4int i=-1;
  G4double oldValue = 0. ;
  G4double tauold ;
  do
    {
      i += 1 ;
      LowEdgeEnergy = timeVector->GetLowEdgeEnergy(i);
      tau = LowEdgeEnergy/ParticleMass ;
      if ( tau <= taulim )
	{
	  Value = clim*std::exp(ppar*std::log(tau/taulim)) ;
	}
      else
	{
	  timelim=clim ;
	  ltaulow = std::log(taulim);
	  ltauhigh = std::log(tau);
	  Value = timelim+ProperTimeIntLog(physicsVector,nbin);
	}
      timeVector->PutValue(i,Value);
      oldValue = Value ;
      tauold = tau ;
    } while (tau<=taulim) ;

  i += 1 ;
  for (G4int j=i; j<TotBin; j++)
    {
      LowEdgeEnergy = timeVector->GetLowEdgeEnergy(j);
      tau = LowEdgeEnergy/ParticleMass ;
      ltaulow = std::log(tauold);
      ltauhigh = std::log(tau);
      Value = oldValue+ProperTimeIntLog(physicsVector,nbin);
      timeVector->PutValue(j,Value);
      oldValue = Value ;
      tauold = tau ;
    }
}

//--------------------------------------------------------------------------------

G4double G4hRDEnergyLoss::RangeIntLin(G4PhysicsVector* physicsVector,
				      G4int nbin)
{
  //  num. integration, linear binning

  G4double dtau,Value,taui,ti,lossi,ci;
  G4bool isOut;
  dtau = (tauhigh-taulow)/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
    {
      taui = taulow + dtau*i ;
      ti = Mass*taui;
      lossi = physicsVector->GetValue(ti,isOut);
      if(i==0)
	ci=0.5;
      else
	{
	  if(i<nbin)
	    ci=1.;
	  else
	    ci=0.5;
	}
      Value += ci/lossi;
    }
  Value *= Mass*dtau;
  return Value;
}

//--------------------------------------------------------------------------------

G4double G4hRDEnergyLoss::RangeIntLog(G4PhysicsVector* physicsVector,
				      G4int nbin)
{
  //  num. integration, logarithmic binning

  G4double ltt,dltau,Value,ui,taui,ti,lossi,ci;
  G4bool isOut;
  ltt = ltauhigh-ltaulow;
  dltau = ltt/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
    {
      ui = ltaulow+dltau*i;
      taui = std::exp(ui);
      ti = Mass*taui;
      lossi = physicsVector->GetValue(ti,isOut);
      if(i==0)
	ci=0.5;
      else
	{
	  if(i<nbin)
	    ci=1.;
	  else
	    ci=0.5;
	}
      Value += ci*taui/lossi;
    }
  Value *= Mass*dltau;
  return Value;
}

//--------------------------------------------------------------------------------

G4double G4hRDEnergyLoss::LabTimeIntLog(G4PhysicsVector* physicsVector,
					G4int nbin)
{
  //  num. integration, logarithmic binning

  G4double ltt,dltau,Value,ui,taui,ti,lossi,ci;
  G4bool isOut;
  ltt = ltauhigh-ltaulow;
  dltau = ltt/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
    {
      ui = ltaulow+dltau*i;
      taui = std::exp(ui);
      ti = ParticleMass*taui;
      lossi = physicsVector->GetValue(ti,isOut);
      if(i==0)
	ci=0.5;
      else
	{
	  if(i<nbin)
	    ci=1.;
	  else
	    ci=0.5;
	}
      Value += ci*taui*(ti+ParticleMass)/(std::sqrt(ti*(ti+2.*ParticleMass))*lossi);
    }
  Value *= ParticleMass*dltau/c_light;
  return Value;
}

//--------------------------------------------------------------------------------

G4double G4hRDEnergyLoss::ProperTimeIntLog(G4PhysicsVector* physicsVector,
					   G4int nbin)
{
  //  num. integration, logarithmic binning

  G4double ltt,dltau,Value,ui,taui,ti,lossi,ci;
  G4bool isOut;
  ltt = ltauhigh-ltaulow;
  dltau = ltt/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
    {
      ui = ltaulow+dltau*i;
      taui = std::exp(ui);
      ti = ParticleMass*taui;
      lossi = physicsVector->GetValue(ti,isOut);
      if(i==0)
	ci=0.5;
      else
	{
	  if(i<nbin)
	    ci=1.;
	  else
	    ci=0.5;
	}
      Value += ci*taui*ParticleMass/(std::sqrt(ti*(ti+2.*ParticleMass))*lossi);
    }
  Value *= ParticleMass*dltau/c_light;
  return Value;
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::BuildRangeCoeffATable( const G4ParticleDefinition& )
{
  // Build tables of coefficients for the energy loss calculation
  //  create table for coefficients "A"

  G4int numOfCouples = (G4int)G4ProductionCutsTable::GetProductionCutsTable()->GetTableSize();

  if(Charge>0.)
    {
      if(thepRangeCoeffATable)
	{ thepRangeCoeffATable->clearAndDestroy();
	  delete thepRangeCoeffATable; }
      thepRangeCoeffATable = new G4PhysicsTable(numOfCouples);
      theRangeCoeffATable = thepRangeCoeffATable ;
      theRangeTable = theRangepTable ;
    }
  else
    {
      if(thepbarRangeCoeffATable)
	{ thepbarRangeCoeffATable->clearAndDestroy();
	  delete thepbarRangeCoeffATable; }
      thepbarRangeCoeffATable = new G4PhysicsTable(numOfCouples);
      theRangeCoeffATable = thepbarRangeCoeffATable ;
      theRangeTable = theRangepbarTable ;
    }
 
  G4double R2 = RTable*RTable ;
  G4double R1 = RTable+1.;
  G4double w = R1*(RTable-1.)*(RTable-1.);
  G4double w1 = RTable/w , w2 = -RTable*R1/w , w3 = R2/w ;
  G4double Ti , Tim , Tip , Ri , Rim , Rip , Value ;
  G4bool isOut;

  //  loop for materials
  for (G4int J=0; J<numOfCouples; J++)
    {
      G4int binmax=TotBin ;
      G4PhysicsLinearVector* aVector = 
	new G4PhysicsLinearVector(0.,binmax, TotBin);
      Ti = LowestKineticEnergy ;
      if (Ti < DBL_MIN) Ti = 1.e-8;
      G4PhysicsVector* rangeVector= (*theRangeTable)[J];

      for ( G4int i=0; i<TotBin; i++)
	{
	  Ri = rangeVector->GetValue(Ti,isOut) ;
	  if (Ti < DBL_MIN) Ti = 1.e-8;
	  if ( i==0 )
	    Rim = 0. ;
	  else
	    {
	      // ---- MGP ---- Modified to avoid a floating point exception
	      // The correction is just a temporary patch, the whole class should be redesigned 
	      // Original: Tim = Ti/RTable  results in 0./0. 
	      if (RTable != 0.)
		{
		  Tim = Ti/RTable ;
		}
	      else
		{
		  Tim = 0.;
		}
	      Rim = rangeVector->GetValue(Tim,isOut);
	    }
	  if ( i==(TotBin-1))
	    Rip = Ri ;
	  else
	    {
	      Tip = Ti*RTable ;
	      Rip = rangeVector->GetValue(Tip,isOut);
	    }
	  Value = (w1*Rip + w2*Ri + w3*Rim)/(Ti*Ti) ; 

	  aVector->PutValue(i,Value);
	  Ti = RTable*Ti ;
	}
  
      theRangeCoeffATable->insert(aVector);
    } 
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::BuildRangeCoeffBTable( const G4ParticleDefinition& )
{
  // Build tables of coefficients for the energy loss calculation
  //  create table for coefficients "B"

  G4int numOfCouples =
        (G4int)G4ProductionCutsTable::GetProductionCutsTable()->GetTableSize();

  if(Charge>0.)
    {
      if(thepRangeCoeffBTable)
	{ thepRangeCoeffBTable->clearAndDestroy();
	  delete thepRangeCoeffBTable; }
      thepRangeCoeffBTable = new G4PhysicsTable(numOfCouples);
      theRangeCoeffBTable = thepRangeCoeffBTable ;
      theRangeTable = theRangepTable ;
    }
  else
    {
      if(thepbarRangeCoeffBTable)
	{ thepbarRangeCoeffBTable->clearAndDestroy();
	  delete thepbarRangeCoeffBTable; }
      thepbarRangeCoeffBTable = new G4PhysicsTable(numOfCouples);
      theRangeCoeffBTable = thepbarRangeCoeffBTable ;
      theRangeTable = theRangepbarTable ;
    }

  G4double R2 = RTable*RTable ;
  G4double R1 = RTable+1.;
  G4double w = R1*(RTable-1.)*(RTable-1.);
  if (w < DBL_MIN) w = DBL_MIN;
  G4double w1 = -R1/w , w2 = R1*(R2+1.)/w , w3 = -R2*R1/w ;
  G4double Ti , Tim , Tip , Ri , Rim , Rip , Value ;
  G4bool isOut;

  //  loop for materials
  for (G4int J=0; J<numOfCouples; J++)
    {
      G4int binmax=TotBin ;
      G4PhysicsLinearVector* aVector =
	new G4PhysicsLinearVector(0.,binmax, TotBin);
      Ti = LowestKineticEnergy ;
      if (Ti < DBL_MIN) Ti = 1.e-8;
      G4PhysicsVector* rangeVector= (*theRangeTable)[J];
   
      for ( G4int i=0; i<TotBin; i++)
	{
	  Ri = rangeVector->GetValue(Ti,isOut) ;
	  if (Ti < DBL_MIN) Ti = 1.e-8;
	  if ( i==0 )
	    Rim = 0. ;
	  else
	    {
	      if (RTable < DBL_MIN) RTable = DBL_MIN;
	      Tim = Ti/RTable ;
	      Rim = rangeVector->GetValue(Tim,isOut);
	    }
	  if ( i==(TotBin-1))
	    Rip = Ri ;
	  else
	    {
	      Tip = Ti*RTable ;
	      Rip = rangeVector->GetValue(Tip,isOut);
	    }
	  if (Ti < DBL_MIN) Ti = DBL_MIN;
	  Value = (w1*Rip + w2*Ri + w3*Rim)/Ti;

	  aVector->PutValue(i,Value);
	  Ti = RTable*Ti ;
	}
      theRangeCoeffBTable->insert(aVector);
    } 
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::BuildRangeCoeffCTable( const G4ParticleDefinition& )
{
  // Build tables of coefficients for the energy loss calculation
  //  create table for coefficients "C"

  G4int numOfCouples =
        (G4int)G4ProductionCutsTable::GetProductionCutsTable()->GetTableSize();

  if(Charge>0.)
    {
      if(thepRangeCoeffCTable)
	{ thepRangeCoeffCTable->clearAndDestroy();
	  delete thepRangeCoeffCTable; }
      thepRangeCoeffCTable = new G4PhysicsTable(numOfCouples);
      theRangeCoeffCTable = thepRangeCoeffCTable ;
      theRangeTable = theRangepTable ;
    }
  else
    {
      if(thepbarRangeCoeffCTable)
	{ thepbarRangeCoeffCTable->clearAndDestroy();
	  delete thepbarRangeCoeffCTable; }
      thepbarRangeCoeffCTable = new G4PhysicsTable(numOfCouples);
      theRangeCoeffCTable = thepbarRangeCoeffCTable ;
      theRangeTable = theRangepbarTable ;
    }
  
  G4double R2 = RTable*RTable ;
  G4double R1 = RTable+1.;
  G4double w = R1*(RTable-1.)*(RTable-1.);
  G4double w1 = 1./w , w2 = -RTable*R1/w , w3 = RTable*R2/w ;
  G4double Ti , Tim , Tip , Ri , Rim , Rip , Value ;
  G4bool isOut;

  //  loop for materials
  for (G4int J=0; J<numOfCouples; J++)
    {
      G4int binmax=TotBin ;
      G4PhysicsLinearVector* aVector =
	new G4PhysicsLinearVector(0.,binmax, TotBin);
      Ti = LowestKineticEnergy ;
      G4PhysicsVector* rangeVector= (*theRangeTable)[J];
   
      for ( G4int i=0; i<TotBin; i++)
	{
	  Ri = rangeVector->GetValue(Ti,isOut) ;
	  if ( i==0 )
	    Rim = 0. ;
	  else
	    {
	      Tim = Ti/RTable ;
	      Rim = rangeVector->GetValue(Tim,isOut);
	    }
	  if ( i==(TotBin-1))
	    Rip = Ri ;
	  else
	    {
	      Tip = Ti*RTable ;
	      Rip = rangeVector->GetValue(Tip,isOut);
	    }
	  Value = w1*Rip + w2*Ri + w3*Rim ;

	  aVector->PutValue(i,Value);
	  Ti = RTable*Ti ;
	}
      theRangeCoeffCTable->insert(aVector);
    } 
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::
BuildInverseRangeTable(const G4ParticleDefinition& aParticleType)
{
  // Build inverse table of the range table

  G4bool b;

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();

  if(&aParticleType == G4Proton::Proton())
    {
      if(theInverseRangepTable)
	{ theInverseRangepTable->clearAndDestroy();
	  delete theInverseRangepTable; }
      theInverseRangepTable = new G4PhysicsTable(numOfCouples);
      theInverseRangeTable = theInverseRangepTable ;
      theRangeTable = theRangepTable ;
      theDEDXTable =  theDEDXpTable ;
      theRangeCoeffATable = thepRangeCoeffATable ;
      theRangeCoeffBTable = thepRangeCoeffBTable ;
      theRangeCoeffCTable = thepRangeCoeffCTable ;
    }

  if(&aParticleType == G4AntiProton::AntiProton())
    {
      if(theInverseRangepbarTable)
	{ theInverseRangepbarTable->clearAndDestroy();
	  delete theInverseRangepbarTable; }
      theInverseRangepbarTable = new G4PhysicsTable(numOfCouples);
      theInverseRangeTable = theInverseRangepbarTable ;
      theRangeTable = theRangepbarTable ;
      theDEDXTable =  theDEDXpbarTable ;
      theRangeCoeffATable = thepbarRangeCoeffATable ;
      theRangeCoeffBTable = thepbarRangeCoeffBTable ;
      theRangeCoeffCTable = thepbarRangeCoeffCTable ;
    }

  // loop for materials 
  for (std::size_t i=0;  i<numOfCouples; ++i)
    {

      G4PhysicsVector* pv = (*theRangeTable)[i];
      std::size_t nbins   = pv->GetVectorLength();
      G4double elow  = pv->GetLowEdgeEnergy(0);
      G4double ehigh = pv->GetLowEdgeEnergy(nbins-1);
      G4double rlow  = pv->GetValue(elow, b);
      G4double rhigh = pv->GetValue(ehigh, b);

      if (rlow <DBL_MIN) rlow = 1.e-8;
      if (rhigh > 1.e16) rhigh = 1.e16;
      if (rhigh < 1.e-8) rhigh =1.e-8;
      G4double tmpTrick = rhigh/rlow;

      //std::cout << nbins << ", elow " << elow << ", ehigh " << ehigh 
      //		<< ", rlow " << rlow << ", rhigh " << rhigh
      //                << ", trick " << tmpTrick << std::endl;

      if (tmpTrick <= 0. || tmpTrick < DBL_MIN) tmpTrick = 1.e-8; 
      if (tmpTrick > 1.e16) tmpTrick = 1.e16; 
     
      rhigh *= std::exp(std::log(tmpTrick)/((G4double)(nbins-1)));

      G4PhysicsLogVector* v = new G4PhysicsLogVector(rlow, rhigh, nbins);

      v->PutValue(0,elow);
      G4double energy1 = elow;
      G4double range1  = rlow;
      G4double energy2 = elow;
      G4double range2  = rlow;
      std::size_t ilow = 0;
      std::size_t ihigh;

      for (std::size_t j=1; j<nbins; ++j) {

	G4double range = v->GetLowEdgeEnergy(j);

	for (ihigh=ilow+1; ihigh<nbins; ihigh++) {
	  energy2 = pv->GetLowEdgeEnergy(ihigh);
	  range2  = pv->GetValue(energy2, b);
	  if(range2 >= range || ihigh == nbins-1) {
	    ilow = ihigh - 1;
	    energy1 = pv->GetLowEdgeEnergy(ilow);
	    range1  = pv->GetValue(energy1, b);
	    break;
	  }
	}

	G4double e = std::log(energy1) + std::log(energy2/energy1)*std::log(range/range1)/std::log(range2/range1);

	v->PutValue(j,std::exp(e));
      }
      theInverseRangeTable->insert(v);

    }
}

//--------------------------------------------------------------------------------

void G4hRDEnergyLoss::InvertRangeVector(G4int materialIndex,
					G4PhysicsLogVector* aVector)
{
  //  invert range vector for a material

  G4double LowEdgeRange,A,B,C,discr,KineticEnergy ;
  G4double Tbin = LowestKineticEnergy/RTable ;
  G4double rangebin = 0.0 ;
  G4int binnumber = -1 ;
  G4bool isOut ;


  //loop for range values
  for( G4int i=0; i<TotBin; i++)
    {
      LowEdgeRange = aVector->GetLowEdgeEnergy(i) ;  //i.e. GetLowEdgeValue(i)

      if( rangebin < LowEdgeRange )
	{
	  do
	    {
	      binnumber += 1 ;
	      Tbin *= RTable ;
	      rangebin = (*theRangeTable)(materialIndex)->GetValue(Tbin,isOut) ;
	    }
	  while ((rangebin < LowEdgeRange) && (binnumber < TotBin )) ;
	}

      if(binnumber == 0)
	KineticEnergy = LowestKineticEnergy ;
      else if(binnumber == TotBin-1)
	KineticEnergy = HighestKineticEnergy ;
      else
	{
	  A = (*(*theRangeCoeffATable)(materialIndex))(binnumber-1) ;
	  B = (*(*theRangeCoeffBTable)(materialIndex))(binnumber-1) ;
	  C = (*(*theRangeCoeffCTable)(materialIndex))(binnumber-1) ;
	  if(A==0.)
	    KineticEnergy = (LowEdgeRange -C )/B ;
	  else
	    {
	      discr = B*B - 4.*A*(C-LowEdgeRange);
	      discr = discr>0. ? std::sqrt(discr) : 0.;
	      KineticEnergy = 0.5*(discr-B)/A ;
	    }
	}

      aVector->PutValue(i,KineticEnergy) ;
    }
}

//------------------------------------------------------------------------------

G4bool G4hRDEnergyLoss::CutsWhereModified()
{
  G4bool wasModified = false;
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  for (G4int j=0; j<numOfCouples; ++j){
    if (theCoupleTable->GetMaterialCutsCouple(j)->IsRecalcNeeded()) {
      wasModified = true;
      break;
    }
  }
  return wasModified;
}

//-----------------------------------------------------------------------
//--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

//G4bool G4hRDEnergyLoss::IsApplicable(const G4ParticleDefinition&
//                                                     particle)
//{
//   return((particle.GetPDGCharge()!= 0.) && (particle.GetLeptonNumber() == 0));
//}
         
//G4double G4hRDEnergyLoss::GetContinuousStepLimit(
//                                               const G4Track& track,
//                                               G4double,
//                                               G4double currentMinimumStep,
//                                               G4double&)
//{
// 
//  G4double Step =
//    GetConstraints(track.GetDynamicParticle(),track.GetMaterial()) ;
//
//  if((Step>0.0)&&(Step<currentMinimumStep))
//     currentMinimumStep = Step ;
//
//  return Step ;
//}
