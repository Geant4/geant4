// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VIMuEnergyLoss.cc,v 1.3 2000-08-15 09:40:41 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      GEANT4 Collaboration
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4VIMuEnergyLoss physics process -----------
//                by Laszlo Urban, September 1997
// **************************************************************
// It is the implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of muons.
// **************************************************************
//
// corrections by L.Urban on 27/05/98 (other corrs come soon!)
// --------------------------------------------------------------
 

#include "G4VIMuEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4Poisson.hh"

// Initialisation of static members **********************************************
//  ( this stuff should be defined later using RW ..........)
// contributing processes : ion.loss,bremsstrahlung,pair production
//          ->NUMBEROFPROCESSES is initialized to 3.
//  YOU DO NOT HAVE TO CHANGE this variable for a 'normal' run.
// You have to change NUMBEROFPROCESSES  
// if you invent a new process contributing to the cont. energy loss,
//   NUMBEROFPROCESSES should be 4 in this case,
//  or for debugging purposes.
//  The NUMBEROFPROCESSES data member can be changed using the (public static)
//  functions Get/Set/Plus/MinusNUMBEROFPROCESSES (see G4VIMuEnergyLoss.hh)

G4int G4VIMuEnergyLoss::NUMBEROFPROCESSES = 3 ;
//G4int G4VIMuEnergyLoss::NUMBEROFPROCESSES = 2 ;
G4PhysicsTable** G4VIMuEnergyLoss::RecorderOfmuplusProcess =
                                           new G4PhysicsTable*[10] ;
G4int       G4VIMuEnergyLoss::CounterOfmuplusProcess = 0 ;

G4PhysicsTable* G4VIMuEnergyLoss::theDEDXmuplusTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::theRangemuplusTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::theInverseRangemuplusTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::theLabTimemuplusTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::theProperTimemuplusTable = NULL ;

G4double G4VIMuEnergyLoss::CutInmupluslossTable = 0. ;
G4double G4VIMuEnergyLoss::CutInmuminuslossTable = 0. ;

G4PhysicsTable* G4VIMuEnergyLoss::themuplusRangeCoeffATable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::themuplusRangeCoeffBTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::themuplusRangeCoeffCTable = NULL ;

G4PhysicsTable** G4VIMuEnergyLoss::RecorderOfmuminusProcess =
                                           new G4PhysicsTable*[10] ;
G4int       G4VIMuEnergyLoss::CounterOfmuminusProcess = 0 ;

G4PhysicsTable* G4VIMuEnergyLoss::theDEDXmuminusTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::theRangemuminusTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::theInverseRangemuminusTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::theLabTimemuminusTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::theProperTimemuminusTable = NULL ;

G4PhysicsTable* G4VIMuEnergyLoss::themuminusRangeCoeffATable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::themuminusRangeCoeffBTable = NULL ;
G4PhysicsTable* G4VIMuEnergyLoss::themuminusRangeCoeffCTable = NULL ;

// constructor and destructor
 
G4VIMuEnergyLoss::G4VIMuEnergyLoss(const G4String& processName)
   : G4IVContinuousDiscreteProcess (processName),
     dToverTini(0.20),   // max.relative range loss in one Step = 20%
     LowestKineticEnergy(1.00*keV),
     HighestKineticEnergy(1000000.*TeV),
     BIGSTEP ( 1.e-10*DBL_MAX ),
     MaxExcitationNumber (1.e6),
     probLimFluct (0.01),
     nmaxDirectFluct (100),
     nmaxCont1(4),
     nmaxCont2(16),
     theElectron ( G4Electron::Electron() ),
     thePositron ( G4Positron::Positron() ),
     theMuonPlus ( G4MuonPlus::MuonPlus() ),
     theMuonMinus ( G4MuonMinus::MuonMinus() )
{
     theLossTable = NULL ;
     theRangeCoeffATable = NULL ;
     theRangeCoeffBTable = NULL ;
     theRangeCoeffCTable = NULL ;
     lastMaterial = NULL ;
     lastCutInRange = 0. ;
}

G4VIMuEnergyLoss::~G4VIMuEnergyLoss() 
{
     if(theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable;
     }

}
 
 
// methods.............................................

   void G4VIMuEnergyLoss::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{
    ParticleMass = aParticleType.GetPDGMass() ; 
//-------------------------------------------------------------
//  calculate data members TotBin,LOGRTable,RTable first
  G4double lrate ;
  G4int nbin ;
//  binning  corresponds to 2.*dToverTini........................
  G4double binning = 2.*dToverTini ;

  lrate = log(HighestKineticEnergy/LowestKineticEnergy) ;
 // nbin = G4int((lrate/log(1.+dToverTini) + lrate/log(1.+2.*dToverTini))/2.);
  nbin = G4int((lrate/log(1.+binning) + lrate/log(1.+2.*binning))/2.);
  nbin = (nbin+25)/50 ; 
  TotBin = 50*nbin ;
  if(TotBin<50)
    TotBin = 50 ;
  if(TotBin>500)
    TotBin = 500 ;
  LOGRTable=lrate/TotBin;
  RTable   =exp(LOGRTable);
//--------------------------------------------------------------------
  G4bool MakeTable ;
  G4double Charge = aParticleType.GetPDGCharge() ;
  CutInRange = aParticleType.GetLengthCuts();
// Create tables only if there is a new cut value !********************************
//   and at the last contributing process only!*****************************
  if( Charge > 0.)
  {
   if(CounterOfmuplusProcess==NUMBEROFPROCESSES)
   {
    if(CutInRange != CutInmupluslossTable)
      MakeTable = true ;
      CutInmupluslossTable = CutInRange ;
   }
   else
   {
      MakeTable = false ;
   }
  }
  else
  {
   if(CounterOfmuminusProcess==NUMBEROFPROCESSES)
   {
    if(CutInRange != CutInmuminuslossTable)
      MakeTable = true ;
      CutInmuminuslossTable = CutInRange ;
   }
   else
   {
      MakeTable = false ;
   }
  }
  
  if( MakeTable )
  {

// Build energy loss table as a sum of the energy loss due to the
//              different processes.                                           
//******************************************************************
//  different processes.                                           

    const G4MaterialTable* theMaterialTable=
                                     G4Material::GetMaterialTable();

//  create table for the total energy loss

    G4int numOfMaterials = theMaterialTable->length();

  // create/fill muplus or muminus tables depending on the charge of the particle

 if( Charge >0.)    
 {
    RecorderOfProcess=RecorderOfmuplusProcess;
    CounterOfProcess=CounterOfmuplusProcess;

    if(CounterOfProcess == NUMBEROFPROCESSES)
    {
  // create tables
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

    if(CounterOfProcess == NUMBEROFPROCESSES)
    {
  // create tables
      if(theDEDXmuminusTable)
      { theDEDXmuminusTable->clearAndDestroy();
        delete theDEDXmuminusTable; }
      theDEDXmuminusTable = new G4PhysicsTable(numOfMaterials);
      theDEDXTable = theDEDXmuminusTable;
    }
  }

  if(CounterOfProcess == NUMBEROFPROCESSES)
  {
 // fill the tables

//  loop for materials
    G4double LowEdgeEnergy , Value ;
    G4bool isOutRange ;
    G4int J;

    G4PhysicsTable* pointer ;


    for (J=0; J<numOfMaterials; J++)
    {
      // create physics vector and fill it

      G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);   

      // loop for the kinetic energy
   
      for (G4int i=0; i<TotBin; i++)
  
      {
        LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;      
//     here comes the sum of the different tables created by the  
//     processes (ionisation,bremsstrahlung,pair production,etc...)              

        Value = 0. ;
    
        for (G4int process=0; process < NUMBEROFPROCESSES; process++)
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
     CounterOfmuplusProcess=0 ;
 else
     CounterOfmuminusProcess=0 ;

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
    const G4double lowestKineticEnergy(1.00*keV);
    const G4double highestKineticEnergy(1000000.*TeV);
    G4EnergyLossTables::Register(&aParticleType,  
      (Charge > 0)? theDEDXmuplusTable: theDEDXmuminusTable,
      (Charge > 0)? theRangemuplusTable: theRangemuminusTable,
      (Charge > 0)? theInverseRangemuplusTable: theInverseRangemuminusTable,
      (Charge > 0)? theLabTimemuplusTable: theLabTimemuminusTable,
      (Charge > 0)? theProperTimemuplusTable: theProperTimemuminusTable,
      lowestKineticEnergy, highestKineticEnergy, 1.,TotBin);
 
            
}
      
  void G4VIMuEnergyLoss::BuildRangeTable(
                             const G4ParticleDefinition& aParticleType)
// Build range table from the energy loss table
{
    G4double Charge = aParticleType.GetPDGCharge() ;

//  create table

    const G4MaterialTable* theMaterialTable=
                                  G4Material::GetMaterialTable();


    G4int numOfMaterials = theMaterialTable->length();

 if( Charge >0.)
 {    
    if(theRangemuplusTable)
    { theRangemuplusTable->clearAndDestroy();
      delete theRangemuplusTable; }
    theRangemuplusTable = new G4PhysicsTable(numOfMaterials);
    theRangeTable = theRangemuplusTable ;
 }
 else
 {   
    if(theRangemuminusTable)
    { theRangemuminusTable->clearAndDestroy();
      delete theRangemuminusTable; }
    theRangemuminusTable = new G4PhysicsTable(numOfMaterials);
    theRangeTable = theRangemuminusTable ;
 }

// loop for materials

    for (G4int J=0;  J<numOfMaterials; J++)
    {
    // create vector
    G4PhysicsLogVector* aVector;

    aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                            HighestKineticEnergy,TotBin);

    // fill the vector ( ranges for the actual material)

    BuildRangeVector(J, aVector);

    // insert vector to the table

    theRangeTable->insert(aVector);

    }

}    

  void G4VIMuEnergyLoss::BuildTimeTables(
                             const G4ParticleDefinition& aParticleType)
// Build time tables from the energy loss table
{


//  create table

    const G4MaterialTable* theMaterialTable=
                                  G4Material::GetMaterialTable();


    G4int numOfMaterials = theMaterialTable->length();

  if(&aParticleType == theMuonPlus)
  {
    if(theLabTimemuplusTable)
    { theLabTimemuplusTable->clearAndDestroy();
      delete theLabTimemuplusTable; }
    theLabTimemuplusTable = new G4PhysicsTable(numOfMaterials);
    theLabTimeTable = theLabTimemuplusTable ;

    if(theProperTimemuplusTable)
    { theProperTimemuplusTable->clearAndDestroy();
      delete theProperTimemuplusTable; }
    theProperTimemuplusTable = new G4PhysicsTable(numOfMaterials);
    theProperTimeTable = theProperTimemuplusTable ;
  }

  if(&aParticleType == theMuonMinus)
  {
    if(theLabTimemuminusTable)
    { theLabTimemuminusTable->clearAndDestroy();
      delete theLabTimemuminusTable; }
    theLabTimemuminusTable = new G4PhysicsTable(numOfMaterials);
    theLabTimeTable = theLabTimemuminusTable ;

    if(theProperTimemuminusTable)
    { theProperTimemuminusTable->clearAndDestroy();
      delete theProperTimemuminusTable; }
    theProperTimemuminusTable = new G4PhysicsTable(numOfMaterials);
    theProperTimeTable = theProperTimemuminusTable ;
  }

// loop for materials

    for (G4int J=0;  J<numOfMaterials; J++)
    {
    // create vector
    G4PhysicsLogVector* aVector;
    G4PhysicsLogVector* bVector;

    aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                            HighestKineticEnergy,TotBin);

    // fill the vector

    BuildLabTimeVector(J, aVector);

    // insert vector to the table

    theLabTimeTable->insert(aVector);


    bVector = new G4PhysicsLogVector(LowestKineticEnergy,
                            HighestKineticEnergy,TotBin);

    // fill the vector

    BuildProperTimeVector(J, bVector);

    // insert vector to the table

    theProperTimeTable->insert(bVector);


    }
}



void G4VIMuEnergyLoss::BuildRangeVector(G4int materialIndex,
                                     G4PhysicsLogVector* rangeVector)
//  create range vector for a material
{
  static G4int nbin;
  const G4double BigRange = DBL_MAX ;
  G4int maxbint=100;
  G4bool isOut;
  G4double tlim=2.*MeV,t1=0.1*MeV,t2=0.025*MeV ;
  G4double loss1,loss2,ca,cb,cba ;
  G4double taulim,rangelim,ltaulim,ltaumax,
           LowEdgeEnergy,tau,Value,tau1,sqtau1 ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];

  const G4MaterialTable* theMaterialTable =
                                G4Material::GetMaterialTable() ;

  // low energy part first...
  loss1 = physicsVector->GetValue(t1,isOut);
  loss2 = physicsVector->GetValue(t2,isOut);
  tau1 = t1/ParticleMass ;
  sqtau1 = sqrt(tau1) ;
  ca = (4.*loss2-loss1)/sqtau1 ;
  cb = (2.*loss1-4.*loss2)/tau1 ;
  cba = cb/ca ;
  taulim = tlim/ParticleMass ;
  ltaulim = log(taulim) ;
  ltaumax = log(HighestKineticEnergy/ParticleMass) ;


  G4int i=-1;
  G4double oldValue = 0. ;
  G4double tauold ;

  do
  {
   i += 1 ;
   LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i);

    tau = LowEdgeEnergy/ParticleMass;
    if ( tau <= tau1 )
    {
      Value = 2.*ParticleMass*log(1.+cba*sqrt(tau))/cb ;
    }
    else
    {
      Value = 2.*ParticleMass*log(1.+cba*sqtau1)/cb ;
      if(tau<=taulim)
      {
        nbin = maxbint ;
        taulow = tau1 ;
        tauhigh = tau ;
        Value += RangeIntLin(physicsVector,nbin);
      }
      else
      {
        taulow = tau1 ;
        tauhigh = taulim ; 
        Value += RangeIntLin(physicsVector,maxbint) ;
        ltaulow = ltaulim ;
        ltauhigh = log(tau) ;
        nbin = maxbint ;
        Value += RangeIntLog(physicsVector,nbin);
      }
    }


      rangeVector->PutValue(i,Value);
      oldValue = Value ;
      tauold = tau ;

  } while (tau<=taulim) ;

  i += 1 ;

  for (G4int j=i; j<TotBin; j++)
  {
   LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(j);
    tau = LowEdgeEnergy/ParticleMass;
      ltaulow = log(tauold);
      ltauhigh = log(tau);
      nbin = maxbint;
      Value = oldValue+RangeIntLog(physicsVector,nbin);

      rangeVector->PutValue(j,Value);
      oldValue = Value ;
      tauold = tau ;

  }

}    

void G4VIMuEnergyLoss::BuildLabTimeVector(G4int materialIndex,
                                     G4PhysicsLogVector* timeVector)
//  create lab time vector for a material
{

  static G4int nbin;
  G4int maxbint=100;
  G4bool isOut;
  G4double tlim=5.*keV,parlowen=0.4,ppar=0.5-parlowen ;
  G4double losslim,clim,taulim,timelim,ltaulim,ltaumax,
           LowEdgeEnergy,tau,Value ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];

  const G4MaterialTable* theMaterialTable =
                                G4Material::GetMaterialTable() ;

  // low energy part first...
  losslim = physicsVector->GetValue(tlim,isOut);
  taulim=tlim/ParticleMass ;
  clim=sqrt(ParticleMass*tlim/2.)/(c_light*losslim*ppar) ;
  ltaulim = log(taulim);
  ltaumax = log(HighestKineticEnergy/ParticleMass) ;

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
      Value = clim*exp(ppar*log(tau/taulim)) ;
    }
    else
    {
      timelim=clim ;
      ltaulow = log(taulim);
      ltauhigh = log(tau);
      nbin = maxbint;
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
      ltaulow = log(tauold);
      ltauhigh = log(tau);
      nbin = maxbint ;
      Value = oldValue+LabTimeIntLog(physicsVector,nbin);

      timeVector->PutValue(j,Value);
      oldValue = Value ;
      tauold = tau ;
  }

}

void G4VIMuEnergyLoss::BuildProperTimeVector(G4int materialIndex,
                                     G4PhysicsLogVector* timeVector)
//  create proper time vector for a material
{

  static G4int nbin;
  G4int maxbint=100;
  G4bool isOut;
  G4double tlim=5.*keV,parlowen=0.4,ppar=0.5-parlowen ;
  G4double losslim,clim,taulim,timelim,ltaulim,ltaumax,
           LowEdgeEnergy,tau,Value ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];

  const G4MaterialTable* theMaterialTable =
                                G4Material::GetMaterialTable() ;

  // low energy part first...
  losslim = physicsVector->GetValue(tlim,isOut);
  taulim=tlim/ParticleMass ;
  clim=sqrt(ParticleMass*tlim/2.)/(c_light*losslim*ppar) ;
  ltaulim = log(taulim);
  ltaumax = log(HighestKineticEnergy/ParticleMass) ;

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
      Value = clim*exp(ppar*log(tau/taulim)) ;
    }
    else
    {
      timelim=clim ;
      ltaulow = log(taulim);
      ltauhigh = log(tau);
      nbin = maxbint;
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
      ltaulow = log(tauold);
      ltauhigh = log(tau);
      nbin = maxbint ;
      Value = oldValue+ProperTimeIntLog(physicsVector,nbin);

      timeVector->PutValue(j,Value);
      oldValue = Value ;
      tauold = tau ;
  }

}


G4double G4VIMuEnergyLoss::RangeIntLin(G4PhysicsVector* physicsVector,
                                    G4int nbin)
//  num. integration, linear binning
{
  G4double dtau,Value,taui,ti,lossi,ci;
  G4bool isOut;

  dtau = (tauhigh-taulow)/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
  {
    taui = taulow + dtau*i ;
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

    Value += ci/lossi;
  }

  Value *= ParticleMass*dtau;

  return Value;

}


G4double G4VIMuEnergyLoss::RangeIntLog(G4PhysicsVector* physicsVector,
                                    G4int nbin)
//  num. integration, logarithmic binning
{
  G4double ltt,dltau,Value,ui,taui,ti,lossi,ci;
  G4bool isOut;

  ltt = ltauhigh-ltaulow;
  dltau = ltt/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
  {
    ui = ltaulow+dltau*i;
    taui = exp(ui);
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

    Value += ci*taui/lossi;
  }

  Value *= ParticleMass*dltau;

  return Value;

}

G4double G4VIMuEnergyLoss::LabTimeIntLog(G4PhysicsVector* physicsVector,
                                    G4int nbin)
//  num. integration, logarithmic binning
{
  G4double ltt,dltau,Value,ui,taui,ti,lossi,ci;
  G4bool isOut;

  ltt = ltauhigh-ltaulow;
  dltau = ltt/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
  {
    ui = ltaulow+dltau*i;
    taui = exp(ui);
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

    Value += ci*taui*(ti+ParticleMass)/(sqrt(ti*(ti+2.*ParticleMass))*lossi);
  }

  Value *= ParticleMass*dltau/c_light;

  return Value;

}

G4double G4VIMuEnergyLoss::ProperTimeIntLog(G4PhysicsVector* physicsVector,
                                    G4int nbin)
//  num. integration, logarithmic binning
{
  G4double ltt,dltau,Value,ui,taui,ti,lossi,ci;
  G4bool isOut;

  ltt = ltauhigh-ltaulow;
  dltau = ltt/nbin;
  Value = 0.;

  for (G4int i=0; i<=nbin; i++)
  {
    ui = ltaulow+dltau*i;
    taui = exp(ui);
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

    Value += ci*taui*ParticleMass/(sqrt(ti*(ti+2.*ParticleMass))*lossi);
  }

  Value *= ParticleMass*dltau/c_light;

  return Value;

}


void G4VIMuEnergyLoss::BuildRangeCoeffATable(
                            const G4ParticleDefinition& aParticleType)
// Build tables of coefficients for the energy loss calculation
{
   G4double Charge = aParticleType.GetPDGCharge() ;

   const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();

//  create table for coefficients "A"

  G4int numOfMaterials = theMaterialTable->length();

  if(Charge>0.)
  {
    if(themuplusRangeCoeffATable)
    { themuplusRangeCoeffATable->clearAndDestroy();
      delete themuplusRangeCoeffATable; }
    themuplusRangeCoeffATable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffATable = themuplusRangeCoeffATable ;
    theRangeTable = theRangemuplusTable ;
  }
  else  
  {
    if(themuminusRangeCoeffATable)
    { themuminusRangeCoeffATable->clearAndDestroy();
      delete themuminusRangeCoeffATable; }
    themuminusRangeCoeffATable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffATable = themuminusRangeCoeffATable ;
    theRangeTable = theRangemuminusTable ;
  }
 
  G4double R2 = RTable*RTable ;
  G4double R1 = RTable+1.;
  G4double w = R1*(RTable-1.)*(RTable-1.);
  G4double w1 = RTable/w , w2 = -RTable*R1/w , w3 = R2/w ;
  G4double Ti , Tim , Tip , Ri , Rim , Rip , Value ;
  G4bool isOut;

  //  loop for materials
  for (G4int J=0; J<numOfMaterials; J++)
  {

  // create vector
   G4int binmax=TotBin ;
   G4PhysicsLinearVector* aVector =  new G4PhysicsLinearVector(0.,binmax, TotBin);

   //  loop for kinetic energy
   
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

     Value = (w1*Rip + w2*Ri + w3*Rim)/(Ti*Ti) ; 

     aVector->PutValue(i,Value);

     Ti = RTable*Ti ;

   }
  
   theRangeCoeffATable->insert(aVector);

  } 
}


void G4VIMuEnergyLoss::BuildRangeCoeffBTable(
                            const G4ParticleDefinition& aParticleType)
// Build tables of coefficients for the energy loss calculation
{
   G4double Charge = aParticleType.GetPDGCharge() ;

   const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();

//  create table for coefficients "B"


  G4int numOfMaterials = theMaterialTable->length();

  if(Charge>0.)
  {
    if(themuplusRangeCoeffBTable)
    { themuplusRangeCoeffBTable->clearAndDestroy();
      delete themuplusRangeCoeffBTable; }
    themuplusRangeCoeffBTable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffBTable = themuplusRangeCoeffBTable ;
    theRangeTable = theRangemuplusTable ;
  }
  else
  {
    if(themuminusRangeCoeffBTable)
    { themuminusRangeCoeffBTable->clearAndDestroy();
      delete themuminusRangeCoeffBTable; }
    themuminusRangeCoeffBTable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffBTable = themuminusRangeCoeffBTable ;
    theRangeTable = theRangemuminusTable ;
  }

  G4double R2 = RTable*RTable ;
  G4double R1 = RTable+1.;
  G4double w = R1*(RTable-1.)*(RTable-1.);
  G4double w1 = -R1/w , w2 = R1*(R2+1.)/w , w3 = -R2*R1/w ;
  G4double Ti , Tim , Tip , Ri , Rim , Rip , Value ;
  G4bool isOut;

  //  loop for materials
  for (G4int J=0; J<numOfMaterials; J++)
  {

  // create  vector
   G4int binmax=TotBin ;

   G4PhysicsLinearVector* aVector =  new G4PhysicsLinearVector(0.,binmax, TotBin);

   //  loop for kinetic energy

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

     Value = (w1*Rip + w2*Ri + w3*Rim)/Ti;  

     aVector->PutValue(i,Value);

     Ti = RTable*Ti ;

   }
  
   theRangeCoeffBTable->insert(aVector);

  } 
}
void G4VIMuEnergyLoss::BuildRangeCoeffCTable(
                            const G4ParticleDefinition& aParticleType)
// Build tables of coefficients for the energy loss calculation
{
   G4double Charge = aParticleType.GetPDGCharge() ;
   const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();

//  create table for coefficients "C"

  G4int numOfMaterials = theMaterialTable->length();

  if(Charge>0.)
  {
    if(themuplusRangeCoeffCTable)
    { themuplusRangeCoeffCTable->clearAndDestroy();
      delete themuplusRangeCoeffCTable; }
    themuplusRangeCoeffCTable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffCTable = themuplusRangeCoeffCTable ;
    theRangeTable = theRangemuplusTable ;
  }
  else
  {
    if(themuminusRangeCoeffCTable)
    { themuminusRangeCoeffCTable->clearAndDestroy();
      delete themuminusRangeCoeffCTable; }
    themuminusRangeCoeffCTable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffCTable = themuminusRangeCoeffCTable ;
    theRangeTable = theRangemuminusTable ;
  }
  
  const G4double BigRange = DBL_MAX ;
  G4double R2 = RTable*RTable ;
  G4double R1 = RTable+1.;
  G4double w = R1*(RTable-1.)*(RTable-1.);
  G4double w1 = 1./w , w2 = -RTable*R1/w , w3 = RTable*R2/w ;
  G4double Ti , Tim , Tip , Ri , Rim , Rip , Value ;
  G4bool isOut;

  //  loop for materials
  for (G4int J=0; J<numOfMaterials; J++)
  {

  // create  vector
   G4int binmax=TotBin ;
   G4PhysicsLinearVector* aVector =  new G4PhysicsLinearVector(0.,binmax, TotBin);

   //  loop for kinetic energy

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

  void G4VIMuEnergyLoss::BuildInverseRangeTable(
                             const G4ParticleDefinition& aParticleType)
// Build inverse table of the range table
{
    G4double SmallestRange,BiggestRange ;
    G4bool isOut ;

//  create table

    const G4MaterialTable* theMaterialTable=
                                  G4Material::GetMaterialTable();

    G4int numOfMaterials = theMaterialTable->length();

  if(&aParticleType == theMuonPlus)
  {
    if(theInverseRangemuplusTable)
    { theInverseRangemuplusTable->clearAndDestroy();
      delete theInverseRangemuplusTable; }
    theInverseRangemuplusTable = new G4PhysicsTable(numOfMaterials);
    theInverseRangeTable = theInverseRangemuplusTable ;
    theRangeTable = theRangemuplusTable ;
    theDEDXTable =  theDEDXmuplusTable ;
    theRangeCoeffATable = themuplusRangeCoeffATable ;
    theRangeCoeffBTable = themuplusRangeCoeffBTable ;
    theRangeCoeffCTable = themuplusRangeCoeffCTable ;
  }

  if(&aParticleType == theMuonMinus)
  {
    if(theInverseRangemuminusTable)
    { theInverseRangemuminusTable->clearAndDestroy();
      delete theInverseRangemuminusTable; }
    theInverseRangemuminusTable = new G4PhysicsTable(numOfMaterials);
    theInverseRangeTable = theInverseRangemuminusTable ;
    theRangeTable = theRangemuminusTable ;
    theDEDXTable =  theDEDXmuminusTable ;
    theRangeCoeffATable = themuminusRangeCoeffATable ;
    theRangeCoeffBTable = themuminusRangeCoeffBTable ;
    theRangeCoeffCTable = themuminusRangeCoeffCTable ;
  }
// loop for materials

    for (G4int J=0;  J<numOfMaterials; J++)
    {
    SmallestRange = (*theRangeTable)(J)->
                       GetValue(LowestKineticEnergy,isOut) ;
    BiggestRange = (*theRangeTable)(J)->
                       GetValue(HighestKineticEnergy,isOut) ;
    // create vector
    G4PhysicsLogVector* aVector;

    aVector = new G4PhysicsLogVector(SmallestRange,
                            BiggestRange,TotBin);

    // fill the vector ( ranges for the actual material)

    InvertRangeVector(J, aVector);

    // insert vector to the table

    theInverseRangeTable->insert(aVector);

    }
}

void G4VIMuEnergyLoss::InvertRangeVector(G4int materialIndex,
                                     G4PhysicsLogVector* aVector)
//  invert range vector for a material
{
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
        discr = discr>0. ? sqrt(discr) : 0.;
        KineticEnergy = 0.5*(discr-B)/A ;
    }
  }

  aVector->PutValue(i,KineticEnergy) ;

 }
}



G4VParticleChange* G4VIMuEnergyLoss::AlongStepDoIt( 
                                    const G4Track& trackData,const G4Step& stepData) 
 // compute the energy loss after a Step 
{
  const G4DynamicParticle* aParticle;
  G4Material* aMaterial;
  G4bool isOut;
  G4double E,finalT,Step,Tbin,rangebin ;
  const G4double smallLoss=DBL_MIN;
  const G4double BigRange = DBL_MAX ;
  G4int index ;
  G4double cc,discr ;
  G4double Charge  ; 

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;

 // get the actual (true) Step length from stepData 
  Step = stepData.GetStepLength() ;


 // there is no loss for Step=0. !
  if( Step == 0.)
    return &aParticleChange ;  

 // get particle and material pointers from trackData 
  aParticle = trackData.GetDynamicParticle() ;
  Charge = aParticle->GetDefinition()->GetPDGCharge() ;
  index = aMaterial->GetIndex() ;
  E = aParticle->GetKineticEnergy() ;

  if(Charge>0.)
  {
    theRangeTable=theRangemuplusTable;
    theRangeCoeffATable = themuplusRangeCoeffATable ;
    theRangeCoeffBTable = themuplusRangeCoeffBTable ;
    theRangeCoeffCTable = themuplusRangeCoeffCTable ;
  }
  else
  {
    theRangeTable=theRangemuminusTable;
    theRangeCoeffATable = themuminusRangeCoeffATable ;
    theRangeCoeffBTable = themuminusRangeCoeffBTable ;
    theRangeCoeffCTable = themuminusRangeCoeffCTable ;
  }

 // 
    ParticleCutInKineticEnergyNow =
          (aParticle->GetDefinition()->GetEnergyCuts())[index] ;

  if(Step >= BigRange)
  {
   finalT = E ;
   fMeanLoss = 0. ;
  }
  else
  // here comes the 'real' energy loss calculation (material is NOT vacuum)
  {

  if( E < LowestKineticEnergy)
    {
      finalT = 0.0;
      fMeanLoss = E ;
    }
  else 
    {
      if( E > HighestKineticEnergy)
      {
        finalT = E - smallLoss ; 
        fMeanLoss = smallLoss ;
      }
      else
      {
            

  //   loss calculation with quadratic interpolation in the table
  if (Step >= (fRangeNow-CutInRange))
  {
    finalT = 0.;
    fMeanLoss = E ;
  }
  else
  {
  //..........................................................................
  // check if the energy bin has changed 

    Tbin = LowestKineticEnergy*exp(EnergyBinNumber*LOGRTable) ;

    rangebin = (*theRangeTable)(index)->GetValue(Tbin,isOut) ;

    if((fRangeNow-Step)<rangebin)
    {
     do
     {
      EnergyBinNumber-- ;
      Tbin /= RTable ;
      rangebin = (*theRangeTable)(index)->GetValue(Tbin,isOut) ;
     }
     while (((fRangeNow-Step)<rangebin)&&(EnergyBinNumber>0)) ;


      RangeCoeffA = (*(*theRangeCoeffATable)(index))(EnergyBinNumber) ;
      RangeCoeffB = (*(*theRangeCoeffBTable)(index))(EnergyBinNumber) ;
      RangeCoeffC = (*(*theRangeCoeffCTable)(index))(EnergyBinNumber) ;
    }

  //..........................................................................
  //  now the energy loss can be calculated
  //  first the mean loss

    cc=Step+RangeCoeffC-fRangeNow ;
    discr =  RangeCoeffB*RangeCoeffB-4.*RangeCoeffA*cc ;
    discr = discr<=0. ? 0. : sqrt(discr) ;
    fMeanLoss = E-0.5*(discr-RangeCoeffB)/RangeCoeffA ;

  //  now the loss with fluctuation
    finalT = E-GetLossWithFluct(aParticle,aMaterial,1.,fMeanLoss,Step) ;
    if (finalT < 0.) finalT = 0. ;
   }
  }
  }
 }


  //  kill the particle if the kinetic energy <= 0  

  if (finalT <= 0. )
    {
     finalT = 0.;
     aParticleChange.SetStatusChange(fStopButAlive);
    } 

  aParticleChange.SetNumberOfSecondaries(0);
  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(E-finalT) ;

  return &aParticleChange ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VIMuEnergyLoss::GetLossWithFluct(const G4DynamicParticle* aParticle,
	                                  G4Material* aMaterial,
	                                 G4double ChargeSquare,
	                                 G4double MeanLoss,
	                                  G4double step )

{
//  calculate actual loss from the mean loss
//  The model used to get the fluctuation is essentially the same as in Glandz in Geant3.

   static const G4double minLoss = 1.*eV ;
   static const G4double probLim = 0.01 ;
   static const G4double sumaLim = -log(probLim) ;
   static const G4double alim=10.;
   static const G4double kappa = 10. ;
   static const G4double factor = twopi_mc2_rcl2 ;


  // check if the material has changed ( cache mechanism)

  if (aMaterial != lastMaterial)
    {
	lastMaterial = aMaterial;
	imat	   = aMaterial->GetIndex();
	f1Fluct	   = aMaterial->GetIonisation()->GetF1fluct();
	f2Fluct	   = aMaterial->GetIonisation()->GetF2fluct();
	e1Fluct	   = aMaterial->GetIonisation()->GetEnergy1fluct();
	e2Fluct	   = aMaterial->GetIonisation()->GetEnergy2fluct();
	e1LogFluct   = aMaterial->GetIonisation()->GetLogEnergy1fluct();
	e2LogFluct   = aMaterial->GetIonisation()->GetLogEnergy2fluct();
	rateFluct	   = aMaterial->GetIonisation()->GetRateionexcfluct();
	ipotFluct	   = aMaterial->GetIonisation()->GetMeanExcitationEnergy();
	ipotLogFluct = aMaterial->GetIonisation()->GetLogMeanExcEnergy();
    }

	G4double threshold,w1,w2,C,
		beta2,suma,e0,loss,lossc ,w,electronDensity;
  G4double a1,a2,a3;
  G4int p1,p2,p3;
  G4int nb;
  G4double Corrfac, na,alfa,rfac,namean,sa,alfa1,ea,sea;
  G4double dp1,dp3;
  G4double siga ;

  // shortcut for very very small loss
  if(MeanLoss < minLoss) return MeanLoss ;

  // get particle data
  G4double Tkin	  = aParticle->GetKineticEnergy();
  ParticleMass = aParticle->GetMass() ;

  threshold =((*G4Electron::Electron()).GetCutsInEnergy())[imat];
  G4double rmass = electron_mass_c2/ParticleMass;
  G4double tau	 = Tkin/ParticleMass, tau1 = tau+1., tau2 = tau*(tau+2.);
  G4double Tm	 = 2.*electron_mass_c2*tau2/(1.+2.*tau1*rmass+rmass*rmass);

  if (Tm <= ipotFluct) Tm = ipotFluct ;

  if(Tm > threshold) Tm = threshold;
  beta2 = tau2/(tau1*tau1);

  // Gaussian fluctuation ?
  if(MeanLoss >= kappa*Tm)
  {
      electronDensity = aMaterial->GetElectronDensity() ;
      siga = sqrt(MeanLoss*Tm*(0.5-0.25*beta2)*step*
             factor*electronDensity*ChargeSquare/beta2) ;
      loss = G4RandGauss::shoot(MeanLoss,siga) ;
      if(loss < 0.) loss = 0. ;
      return loss ;
  }

  w1 = Tm/ipotFluct;
  w2 = log(2.*electron_mass_c2*tau2);

  C = MeanLoss*(1.-rateFluct)/(w2-ipotLogFluct-beta2);

  a1 = C*f1Fluct*(w2-e1LogFluct-beta2)/e1Fluct;
  a2 = C*f2Fluct*(w2-e2LogFluct-beta2)/e2Fluct;
  if(Tm > ipotFluct)
	a3 = rateFluct*MeanLoss*(Tm-ipotFluct)/(ipotFluct*Tm*log(w1));
  else
  {
	   a1 /= 1.-rateFluct ;
	a2 /= 1.-rateFluct ;
	a3	 = 0. ;
  }

  suma = a1+a2+a3;

  loss = 0. ;

  if(suma < sumaLim)		 // very small Step
    {
		  e0 = aMaterial->GetIonisation()->GetEnergy0fluct();

	if(Tm == ipotFluct)
      {
		a3 = MeanLoss/e0;

	if(a3>alim)
        {
	  siga=sqrt(a3) ;
	  p3 = G4std::max(0,int(G4RandGauss::shoot(a3,siga)+0.5));
        }
	  p3 = G4Poisson(a3);

	loss = p3*e0 ;

	if(p3 > 0)
		loss += (1.-2.*G4UniformRand())*e0 ;

	}
      else
      {
		Tm = Tm-ipotFluct+e0 ;
	a3 = MeanLoss*(Tm-e0)/(Tm*e0*log(Tm/e0));

	if(a3>alim)
        {
	  siga=sqrt(a3) ;
		p3 = G4std::max(0,int(G4RandGauss::shoot(a3,siga)+0.5));
	}
	else
	  p3 = G4Poisson(a3);

	if(p3 > 0)
        {
	  w = (Tm-e0)/Tm ;
	  if(p3 > nmaxCont2)
	  {
	dp3 = G4float(p3) ;
		Corrfac = dp3/G4float(nmaxCont2) ;
	p3 = nmaxCont2 ;
	  }
	  else
	Corrfac = 1. ;

	  for(G4int i=0; i<p3; i++) loss += 1./(1.-w*G4UniformRand()) ;
	  loss *= e0*Corrfac ;
        }
      }
    }

	else			// not so small Step

    {
	  // excitation type 1
	if(a1>alim)
      {
		siga=sqrt(a1) ;
	p1 = G4std::max(0,int(G4RandGauss::shoot(a1,siga)+0.5));
      }
      else
	p1 = G4Poisson(a1);

	// excitation type 2
	if(a2>alim)
      {
		siga=sqrt(a2) ;
	p2 = G4std::max(0,int(G4RandGauss::shoot(a2,siga)+0.5));
      }
      else
	p2 = G4Poisson(a2);

	loss = p1*e1Fluct+p2*e2Fluct;

	// smearing to avoid unphysical peaks
	if(p2 > 0)
	loss += (1.-2.*G4UniformRand())*e2Fluct;
	else if (loss>0.)
	loss += (1.-2.*G4UniformRand())*e1Fluct;

	// ionisation .......................................
	if(a3 > 0.)
     {
	 if(a3>alim)
      {
		siga=sqrt(a3) ;
	p3 = G4std::max(0,int(G4RandGauss::shoot(a3,siga)+0.5));
      }
      else

	p3 = G4Poisson(a3);

	lossc = 0.;
	if(p3 > 0)
      {
		na = 0.;
	alfa = 1.;
	if (p3 > nmaxCont2)
        {
	  dp3		= G4float(p3);
		rfac		= dp3/(G4float(nmaxCont2)+dp3);
	namean	= G4float(p3)*rfac;
		sa		= G4float(nmaxCont1)*rfac;
		 na		= G4RandGauss::shoot(namean,sa);
	if (na > 0.)
	  {
		  alfa	 = w1*G4float(nmaxCont2+p3)/(w1*G4float(nmaxCont2)+G4float(p3));
		 alfa1  = alfa*log(alfa)/(alfa-1.);
	ea	   = na*ipotFluct*alfa1;
		sea	   = ipotFluct*sqrt(na*(alfa-alfa1*alfa1));
		 lossc += G4RandGauss::shoot(ea,sea);
	}
        }

	nb = G4int(G4float(p3)-na);
	if (nb > 0)
        {
	  w2 = alfa*ipotFluct;
	  w  = (Tm-w2)/Tm;
		for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());

	}
      }		
	loss += lossc;
     }
    }

  return loss ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




