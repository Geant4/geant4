// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VIhEnergyLoss.cc,v 1.3 2000-08-15 09:42:45 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// -----------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      GEANT4 Collaboration
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4VIhEnergyLoss physics process -----------
//                by Laszlo Urban, 30 May 1997 
//
// **************************************************************
// It is the first implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of charged hadrons.
// **************************************************************
//
// 7/10/98: bug fixes + some cleanup , L.Urban 
// 26/10/98 : cleanup + TOF tables, L.Urban
// --------------------------------------------------------------

#include "G4VIhEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4Poisson.hh"

// Initialisation of static members ******************************************
// contributing processes : ion.loss ->NumberOfProcesses is initialized
//   to 1 . YOU DO NOT HAVE TO CHANGE this variable for a 'normal' run.
// You have to change NumberOfProcesses  
// if you invent a new process contributing to the cont. energy loss,
//   NumberOfProcesses should be 2 in this case,
//  or for debugging purposes.
//  The NumberOfProcesses data member can be changed using the (public static)
//  functions Get/Set/Plus/MinusNumberOfProcesses (see G4VIhEnergyLoss.hh)

G4int G4VIhEnergyLoss::NumberOfProcesses = 1 ;

G4int            G4VIhEnergyLoss::CounterOfProcess = 0 ;
G4PhysicsTable** G4VIhEnergyLoss::RecorderOfProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4VIhEnergyLoss::CounterOfpProcess = 0 ;
G4PhysicsTable** G4VIhEnergyLoss::RecorderOfpProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4VIhEnergyLoss::CounterOfpbarProcess = 0 ;
G4PhysicsTable** G4VIhEnergyLoss::RecorderOfpbarProcess =
                                           new G4PhysicsTable*[10] ;

G4PhysicsTable* G4VIhEnergyLoss::theDEDXpTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theDEDXpbarTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theRangepTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theRangepbarTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theInverseRangepTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theInverseRangepbarTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theLabTimepTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theLabTimepbarTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theProperTimepTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theProperTimepbarTable = NULL ;

G4PhysicsTable* G4VIhEnergyLoss::thepRangeCoeffATable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::thepRangeCoeffBTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::thepRangeCoeffCTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::thepbarRangeCoeffATable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::thepbarRangeCoeffBTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::thepbarRangeCoeffCTable = NULL ;

G4PhysicsTable* G4VIhEnergyLoss::theDEDXTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theRangeTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theInverseRangeTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theLabTimeTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theProperTimeTable = NULL ;

G4PhysicsTable* G4VIhEnergyLoss::theRangeCoeffATable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theRangeCoeffBTable = NULL ;
G4PhysicsTable* G4VIhEnergyLoss::theRangeCoeffCTable = NULL ;

const G4Proton*     G4VIhEnergyLoss::theProton    =G4Proton::Proton() ;
const G4AntiProton* G4VIhEnergyLoss::theAntiProton=G4AntiProton::AntiProton() ;

G4double G4VIhEnergyLoss::ParticleMass;

G4double G4VIhEnergyLoss::Mass,
         G4VIhEnergyLoss::taulow, 
         G4VIhEnergyLoss::tauhigh, 
         G4VIhEnergyLoss::ltaulow, 
         G4VIhEnergyLoss::ltauhigh; 

G4double G4VIhEnergyLoss::CutInRange = 0. ;

G4double G4VIhEnergyLoss::dRoverRange = 0.20 ;
G4double G4VIhEnergyLoss::finalRange = 200.*micrometer ;

G4bool   G4VIhEnergyLoss::rndmStepFlag   = false ;
G4bool   G4VIhEnergyLoss::EnlossFlucFlag = true ;

G4double G4VIhEnergyLoss::LowestKineticEnergy;
G4double G4VIhEnergyLoss::HighestKineticEnergy;
G4int G4VIhEnergyLoss::TotBin;
G4double G4VIhEnergyLoss::RTable,G4VIhEnergyLoss::LOGRTable;

// constructor and destructor
 
G4VIhEnergyLoss::G4VIhEnergyLoss(const G4String& processName)
   : G4IVContinuousDiscreteProcess (processName),
     theLossTable (NULL),
     lastMaterial (NULL),
     MaxExcitationNumber (1.e6),
     probLimFluct (0.01),
     nmaxDirectFluct (100),
     nmaxCont1(4),
     nmaxCont2(16) 
{ }

G4VIhEnergyLoss::~G4VIhEnergyLoss() 
{
     if(theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable;
     }
}
 
void G4VIhEnergyLoss::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{
  //  calculate data members TotBin,LOGRTable,RTable first
  G4double lrate ;
  G4int nbin ;
  G4double binning = 2.*dRoverRange;
  lrate = log(HighestKineticEnergy/LowestKineticEnergy) ;
  nbin = G4int((lrate/log(1.+binning) + lrate/log(1.+2.*binning))/2.);
  nbin = (nbin+25)/50 ;
  TotBin = 50*nbin ;
  if(TotBin<50)
    TotBin = 50 ;
  if(TotBin>500)
    TotBin = 500 ;
  LOGRTable=lrate/TotBin;
  RTable   =exp(LOGRTable);

  // create table if there is no table or there is a new cut value
  G4bool MakeTable = false ;
     
  ParticleMass = aParticleType.GetPDGMass() ;

  G4double newCutInRange = aParticleType.GetLengthCuts();

  // create/fill proton or antiproton tables depending on the charge 
  G4double Charge = aParticleType.GetPDGCharge();

  if (Charge>0.) {theDEDXTable= theDEDXpTable;}
  else           {theDEDXTable= theDEDXpbarTable;}

  if ((CutInRange != newCutInRange) || (theDEDXTable==NULL)) 
  {
    MakeTable = true ;
    CutInRange = newCutInRange ;
  }
  
  if( MakeTable )
  {

  // Build energy loss table as a sum of the energy loss due to the
  //              different processes.                                           
    const G4MaterialTable* theMaterialTable=
                                     G4Material::GetMaterialTable();
    G4int numOfMaterials = theMaterialTable->length();

    if( Charge >0.)    
    {
      RecorderOfProcess=RecorderOfpProcess;
      CounterOfProcess=CounterOfpProcess;

      if(CounterOfProcess == NumberOfProcesses)
      {
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

      if(CounterOfProcess == NumberOfProcesses)
      {
        if(theDEDXpbarTable)
        { theDEDXpbarTable->clearAndDestroy();
          delete theDEDXpbarTable; }
        theDEDXpbarTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXpbarTable;
      }
    }

    if(CounterOfProcess == NumberOfProcesses)
    {
      //  loop for materials
      G4double LowEdgeEnergy , Value ;
      G4bool isOutRange ;
      G4PhysicsTable* pointer ;

      for (G4int J=0; J<numOfMaterials; J++)
      { 
        // create physics vector and fill it
        G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);   

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
  const G4double lowestKineticEnergy(1.00*keV);
  const G4double highestKineticEnergy(100.*TeV);

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
    lowestKineticEnergy, highestKineticEnergy,
    proton_mass_c2/aParticleType.GetPDGMass(),TotBin);
}
      
void G4VIhEnergyLoss::BuildRangeTable(
                             const G4ParticleDefinition& aParticleType)
// Build range table from the energy loss table
{
   G4double Charge = aParticleType.GetPDGCharge() ;
   Mass = proton_mass_c2; 
   const G4MaterialTable* theMaterialTable=
                                 G4Material::GetMaterialTable();
   G4int numOfMaterials = theMaterialTable->length();

   if( Charge >0.)
   {    
     if(theRangepTable)
     { theRangepTable->clearAndDestroy();
       delete theRangepTable; }
     theRangepTable = new G4PhysicsTable(numOfMaterials);
     theRangeTable = theRangepTable ;
   }
   else
   {   
     if(theRangepbarTable)
     { theRangepbarTable->clearAndDestroy();
       delete theRangepbarTable; }
     theRangepbarTable = new G4PhysicsTable(numOfMaterials);
     theRangeTable = theRangepbarTable ;
   }

   // loop for materials

   for (G4int J=0;  J<numOfMaterials; J++)
   {
     G4PhysicsLogVector* aVector;
     aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                              HighestKineticEnergy,TotBin);

     BuildRangeVector(J, aVector);
     theRangeTable->insert(aVector);
   }
}    

void G4VIhEnergyLoss::BuildTimeTables(
                             const G4ParticleDefinition& aParticleType)
{
  const G4MaterialTable* theMaterialTable=
                                 G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();
  
  if(&aParticleType == G4Proton::Proton())
  {
    if(theLabTimepTable)
    { theLabTimepTable->clearAndDestroy();
      delete theLabTimepTable; }
    theLabTimepTable = new G4PhysicsTable(numOfMaterials);
    theLabTimeTable = theLabTimepTable ;

    if(theProperTimepTable)
    { theProperTimepTable->clearAndDestroy();
      delete theProperTimepTable; }
    theProperTimepTable = new G4PhysicsTable(numOfMaterials);
    theProperTimeTable = theProperTimepTable ;
  }

  if(&aParticleType == G4AntiProton::AntiProton())
  {
    if(theLabTimepbarTable)
    { theLabTimepbarTable->clearAndDestroy();
      delete theLabTimepbarTable; }
    theLabTimepbarTable = new G4PhysicsTable(numOfMaterials);
    theLabTimeTable = theLabTimepbarTable ;

    if(theProperTimepbarTable)
    { theProperTimepbarTable->clearAndDestroy();
      delete theProperTimepbarTable; }
    theProperTimepbarTable = new G4PhysicsTable(numOfMaterials);
    theProperTimeTable = theProperTimepbarTable ;
  }

  for (G4int J=0;  J<numOfMaterials; J++)
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

void G4VIhEnergyLoss::BuildRangeVector(G4int materialIndex,
                                     G4PhysicsLogVector* rangeVector)
//  create range vector for a material
{
  G4int nbin=100;
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
  tau1 = t1/Mass ;
  sqtau1 = sqrt(tau1) ;
  ca = (4.*loss2-loss1)/sqtau1 ;
  cb = (2.*loss1-4.*loss2)/tau1 ;
  cba = cb/ca ;
  taulim = tlim/Mass ;
  ltaulim = log(taulim) ;
  ltaumax = log(HighestKineticEnergy/Mass) ;

  G4int i=-1;
  G4double oldValue = 0. ;
  G4double tauold ;
  do
  {
    i += 1 ;
    LowEdgeEnergy = rangeVector->GetLowEdgeEnergy(i);
    tau = LowEdgeEnergy/Mass;
    if ( tau <= tau1 )
    {
      Value = 2.*Mass*log(1.+cba*sqrt(tau))/cb ;
    }
    else
    {
      Value = 2.*Mass*log(1.+cba*sqtau1)/cb ;
      if(tau<=taulim)
      {
        taulow = tau1 ;
        tauhigh = tau ;
        Value += RangeIntLin(physicsVector,nbin);
      }
      else
      {
        taulow = tau1 ;
        tauhigh = taulim ; 
        Value += RangeIntLin(physicsVector,nbin) ;
        ltaulow = ltaulim ;
        ltauhigh = log(tau) ;
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
    tau = LowEdgeEnergy/Mass;
    ltaulow = log(tauold);
    ltauhigh = log(tau);
    Value = oldValue+RangeIntLog(physicsVector,nbin);
    rangeVector->PutValue(j,Value);
    oldValue = Value ;
    tauold = tau ;
  }
}    

void G4VIhEnergyLoss::BuildLabTimeVector(G4int materialIndex,
                                     G4PhysicsLogVector* timeVector)
//  create lab time vector for a material
{

  G4int nbin=100;
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
    Value = oldValue+LabTimeIntLog(physicsVector,nbin);
    timeVector->PutValue(j,Value);
    oldValue = Value ;
    tauold = tau ;
  }
}

void G4VIhEnergyLoss::BuildProperTimeVector(G4int materialIndex,
                                     G4PhysicsLogVector* timeVector)
//  create proper time vector for a material
{
  G4int nbin=100;
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
    Value = oldValue+ProperTimeIntLog(physicsVector,nbin);
    timeVector->PutValue(j,Value);
    oldValue = Value ;
    tauold = tau ;
  }
}

G4double G4VIhEnergyLoss::RangeIntLin(G4PhysicsVector* physicsVector,
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


G4double G4VIhEnergyLoss::RangeIntLog(G4PhysicsVector* physicsVector,
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

G4double G4VIhEnergyLoss::LabTimeIntLog(G4PhysicsVector* physicsVector,
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

G4double G4VIhEnergyLoss::ProperTimeIntLog(G4PhysicsVector* physicsVector,
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

void G4VIhEnergyLoss::BuildRangeCoeffATable(
                            const G4ParticleDefinition& aParticleType)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "A"
{
  const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();
  G4double Charge = aParticleType.GetPDGCharge() ;

  if(Charge>0.)
  {
    if(thepRangeCoeffATable)
    { thepRangeCoeffATable->clearAndDestroy();
      delete thepRangeCoeffATable; }
    thepRangeCoeffATable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffATable = thepRangeCoeffATable ;
    theRangeTable = theRangepTable ;
  }
  else  
  {
    if(thepbarRangeCoeffATable)
    { thepbarRangeCoeffATable->clearAndDestroy();
      delete thepbarRangeCoeffATable; }
    thepbarRangeCoeffATable = new G4PhysicsTable(numOfMaterials);
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
  for (G4int J=0; J<numOfMaterials; J++)
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
      Value = (w1*Rip + w2*Ri + w3*Rim)/(Ti*Ti) ; 

      aVector->PutValue(i,Value);
      Ti = RTable*Ti ;
    }
  
    theRangeCoeffATable->insert(aVector);
  } 
}


void G4VIhEnergyLoss::BuildRangeCoeffBTable(
                            const G4ParticleDefinition& aParticleType)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "B"
{
  const G4MaterialTable* theMaterialTable=
                               G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();
  G4double Charge = aParticleType.GetPDGCharge() ;

  if(Charge>0.)
  {
    if(thepRangeCoeffBTable)
    { thepRangeCoeffBTable->clearAndDestroy();
      delete thepRangeCoeffBTable; }
    thepRangeCoeffBTable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffBTable = thepRangeCoeffBTable ;
    theRangeTable = theRangepTable ;
  }
  else
  {
    if(thepbarRangeCoeffBTable)
    { thepbarRangeCoeffBTable->clearAndDestroy();
      delete thepbarRangeCoeffBTable; }
    thepbarRangeCoeffBTable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffBTable = thepbarRangeCoeffBTable ;
    theRangeTable = theRangepbarTable ;
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
      Value = (w1*Rip + w2*Ri + w3*Rim)/Ti;

      aVector->PutValue(i,Value);
      Ti = RTable*Ti ;
    }
    theRangeCoeffBTable->insert(aVector);
  } 
}

void G4VIhEnergyLoss::BuildRangeCoeffCTable(
                            const G4ParticleDefinition& aParticleType)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "C"
{
  const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();
  G4double Charge = aParticleType.GetPDGCharge() ;

  if(Charge>0.)
  {
    if(thepRangeCoeffCTable)
    { thepRangeCoeffCTable->clearAndDestroy();
      delete thepRangeCoeffCTable; }
    thepRangeCoeffCTable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffCTable = thepRangeCoeffCTable ;
    theRangeTable = theRangepTable ;
  }
  else
  {
    if(thepbarRangeCoeffCTable)
    { thepbarRangeCoeffCTable->clearAndDestroy();
      delete thepbarRangeCoeffCTable; }
    thepbarRangeCoeffCTable = new G4PhysicsTable(numOfMaterials);
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
  for (G4int J=0; J<numOfMaterials; J++)
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

void G4VIhEnergyLoss::BuildInverseRangeTable(
                             const G4ParticleDefinition& aParticleType)
// Build inverse table of the range table
{
  G4double SmallestRange,BiggestRange ;
  G4bool isOut ;
  const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();
  if(&aParticleType == G4Proton::Proton())
  {
    if(theInverseRangepTable)
    { theInverseRangepTable->clearAndDestroy();
      delete theInverseRangepTable; }
    theInverseRangepTable = new G4PhysicsTable(numOfMaterials);
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
    theInverseRangepbarTable = new G4PhysicsTable(numOfMaterials);
    theInverseRangeTable = theInverseRangepbarTable ;
    theRangeTable = theRangepbarTable ;
    theDEDXTable =  theDEDXpbarTable ;
    theRangeCoeffATable = thepbarRangeCoeffATable ;
    theRangeCoeffBTable = thepbarRangeCoeffBTable ;
    theRangeCoeffCTable = thepbarRangeCoeffCTable ;
  }

  // loop for materials 
  for (G4int J=0;  J<numOfMaterials; J++)
  {
    SmallestRange = (*theRangeTable)(J)->
                       GetValue(LowestKineticEnergy,isOut) ;
    BiggestRange = (*theRangeTable)(J)->
                       GetValue(HighestKineticEnergy,isOut) ;
    G4PhysicsLogVector* aVector;
    aVector = new G4PhysicsLogVector(SmallestRange,
                            BiggestRange,TotBin);

    InvertRangeVector(J, aVector);

    theInverseRangeTable->insert(aVector);
  }
}

void G4VIhEnergyLoss::InvertRangeVector(G4int materialIndex,
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

G4double G4VIhEnergyLoss::GetConstraints(const G4DynamicParticle *aParticle,
                                              G4Material *aMaterial)
{
  // returns the Step limit = range !
  // dRoverRange is the max. allowed relative range loss in one step
  // it calculates dEdx and the range as well....

  G4double KineticEnergy,StepLimit;

  G4double Charge = aParticle->GetDefinition()->GetPDGCharge() ;

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

  // compute the  Step limit ..............

  StepLimit = fRangeNow ;

  return StepLimit ;
}

G4VParticleChange* G4VIhEnergyLoss::AlongStepDoIt( 
                              const G4Track& trackData,const G4Step& stepData) 
 // compute the energy loss after a step 
{
  const G4DynamicParticle* aParticle;
  G4Material* aMaterial;
  G4bool isOut;
  G4double E,finalT,Step,Charge,ChargeSquare,MeanLoss ;
  G4int index ;

  // do not track further if kin.energy < 1. eV
  const G4double MinKineticEnergy = 1.*eV ;
 
  const G4double stepoverrangelimit = 0.02 ;

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;

  // get the actual (true) Step length from stepData 
  Step = stepData.GetStepLength() ;

  // there is no loss for Step=0. !
  if( Step == 0.)
    return &aParticleChange ;  

  aParticle = trackData.GetDynamicParticle() ;
  Charge = aParticle->GetDefinition()->GetPDGCharge() ;
  ChargeSquare = Charge*Charge ;

  index = aMaterial->GetIndex() ;
  E = aParticle->GetKineticEnergy() ;

  if(E < MinKineticEnergy) 
  {
    finalT = 0.;
    MeanLoss = E ;
  }
  else
  {
    if(Step >= fRangeNow )
    {
      finalT = 0.0;
      MeanLoss = E ;
    }
    else if( E > HighestKineticEnergy) 
    {
      MeanLoss = Step*fdEdx ;
      if(MeanLoss > E) MeanLoss = E ;

      finalT = E - MeanLoss ;
    }
    else
    {
      if(Step>stepoverrangelimit*fRangeNow)
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
      else
      {
        MeanLoss = Step*fdEdx ;
      }
      finalT = E - MeanLoss ;
    }
  } 

  //  now the loss with fluctuation
  if ((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowestKineticEnergy))
  {
    MeanLoss /= ChargeSquare ;
    finalT = E-GetLossWithFluct(aParticle,aMaterial,ChargeSquare,MeanLoss,Step)*ChargeSquare ;
    if (finalT < 0.) finalT = 0. ;
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

  aParticleChange.SetNumberOfSecondaries(0);
  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(E-finalT) ;

  return &aParticleChange ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VIhEnergyLoss::GetLossWithFluct(const G4DynamicParticle* aParticle,
				 		G4Material* aMaterial,
				 		G4double ChargeSquare,
				 		G4double	   MeanLoss,
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
		na			= G4RandGauss::shoot(namean,sa);
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







