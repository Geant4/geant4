// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hLowEnergyLoss.cc,v 1.8 2000-11-05 12:17:06 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      GEANT4 Collaboration
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
// 31/03/00 : rename to lowenergy as G4hLowEnergyLoss.cc V.Ivanchenko
// 05/11/00 : new method to calculate particle ranges
// --------------------------------------------------------------
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hLowEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4Poisson.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Initialisation of static members ******************************************
// contributing processes : ion.loss ->NumberOfProcesses is initialized
//   to 1 . YOU DO NOT HAVE TO CHANGE this variable for a 'normal' run.
// You have to change NumberOfProcesses  
// if you invent a new process contributing to the cont. energy loss,
//   NumberOfProcesses should be 2 in this case,
//  or for debugging purposes.
//  The NumberOfProcesses data member can be changed using the (public static)
//  functions Get/Set/Plus/MinusNumberOfProcesses (see G4hLowEnergyLoss.hh)

G4int G4hLowEnergyLoss::NumberOfProcesses = 1 ;

G4int            G4hLowEnergyLoss::CounterOfProcess = 0 ;
G4PhysicsTable** G4hLowEnergyLoss::RecorderOfProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4hLowEnergyLoss::CounterOfpProcess = 0 ;
G4PhysicsTable** G4hLowEnergyLoss::RecorderOfpProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4hLowEnergyLoss::CounterOfpbarProcess = 0 ;
G4PhysicsTable** G4hLowEnergyLoss::RecorderOfpbarProcess =
                                           new G4PhysicsTable*[10] ;

G4PhysicsTable* G4hLowEnergyLoss::theDEDXpTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theDEDXpbarTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theRangepTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theRangepbarTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theInverseRangepTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theInverseRangepbarTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theLabTimepTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theLabTimepbarTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theProperTimepTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theProperTimepbarTable = NULL ;

G4PhysicsTable* G4hLowEnergyLoss::thepRangeCoeffATable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::thepRangeCoeffBTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::thepRangeCoeffCTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::thepbarRangeCoeffATable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::thepbarRangeCoeffBTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::thepbarRangeCoeffCTable = NULL ;

G4PhysicsTable* G4hLowEnergyLoss::theDEDXTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theRangeTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theInverseRangeTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theLabTimeTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theProperTimeTable = NULL ;

G4PhysicsTable* G4hLowEnergyLoss::theRangeCoeffATable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theRangeCoeffBTable = NULL ;
G4PhysicsTable* G4hLowEnergyLoss::theRangeCoeffCTable = NULL ;

//const G4Proton* G4hLowEnergyLoss::theProton=G4Proton::Proton() ;
//const G4AntiProton* G4hLowEnergyLoss::theAntiProton=G4AntiProton::AntiProton() ;

G4double G4hLowEnergyLoss::ParticleMass;
G4double G4hLowEnergyLoss::ptableElectronCutInRange = 0.0*mm ;
G4double G4hLowEnergyLoss::pbartableElectronCutInRange = 0.0*mm ;

G4double G4hLowEnergyLoss::Mass,
         G4hLowEnergyLoss::taulow, 
         G4hLowEnergyLoss::tauhigh, 
         G4hLowEnergyLoss::ltaulow, 
         G4hLowEnergyLoss::ltauhigh; 

G4double G4hLowEnergyLoss::dRoverRange = 0.20 ;
G4double G4hLowEnergyLoss::finalRange = 200.*micrometer ;

G4double     G4hLowEnergyLoss::c1lim = dRoverRange ;
G4double     G4hLowEnergyLoss::c2lim = 2.*(1.-dRoverRange)*finalRange ;
G4double     G4hLowEnergyLoss::c3lim = -(1.-dRoverRange)*finalRange*finalRange;

G4double         G4hLowEnergyLoss::Charge ;   

G4bool   G4hLowEnergyLoss::rndmStepFlag   = false ;
G4bool   G4hLowEnergyLoss::EnlossFlucFlag = true ;

G4double G4hLowEnergyLoss::LowestKineticEnergy,G4hLowEnergyLoss::HighestKineticEnergy;	
G4int    G4hLowEnergyLoss::TotBin ;
G4double G4hLowEnergyLoss::RTable,G4hLowEnergyLoss::LOGRTable;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// constructor and destructor
 
G4hLowEnergyLoss::G4hLowEnergyLoss(const G4String& processName)
   : G4VContinuousDiscreteProcess (processName),
     theLossTable (NULL),
     MinKineticEnergy(0.0), 
     linLossLimit(0.05),
     lastMaterial (NULL),
     MaxExcitationNumber (1.e6),
     probLimFluct (0.01),
     nmaxDirectFluct (100),
     nmaxCont1(4),
     nmaxCont2(16) 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hLowEnergyLoss::~G4hLowEnergyLoss() 
{
     if(theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable;
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
void G4hLowEnergyLoss::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{
  //  calculate data members TotBin,LOGRTable,RTable first
 
  G4double binning = dRoverRange;
  G4double lrate = log(HighestKineticEnergy/LowestKineticEnergy);
  G4int    nbin =  G4int(lrate/log(1.+binning) + 0.5 );
  nbin = (nbin+25)/50;
  LOGRTable=lrate/TotBin;
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

      if(CounterOfProcess == NumberOfProcesses)
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

      if(CounterOfProcess == NumberOfProcesses)
      {
        if(theDEDXpbarTable)
        { theDEDXpbarTable->clearAndDestroy();
          delete theDEDXpbarTable; }
        theDEDXpbarTable = new G4PhysicsTable(numOfMaterials);
        theDEDXTable = theDEDXpbarTable;
        pbartableElectronCutInRange = ElectronCutInRange ;
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
    proton_mass_c2/aParticleType.GetPDGMass(),TotBin);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
      
void G4hLowEnergyLoss::BuildRangeTable(
                             const G4ParticleDefinition& aParticleType)
// Build range table from the energy loss table
{
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyLoss::BuildTimeTables(
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyLoss::BuildRangeVector(G4int materialIndex,
                                        G4PhysicsLogVector* rangeVector)
//  create range vector for a material
{
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyLoss::BuildLabTimeVector(G4int materialIndex,
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyLoss::BuildProperTimeVector(G4int materialIndex,
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyLoss::RangeIntLin(G4PhysicsVector* physicsVector,
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyLoss::RangeIntLog(G4PhysicsVector* physicsVector,
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyLoss::LabTimeIntLog(G4PhysicsVector* physicsVector,
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyLoss::ProperTimeIntLog(G4PhysicsVector* physicsVector,
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyLoss::BuildRangeCoeffATable(
                            const G4ParticleDefinition& aParticleType)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "A"
{
  const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyLoss::BuildRangeCoeffBTable(
                            const G4ParticleDefinition& aParticleType)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "B"
{
  const G4MaterialTable* theMaterialTable=
                               G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyLoss::BuildRangeCoeffCTable(
                            const G4ParticleDefinition& aParticleType)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "C"
{
  const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyLoss::BuildInverseRangeTable(
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyLoss::InvertRangeVector(G4int materialIndex,
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



