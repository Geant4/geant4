// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuEnergyLoss.cc,v 1.6 1999-06-18 11:30:47 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4MuEnergyLoss physics process -----------
//                by Laszlo Urban, September 1997
// **************************************************************
// It is the implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of muons.
// **************************************************************
//
// corrections by L.Urban on 27/05/98 (other corrs come soon!)
// cleanup  L.Urban on 23/10/98
// --------------------------------------------------------------
 

#include "G4MuEnergyLoss.hh"
#include "G4EnergyLossTables.hh"

// Initialisation of static members *******************************************
// contributing processes : ion.loss,bremsstrahlung,pair production
//          ->NbOfProcesses is initialized to 3.
//  YOU DO NOT HAVE TO CHANGE this variable for a 'normal' run.
// You have to change NbOfProcesses  
// if you invent a new process contributing to the cont. energy loss,
//   NbOfProcesses should be 4 in this case,
//  or for debugging purposes.
//  The NbOfProcesses data member can be changed using the (public static)
//  functions Get/Set/Plus/MinusNbOfProcesses (see G4MuEnergyLoss.hh)

G4int G4MuEnergyLoss::NbOfProcesses = 3 ; // !!!!!!!!!!!!!!!

G4PhysicsTable** G4MuEnergyLoss::RecorderOfmuplusProcess  =
                                           new G4PhysicsTable*[10] ;
G4PhysicsTable** G4MuEnergyLoss::RecorderOfmuminusProcess =
                                           new G4PhysicsTable*[10] ;
G4int       G4MuEnergyLoss::CounterOfmuplusProcess  = 0 ;
G4int       G4MuEnergyLoss::CounterOfmuminusProcess = 0 ;

G4bool           G4MuEnergyLoss::rndmStepFlag   = false;
G4bool           G4MuEnergyLoss::EnlossFlucFlag = true;
G4double         G4MuEnergyLoss::dRoverRange    = 20*perCent;
G4double         G4MuEnergyLoss::finalRange     = 200*micrometer;

G4PhysicsTable* G4MuEnergyLoss::theDEDXmuplusTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::theRangemuplusTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::theInverseRangemuplusTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::theLabTimemuplusTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::theProperTimemuplusTable = NULL ;

G4PhysicsTable* G4MuEnergyLoss::theDEDXmuminusTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::theRangemuminusTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::theInverseRangemuminusTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::theLabTimemuminusTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::theProperTimemuminusTable = NULL ;

G4PhysicsTable* G4MuEnergyLoss::themuplusRangeCoeffATable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::themuplusRangeCoeffBTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::themuplusRangeCoeffCTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::themuminusRangeCoeffATable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::themuminusRangeCoeffBTable = NULL ;
G4PhysicsTable* G4MuEnergyLoss::themuminusRangeCoeffCTable = NULL ;
 
G4EnergyLossMessenger* G4MuEnergyLoss::eLossMessenger = NULL ;

// constructor and destructor
 
G4MuEnergyLoss::G4MuEnergyLoss(const G4String& processName)
   : G4VContinuousDiscreteProcess (processName),
     LowestKineticEnergy(1.00*keV),
     HighestKineticEnergy(1000000.*TeV),
     MaxExcitationNumber (1.e6),
     probLimFluct (0.01),
     nmaxDirectFluct (100),
     nmaxCont1(4),
     nmaxCont2(16),
     theLossTable(NULL),
     theRangeCoeffATable(NULL),
     theRangeCoeffBTable(NULL),
     theRangeCoeffCTable(NULL),
     lastMaterial(NULL),
     lastgammaCutInRange(0.),
     lastelectronCutInRange(0.),
     theElectron ( G4Electron::Electron() ),
     thePositron ( G4Positron::Positron() ),
     theMuonPlus ( G4MuonPlus::MuonPlus() ),
     theMuonMinus ( G4MuonMinus::MuonMinus() )
{ }

G4MuEnergyLoss::~G4MuEnergyLoss() 
{
     if(theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable;
     }

}
 
   void G4MuEnergyLoss::BuildDEDXTable(
                         const G4ParticleDefinition& aParticleType)
{
  //  calculate data members TotBin,LOGRTable,RTable first
  G4double binning = dRoverRange;
  G4double lrate = log(HighestKineticEnergy/LowestKineticEnergy);
  G4int    nbin =  G4int(lrate/log(1.+binning) + 0.5 );
  nbin = (nbin+25)/50;
  TotBin =50*nbin ;
  if (TotBin<50) TotBin = 50;
  if (TotBin>500) TotBin = 500;
  LOGRTable=lrate/TotBin;
  RTable   =exp(LOGRTable);

  G4bool MakeTable ;
  ParticleMass = aParticleType.GetPDGMass() ; 
  G4double Charge = aParticleType.GetPDGCharge() ;
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
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);   
        // loop for the kinetic energy
        for (G4int i=0; i<TotBin; i++)
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

    lastgammaCutInRange = gammaCutInRange ;
    lastelectronCutInRange = electronCutInRange ;
  }
}
      
  void G4MuEnergyLoss::BuildRangeTable(
                             const G4ParticleDefinition& aParticleType)
// Build range table from the energy loss table
{
   G4double Charge = aParticleType.GetPDGCharge() ;
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

   for (G4int J=0;  J<numOfMaterials; J++)
   {
     G4PhysicsLogVector* aVector;

     aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                            HighestKineticEnergy,TotBin);

     BuildRangeVector(J, aVector);

     theRangeTable->insert(aVector);

   }

}    

  void G4MuEnergyLoss::BuildTimeTables(
                             const G4ParticleDefinition& aParticleType)
// Build time tables from the energy loss table
{
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

void G4MuEnergyLoss::BuildRangeVector(G4int materialIndex,
                                     G4PhysicsLogVector* rangeVector)
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

void G4MuEnergyLoss::BuildLabTimeVector(G4int materialIndex,
                                     G4PhysicsLogVector* timeVector)
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

void G4MuEnergyLoss::BuildProperTimeVector(G4int materialIndex,
                                     G4PhysicsLogVector* timeVector)
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


G4double G4MuEnergyLoss::RangeIntLin(G4PhysicsVector* physicsVector,
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


G4double G4MuEnergyLoss::RangeIntLog(G4PhysicsVector* physicsVector,
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

G4double G4MuEnergyLoss::LabTimeIntLog(G4PhysicsVector* physicsVector,
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

G4double G4MuEnergyLoss::ProperTimeIntLog(G4PhysicsVector* physicsVector,
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

void G4MuEnergyLoss::BuildRangeCoeffATable(
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
   G4PhysicsLinearVector* aVector=new G4PhysicsLinearVector(0.,binmax, TotBin);

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


void G4MuEnergyLoss::BuildRangeCoeffBTable(
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
   G4int binmax=TotBin ;

   G4PhysicsLinearVector* aVector=new G4PhysicsLinearVector(0.,binmax, TotBin);

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
void G4MuEnergyLoss::BuildRangeCoeffCTable(
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
   G4PhysicsLinearVector* aVector=new G4PhysicsLinearVector(0.,binmax, TotBin);

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

  void G4MuEnergyLoss::BuildInverseRangeTable(
                             const G4ParticleDefinition& aParticleType)
// Build inverse table of the range table
{
  G4double SmallestRange,BiggestRange ;
  G4bool isOut ;

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

void G4MuEnergyLoss::InvertRangeVector(G4int materialIndex,
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

G4double G4MuEnergyLoss::GetConstraints(const G4DynamicParticle *aParticle,
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


  G4double c1=dRoverRange, c2=2.*(1.-dRoverRange)*finalRange ,
                 c3=-(1.-dRoverRange)*finalRange*finalRange ;

  G4double Thigh = HighestKineticEnergy/RTable ;
  KineticEnergy = aParticle->GetKineticEnergy();

  bin = G4int(log(KineticEnergy/LowestKineticEnergy)/LOGRTable) ;
  EnergyBinNumber = bin ;

  index = aMaterial->GetIndex() ;

  if( KineticEnergy < LowestKineticEnergy )
    {
      fdEdx = sqrt(KineticEnergy/LowestKineticEnergy)*
             (*theDEDXTable)(index)->GetValue(LowestKineticEnergy,isOutRange) ;

      fRangeNow = sqrt(KineticEnergy/LowestKineticEnergy)*
             (*theRangeTable)(index)->GetValue(LowestKineticEnergy,isOutRange) ;

      StepLimit = fRangeNow ;
    }
  else if (KineticEnergy > Thigh )
    {
      fdEdx = (*theDEDXTable)(index)->GetValue(Thigh,isOutRange);
      fRangeNow = (*theRangeTable)(index)->GetValue(Thigh,isOutRange);
      if (fdEdx > 0.) fRangeNow += (KineticEnergy-Thigh)/fdEdx;
      StepLimit = c1*fRangeNow;
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
         StepLimit = c1*fRangeNow+c2+c3/fRangeNow ;

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

G4VParticleChange* G4MuEnergyLoss::AlongStepDoIt( 
                              const G4Track& trackData,const G4Step& stepData) 
{
 // compute the energy loss after a Step

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
   const G4double linLossLimit = 0.02 ;
  
  G4double MeanLoss, finalT;
 
  if (E < MinKineticEnergy)    finalT = 0.; 
  else if ( E<= LowestKineticEnergy)
  {
    if (Step >= fRangeNow)  finalT = 0.;
    else finalT = E - Step*fdEdx ;
  }
   
  else if (E>=HighestKineticEnergy) finalT = E - Step*fdEdx;

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
  if ((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowestKineticEnergy))

  {
    finalT = E-GetLossWithFluct(aParticle,aMaterial,MeanLoss);
    if (finalT < 0.) finalT = E-MeanLoss;
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

G4double G4MuEnergyLoss::GetLossWithFluct(const G4DynamicParticle* aParticle,
                                               G4Material* aMaterial,
                                               G4double    MeanLoss)
//  calculate actual loss from the mean loss
//  The model used to get the fluctuation is the same as in Glandz in Geant3.
{
  static const G4double Tlow=10.*keV ;

  // check if the material has changed ( cache mechanism)

  if (aMaterial != lastMaterial)
    {
      lastMaterial = aMaterial;
      imat         = aMaterial->GetIndex();
      f1Fluct      = aMaterial->GetIonisation()->GetF1fluct();
      f2Fluct      = aMaterial->GetIonisation()->GetF2fluct();
      e1Fluct      = aMaterial->GetIonisation()->GetEnergy1fluct();
      e2Fluct      = aMaterial->GetIonisation()->GetEnergy2fluct();
      e1LogFluct   = aMaterial->GetIonisation()->GetLogEnergy1fluct();
      e2LogFluct   = aMaterial->GetIonisation()->GetLogEnergy2fluct();
      rateFluct    = aMaterial->GetIonisation()->GetRateionexcfluct();
      ipotFluct    = aMaterial->GetIonisation()->GetMeanExcitationEnergy();
      ipotLogFluct = aMaterial->GetIonisation()->GetLogMeanExcEnergy();
    }

  G4double threshold,w1,w2,w3,lnw3,C,prob,
           beta2,suma,e0,Em,loss,lossc ,w;
  G4double a1,a2,a3;
  long p1,p2,p3;
  G4int nb;
  G4double Corrfac, na,alfa,rfac,namean,sa,alfa1,ea,sea;
  G4double dp1,dnmaxDirectFluct,dp3,dnmaxCont2;

  // get particle data
  G4double Tkin   = aParticle->GetKineticEnergy();
  G4double charge = aParticle->GetDefinition()->GetPDGCharge();
  threshold =((*G4Electron::Electron()).GetCutsInEnergy())[imat];

  G4double rmass = electron_mass_c2/ParticleMass;
  G4double tau   = Tkin/ParticleMass, tau1 = tau+1., tau2 = tau*(tau+2.);
  G4double Tm    = 2.*electron_mass_c2*tau2/(1.+2.*tau1*rmass+rmass*rmass)
                  -ipotFluct;
  if (Tm < 0.) Tm = 0.;
  else if (Tm > threshold) Tm = threshold;

  w1 = Tm+ipotFluct;
  w2 = w1/ipotFluct;
  w3 = 2.*electron_mass_c2*tau2;
  lnw3 = log(w3);
  beta2 = tau2/(tau1*tau1);

  C = (1.-rateFluct)*MeanLoss/(lnw3-ipotLogFluct-beta2);

  a1 = C*f1Fluct*(lnw3-e1LogFluct-beta2)/e1Fluct;
  a2 = C*f2Fluct*(lnw3-e2LogFluct-beta2)/e2Fluct;
  if (Tm > 0.) a3 = rateFluct*MeanLoss*Tm/(ipotFluct*w1*log(w2));
  else { a1 /= rateFluct; a2 /= rateFluct; a3 = 0.;}
  suma = a1+a2+a3;
 
  //no fluctuation if the loss is too big
  if (suma > MaxExcitationNumber) return MeanLoss;

  suma<50.? prob = exp(-suma) : prob = 0.;

  if (prob > probLimFluct)         // very small Step
    {
      e0 = aMaterial->GetIonisation()->GetEnergy0fluct();
      if (Tm <= 0.)
        {
          a1 = MeanLoss/e0;
          p1 = RandPoisson::shoot(a1);
          loss = p1*e0 ;
        }
     else
        {
          Em = Tm+e0;
          a1 = MeanLoss*(Em-e0)/(Em*e0*log(Em/e0));
          p1 = RandPoisson::shoot(a1);
          w  = (Em-e0)/Em;
          // just to save time
          if (p1 > nmaxDirectFluct)
            {
              dp1 = p1;
              dnmaxDirectFluct=nmaxDirectFluct;
              Corrfac = dp1/dnmaxDirectFluct;
              p1 = nmaxDirectFluct;
            }
          else Corrfac = 1.;

          loss = 0.;
          for (long i=0; i<p1; i++) loss += 1./(1.-w*G4UniformRand());
          loss *= (e0*Corrfac);
        }
    }

  else                              // not so small Step
    {
      p1 = RandPoisson::shoot(a1);
      p2 = RandPoisson::shoot(a2);
      loss = p1*e1Fluct+p2*e2Fluct;
      if (loss>0.) loss += (1.-2.*G4UniformRand())*e1Fluct;
      p3 = RandPoisson::shoot(a3);

      lossc = 0.; na = 0.; alfa = 1.;
      if (p3 > nmaxCont2)
        {
          dp3        = p3;
          dnmaxCont2 = nmaxCont2;
          rfac       = dp3/(dnmaxCont2+dp3);
          namean     = p3*rfac;
          sa         = nmaxCont1*rfac;
          na         = RandGauss::shoot(namean,sa);
          if (na > 0.)
            {
              alfa   = w2*(nmaxCont2+p3)/(w2*nmaxCont2+p3);
              alfa1  = alfa*log(alfa)/(alfa-1.);
              ea     = na*ipotFluct*alfa1;
              sea    = ipotFluct*sqrt(na*(alfa-alfa1*alfa1));
              lossc += RandGauss::shoot(ea,sea);
            }
        }

      nb = G4int(p3-na);
      if (nb > 0)
        {
          w2 = alfa*ipotFluct;
          w  = (w1-w2)/w1;
          for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());
        }

      loss += lossc; 
    }

  return loss ;
}


