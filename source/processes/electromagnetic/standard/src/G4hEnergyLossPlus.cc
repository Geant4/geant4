// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hEnergyLossPlus.cc,v 1.8 1999-07-30 10:16:05 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// -----------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hEnergyLossPlus physics process -----------
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
// 01/03/99 : creation of sub-cutoff delta rays, L.Urban
// 28/04/99 : bug fixed in DoIt , L.Urban
// --------------------------------------------------------------

#include "G4hEnergyLossPlus.hh"
#include "G4EnergyLossTables.hh"

// Initialisation of static members ******************************************
// contributing processes : ion.loss ->NumberOfProcesses is initialized
//   to 1 . YOU DO NOT HAVE TO CHANGE this variable for a 'normal' run.
// You have to change NumberOfProcesses  
// if you invent a new process contributing to the cont. energy loss,
//   NumberOfProcesses should be 2 in this case,
//  or for debugging purposes.
//  The NumberOfProcesses data member can be changed using the (public static)
//  functions Get/Set/Plus/MinusNumberOfProcesses (see G4hEnergyLossPlus.hh)

G4int G4hEnergyLossPlus::NumberOfProcesses = 1 ;

G4int            G4hEnergyLossPlus::CounterOfProcess = 0 ;
G4PhysicsTable** G4hEnergyLossPlus::RecorderOfProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4hEnergyLossPlus::CounterOfpProcess = 0 ;
G4PhysicsTable** G4hEnergyLossPlus::RecorderOfpProcess =
                                           new G4PhysicsTable*[10] ;

G4int            G4hEnergyLossPlus::CounterOfpbarProcess = 0 ;
G4PhysicsTable** G4hEnergyLossPlus::RecorderOfpbarProcess =
                                           new G4PhysicsTable*[10] ;

G4PhysicsTable* G4hEnergyLossPlus::theDEDXpTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theDEDXpbarTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theRangepTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theRangepbarTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theInverseRangepTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theInverseRangepbarTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theLabTimepTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theLabTimepbarTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theProperTimepTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theProperTimepbarTable = NULL ;

G4PhysicsTable* G4hEnergyLossPlus::thepRangeCoeffATable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::thepRangeCoeffBTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::thepRangeCoeffCTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::thepbarRangeCoeffATable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::thepbarRangeCoeffBTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::thepbarRangeCoeffCTable = NULL ;

G4PhysicsTable* G4hEnergyLossPlus::theDEDXTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theRangeTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theInverseRangeTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theLabTimeTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theProperTimeTable = NULL ;

G4PhysicsTable* G4hEnergyLossPlus::theRangeCoeffATable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theRangeCoeffBTable = NULL ;
G4PhysicsTable* G4hEnergyLossPlus::theRangeCoeffCTable = NULL ;

const G4Proton* G4hEnergyLossPlus::theProton=G4Proton::Proton() ;
const G4AntiProton* G4hEnergyLossPlus::theAntiProton=G4AntiProton::AntiProton() ;

G4double G4hEnergyLossPlus::ParticleMass;
G4double G4hEnergyLossPlus::ptableElectronCutInRange = 0.0*mm ;
G4double G4hEnergyLossPlus::pbartableElectronCutInRange = 0.0*mm ;

G4double G4hEnergyLossPlus::Mass,
         G4hEnergyLossPlus::taulow, 
         G4hEnergyLossPlus::tauhigh, 
         G4hEnergyLossPlus::ltaulow, 
         G4hEnergyLossPlus::ltauhigh; 

G4double G4hEnergyLossPlus::dRoverRange = 0.20 ;
G4double G4hEnergyLossPlus::finalRange = 200.*micrometer ;

G4double     G4hEnergyLossPlus::c1lim = dRoverRange ;
G4double     G4hEnergyLossPlus::c2lim = 2.*(1.-dRoverRange)*finalRange ;
G4double     G4hEnergyLossPlus::c3lim = -(1.-dRoverRange)*finalRange*finalRange;

G4double         G4hEnergyLossPlus::MinDeltaCutInRange = 0.1*mm ;
G4double*        G4hEnergyLossPlus::MinDeltaEnergy     = NULL   ;

G4double         G4hEnergyLossPlus::Charge ;   

G4bool   G4hEnergyLossPlus::rndmStepFlag   = false ;
G4bool   G4hEnergyLossPlus::EnlossFlucFlag = true ;

G4double G4hEnergyLossPlus::LowestKineticEnergy;
G4double G4hEnergyLossPlus::HighestKineticEnergy;
G4int G4hEnergyLossPlus::TotBin;
G4double G4hEnergyLossPlus::RTable,G4hEnergyLossPlus::LOGRTable;

G4double G4hEnergyLossPlus::c0N       = 9.0e-21*MeV*MeV*mm*mm ;
G4double G4hEnergyLossPlus::c1N       = 25.0e-21*keV*mm*mm    ;
G4double G4hEnergyLossPlus::c2N       = 13.25e-21*keV*mm*mm   ;
G4double G4hEnergyLossPlus::c3N       = 0.500e-21*mm*mm       ;
G4int    G4hEnergyLossPlus::Ndeltamax = 100                   ;

// constructor and destructor
 
G4hEnergyLossPlus::G4hEnergyLossPlus(const G4String& processName)
   : G4VContinuousDiscreteProcess (processName),
     theLossTable (NULL),
     MinKineticEnergy(1.*eV), 
     linLossLimit(0.05),
     lastMaterial (NULL),
     MaxExcitationNumber (1.e6),
     probLimFluct (0.01),
     nmaxDirectFluct (100),
     nmaxCont1(4),
     nmaxCont2(16) 
{ }

G4hEnergyLossPlus::~G4hEnergyLossPlus() 
{
     if(theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable;
     }
     if(MinDeltaEnergy) delete MinDeltaEnergy ;
}
 
void G4hEnergyLossPlus::BuildDEDXTable(
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

  if(aParticleType.GetParticleName()=="proton")
  {
    // create array for the min. delta cuts in kinetic energy
    G4cout << endl;
    G4cout.precision(5) ;
    G4cout << "hIoni+ Minimum Delta cut in range=" << MinDeltaCutInRange/mm
           << "  mm." << endl;
    G4cout << " min. delta energies (keV) " << endl;
    G4cout << "   material         min.delta energy " << endl;
    G4cout << endl;

    if(MinDeltaEnergy) delete MinDeltaEnergy ;
    MinDeltaEnergy = new G4double [numOfMaterials] ;
    G4double Tlowerlimit = 1.*keV ;
    for(G4int mat=0; mat<numOfMaterials; mat++)
    {
      MinDeltaEnergy[mat] = G4EnergyLossTables::GetPreciseEnergyFromRange(
                            G4Electron::Electron(),MinDeltaCutInRange,
                                       (*theMaterialTable)(mat)) ;
      if(MinDeltaEnergy[mat]<Tlowerlimit) MinDeltaEnergy[mat]=Tlowerlimit ;
      G4cout << setw(20) << (*theMaterialTable)(mat)->GetName()
             << setw(15) << MinDeltaEnergy[mat]/keV << endl;
    }
  }
}
      
void G4hEnergyLossPlus::BuildRangeTable(
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

void G4hEnergyLossPlus::BuildTimeTables(
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

void G4hEnergyLossPlus::BuildRangeVector(G4int materialIndex,
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

void G4hEnergyLossPlus::BuildLabTimeVector(G4int materialIndex,
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

void G4hEnergyLossPlus::BuildProperTimeVector(G4int materialIndex,
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

G4double G4hEnergyLossPlus::RangeIntLin(G4PhysicsVector* physicsVector,
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


G4double G4hEnergyLossPlus::RangeIntLog(G4PhysicsVector* physicsVector,
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

G4double G4hEnergyLossPlus::LabTimeIntLog(G4PhysicsVector* physicsVector,
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

G4double G4hEnergyLossPlus::ProperTimeIntLog(G4PhysicsVector* physicsVector,
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

void G4hEnergyLossPlus::BuildRangeCoeffATable(
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


void G4hEnergyLossPlus::BuildRangeCoeffBTable(
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

void G4hEnergyLossPlus::BuildRangeCoeffCTable(
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

void G4hEnergyLossPlus::BuildInverseRangeTable(
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

void G4hEnergyLossPlus::InvertRangeVector(G4int materialIndex,
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

G4double G4hEnergyLossPlus::GetConstraints(const G4DynamicParticle *aParticle,
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

G4VParticleChange* G4hEnergyLossPlus::AlongStepDoIt( 
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

    else if(( E > HighestKineticEnergy)||( E <= LowestKineticEnergy))
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

  //   subcutoff delta ray production start                          
  G4double MinDeltaEnergyNow,Tc,TmintoProduceDelta,w,ww ;
  G4double rcut,T0,presafety,postsafety,
           delta,fragment,Tmax,mass ;
  G4double frperstep,x1,y1,z1,dx,dy,dz,dTime,time0,DeltaTime;
  G4double epsil = MinKineticEnergy/2. ;

  MinDeltaEnergyNow = MinDeltaEnergy[index] ;
  Tc=G4Electron::Electron()->GetCutsInEnergy()[index];
  const G4ParticleDefinition* aParticleType=aParticle->GetDefinition() ;
  mass=aParticleType->GetPDGMass() ;
  w=mass+electron_mass_c2 ;
  ww=2.*mass-MinDeltaEnergyNow ;
  TmintoProduceDelta=0.5*(sqrt(ww*ww+2.*w*w*MinDeltaEnergyNow/
                       electron_mass_c2)-ww) ;

  if((E > TmintoProduceDelta) && (MeanLoss > MinDeltaEnergyNow)
                                   && (finalT > MinKineticEnergy))
  {
    // max. possible delta energy 
    Tmax = 2.*electron_mass_c2*E*(E+2.*mass)/
           (mass*mass+2.*electron_mass_c2*(E+mass)+
            electron_mass_c2*electron_mass_c2) ;

    rcut=G4Electron::Electron()->GetCuts();

    if(Tc > Tmax) Tc=Tmax ;
    // generate subcutoff delta rays only if Tc>MinDeltaEnergyNow!
    if((Tc > MinDeltaEnergyNow) && (Tmax > MinDeltaEnergyNow))
    {
      presafety  = stepData.GetPreStepPoint()->GetSafety() ;
      postsafety = stepData.GetPostStepPoint()->GetSafety() ;

      if((presafety>=rcut)&&(postsafety>=rcut))
      {
        fragment = 0. ;
      } 
      else
      {
        x1=stepData.GetPreStepPoint()->GetPosition().x();
        y1=stepData.GetPreStepPoint()->GetPosition().y();
        z1=stepData.GetPreStepPoint()->GetPosition().z();
        dx=stepData.GetPostStepPoint()->GetPosition().x()-x1 ;
        dy=stepData.GetPostStepPoint()->GetPosition().y()-y1 ;
        dz=stepData.GetPostStepPoint()->GetPosition().z()-z1 ;
        time0=stepData.GetPreStepPoint()->GetGlobalTime();
        dTime=stepData.GetPostStepPoint()->GetGlobalTime()-time0;

        if((presafety<rcut)&&(postsafety<rcut))
        {
          fragment = Step ;
          frperstep=1. ;
        }
        else if(presafety<rcut)
        {
          delta=presafety*Step/(postsafety-presafety) ;
          fragment=rcut*(Step+delta)/postsafety-delta ;
          frperstep=fragment/Step;
        }
        else if(postsafety<rcut)
        {
          delta=postsafety*Step/(presafety-postsafety) ;
          fragment=rcut*(Step+delta)/presafety-delta ;
          x1 += dx;
          y1 += dy;
          z1 += dz;  
          time0 += dTime ;

          frperstep=-fragment/Step;
        }
      }

      if(fragment>0.)
      {
        T0=G4EnergyLossTables::GetPreciseEnergyFromRange(
                                             G4Electron::Electron(),
                                             min(presafety,postsafety),
                                             aMaterial) ;

        // absolute lower limit for T0
        if(T0<MinDeltaEnergyNow) T0=MinDeltaEnergyNow ;

        // compute nb of delta rays to be generated
        G4int N=int(fragment*(c0N/(E*T0)+c1N/T0-(c2N+c3N*T0)/Tc)* 
                (aMaterial->GetTotNbOfElectPerVolume())+0.5) ;

        if(N > 0)
        {
          G4double Tkin,Etot,P,T,p,costheta,sintheta,phi,dirx,diry,dirz,
                   Pnew,Px,Py,Pz,delToverTc,
                   delTkin,delLoss,rate,
                   urandom ;
          G4ThreeVector ParticleDirection ;
          G4StepPoint *point ;
  
          Tkin = E ;
          Etot = Tkin+mass ;
          P    = sqrt(Tkin*(Etot+mass)) ;

          aParticleChange.SetNumberOfSecondaries(N);
          G4int subdelta = 0;
          do {
               subdelta += 1 ;

               Tmax = 2.*electron_mass_c2*Tkin*(Tkin+2.*mass)/
                      (mass*mass+2.*electron_mass_c2*(Tkin+mass)+
                        electron_mass_c2*electron_mass_c2) ;

               if(Tc>Tmax) Tc = Tmax ;

               //check if there is enough energy ....
               if((Tkin>TmintoProduceDelta)&&(Tc > T0)&&(MeanLoss>0.))
               {
                 delToverTc=1.-T0/Tc ;
                 T=T0/(1.-delToverTc*G4UniformRand()) ;
                 if(T > MeanLoss) T=MeanLoss ;
                 MeanLoss -= T ;
                 p=sqrt(T*(T+2.*electron_mass_c2)) ;

                 costheta = T*(Etot+electron_mass_c2)/(P*p) ;
                 if(costheta<-1.) costheta=-1.;
                 if(costheta> 1.) costheta= 1.;

                 phi=twopi*G4UniformRand() ;
                 sintheta=sqrt(1.-costheta*costheta);
                 dirx=sintheta*cos(phi);
                 diry=sintheta*sin(phi);
                 dirz=costheta;
               }
               else
               {
                 T=epsil ;
                 p=sqrt(T*(T+2.*electron_mass_c2)) ;
                 dirx=0.;
                 diry=0.;
                 dirz=1.;
               } 

               urandom = G4UniformRand() ;
               // distribute x,y,z along Pre-Post !
               G4double xd,yd,zd ;
               xd=x1+frperstep*dx*urandom ;
               yd=y1+frperstep*dy*urandom ;
               zd=z1+frperstep*dz*urandom ;
               G4ThreeVector DeltaPosition(xd,yd,zd) ;
               DeltaTime=time0+frperstep*dTime*urandom ;
               ParticleDirection=stepData.GetPostStepPoint()->
                                   GetMomentumDirection() ;

               G4ThreeVector DeltaDirection(dirx,diry,dirz) ;
               DeltaDirection.rotateUz(ParticleDirection);

               G4DynamicParticle* theDelta = new G4DynamicParticle ;
               theDelta->SetDefinition(G4Electron::Electron());
               theDelta->SetKineticEnergy(T);

               theDelta->SetMomentumDirection(DeltaDirection.x(),
                              DeltaDirection.y(),DeltaDirection.z());

               // update initial particle,fill ParticleChange
               Tkin -= T ;
               Px =(P*ParticleDirection.x()-p*DeltaDirection.x()) ;
               Py =(P*ParticleDirection.y()-p*DeltaDirection.y()) ;
               Pz =(P*ParticleDirection.z()-p*DeltaDirection.z()) ;
               Pnew = sqrt(Px*Px+Py*Py+Pz*Pz) ;
               Px /= Pnew ;
               Py /= Pnew ;
               Pz /= Pnew ;
               P  = Pnew ;
               G4ThreeVector ParticleDirectionnew(Px,Py,Pz) ;
               ParticleDirection = ParticleDirectionnew;

               G4Track* deltaTrack =
                        new G4Track(theDelta,DeltaTime,DeltaPosition);
               deltaTrack->
                SetTouchable(stepData.GetPostStepPoint()->GetTouchable()) ;

               deltaTrack->SetParentID(trackData.GetTrackID()) ;

               aParticleChange.AddSecondary(deltaTrack) ;

             } while (subdelta<N) ;

             // update the particle direction and kinetic energy
             aParticleChange.SetMomentumChange(Px,Py,Pz) ;
             E = Tkin ;
           }

         }
       }
     }

     finalT = E - MeanLoss ;
     if(finalT < MinKineticEnergy) finalT = 0. ;

  //  now the loss with fluctuation
  if((EnlossFlucFlag) && (finalT > 0.) && (finalT < E)&&(E > LowestKineticEnergy))
  {
    MeanLoss /= ChargeSquare ;
    finalT = E-GetLossWithFluct(aParticle,aMaterial,MeanLoss)*ChargeSquare ;
    if (finalT < 0.) finalT = E-MeanLoss ;
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

  // aParticleChange.SetNumberOfSecondaries(0);
  aParticleChange.SetEnergyChange( finalT ) ;
  aParticleChange.SetLocalEnergyDeposit(E-finalT) ;

  return &aParticleChange ;
}

G4double G4hEnergyLossPlus::GetLossWithFluct(const G4DynamicParticle* aParticle,
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
  G4double siga ;
  static const G4double alim=10.;

  // get particle data
  G4double Tkin   = aParticle->GetKineticEnergy();
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
  if (suma > MaxExcitationNumber)  return  MeanLoss;

  suma<50.? prob = exp(-suma) : prob = 0.;

  if (prob > probLimFluct)         // very small Step
    {
      e0 = aMaterial->GetIonisation()->GetEnergy0fluct();
      if (Tm <= 0.)
        {
          a1 = MeanLoss/e0;
          if(a1>alim)
          {
            siga=sqrt(a1) ;
            p1 = max(0,int(RandGauss::shoot(a1,siga)+0.5));
          }
          else
            p1 = RandPoisson::shoot(a1);
          loss = p1*e0 ;
        }
     else
        {
          Em = Tm+e0;
          a1 = MeanLoss*(Em-e0)/(Em*e0*log(Em/e0));
          if(a1>alim)
          {
            siga=sqrt(a1) ;
            p1 = max(0,int(RandGauss::shoot(a1,siga)+0.5));
          }
          else
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
      if(a1>alim)
      {
        siga=sqrt(a1) ;
        p1 = max(0,int(RandGauss::shoot(a1,siga)+0.5));
      }
      else
       p1 = RandPoisson::shoot(a1);
      if(a2>alim)
      {
        siga=sqrt(a2) ;
        p2 = max(0,int(RandGauss::shoot(a2,siga)+0.5));
      }
      else
        p2 = RandPoisson::shoot(a2);
      loss = p1*e1Fluct+p2*e2Fluct;
      if (loss>0.) loss += (1.-2.*G4UniformRand())*e1Fluct;
      if(a3>alim)
      {
        siga=sqrt(a3) ;
        p3 = max(0,int(RandGauss::shoot(a3,siga)+0.5));
      }
      else
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


