// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIenergyLoss.cc,v 1.2 1999-04-16 09:02:01 grichine Exp $
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
//      ---------- G4PAIenergyLoss physics process -----------
//                by V. Grichine, 30 Nov 1997 
// **************************************************************
// It is the first implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the energy loss of charged hadrons.
// **************************************************************
//
// corrected by V. Grichine on 24/11/97
// corrected by L. Urban    on 27/05/98   ( other corrections come soon!)
//
 

#include "G4PAIenergyLoss.hh"
#include "G4PAIonisation.hh"
#include "G4EnergyLossTables.hh"

////////////////////////////////////////////////////////////////////////////
//
// Initialisation of static members 
//  ( this stuff should be defined later using RW ..........)
// contributing processes : ion.loss ->NUMBEROFPROCESSES is initialized
//   to 1 . YOU DO NOT HAVE TO CHANGE this variable for a 'normal' run.
// You have to change NUMBEROFPROCESSES  
// if you invent a new process contributing to the cont. energy loss,
//   NUMBEROFPROCESSES should be 2 in this case,
//  or for debugging purposes.
//  The NUMBEROFPROCESSES data member can be changed using the (public static)
//  functions Get/Set/Plus/MinusNUMBEROFPROCESSES (see G4hEnergyLoss.hh)

G4int G4PAIenergyLoss::NUMBEROFPROCESSES = 1 ;
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

G4PhysicsTable* G4PAIenergyLoss::theDEDXTable = NULL;

// G4PhysicsTable* G4PAIenergyLoss::fPAItransferBank = NULL ;

G4double G4PAIenergyLoss::Mass,
         G4PAIenergyLoss::taulow,
         G4PAIenergyLoss::tauhigh,
         G4PAIenergyLoss::ltaulow,
         G4PAIenergyLoss::ltauhigh;

G4double G4PAIenergyLoss::CutInRange = 0;

const G4double G4PAIenergyLoss::LowestKineticEnergy= 1.00*keV ;
const G4double G4PAIenergyLoss::HighestKineticEnergy= 100.*TeV ;
G4int G4PAIenergyLoss::TotBin ;
G4double G4PAIenergyLoss::RTable,G4PAIenergyLoss::LOGRTable;


// constructor and destructor
 
G4PAIenergyLoss::G4PAIenergyLoss(const G4String& processName)
   : G4VContinuousDiscreteProcess (processName),
     dToverTini(0.20),   // max.relative range loss in one Step = 20%
     // LowestKineticEnergy(1.00*keV),
     // HighestKineticEnergy(100.*TeV),
     MaxExcitationNumber (1.e6),
     probLimFluct (0.01),
     nmaxDirectFluct (100),
     nmaxCont1(4),
     nmaxCont2(16),
     theElectron ( G4Electron::Electron() ),
     theProton ( G4Proton::Proton() ),
     theAntiProton ( G4AntiProton::AntiProton() )
{
     theLossTable = NULL ;
     lastMaterial = NULL ;

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

    if(CounterOfProcess == NUMBEROFPROCESSES)
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

    if(CounterOfProcess == NUMBEROFPROCESSES)
    {
  // create tables
      if(theDEDXpbarTable)
      { theDEDXpbarTable->clearAndDestroy();
        delete theDEDXpbarTable; }
      theDEDXpbarTable = new G4PhysicsTable(numOfMaterials);
      theDEDXTable = theDEDXpbarTable;
    }
  }


  if(CounterOfProcess == NUMBEROFPROCESSES)
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
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);   

      // loop for the kinetic energy
   
      for (G4int i=0; i<TotBin; i++)
  
      {
        LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;      
//     here comes the sum of the different tables created by the  
//     processes (ionisation,etc...)              

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
 
  //  reset counter to zero 

  if( Charge >0.)    CounterOfpProcess=0 ;
  else  CounterOfpbarProcess=0 ;

    // Build range table

    BuildRangeTable( aParticleType);  
 
    // Build coeff tables for the energy loss calculation

    BuildRangeCoeffATable( aParticleType);
    BuildRangeCoeffBTable( aParticleType);
    BuildRangeCoeffCTable( aParticleType);
 
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

/////////////////////////////////////////////////////////////////////////
//
// Build range table from the energy loss table
      
void G4PAIenergyLoss::BuildRangeTable( const G4ParticleDefinition& aParticleType)

{
    G4double Charge = aParticleType.GetPDGCharge() ;

    Mass = proton_mass_c2; 

//  create table

    const G4MaterialTable* theMaterialTable=
                                  G4Material::GetMaterialTable();


    G4int numOfMaterials = theMaterialTable->length();

  G4PhysicsTable* theRangeTable;

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

///////////////////////////////////////////////////////////////////
//
//  create range vector for a material
//

void G4PAIenergyLoss::BuildRangeVector(G4int materialIndex,
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
  tau1 = t1/Mass ;
  sqtau1 = sqrt(tau1) ;
  ca = (4.*loss2-loss1)/sqtau1 ;
  cb = (2.*loss1-4.*loss2)/tau1 ;
  cba = cb/ca ;
  taulim = tlim/Mass ;
  ltaulim = log(taulim) ;
  ltaumax = log(HighestKineticEnergy/Mass) ;

  // loop for kinetic energy

  for (G4int i=0; i<TotBin; i++)
  {
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
        nbin = (G4int)(maxbint*(tau-tau1)/(taulim-tau1)) ;
        if(nbin<1) nbin = 1;
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
        nbin = (G4int)(maxbint*(ltauhigh-ltaulow)/(ltaumax-ltaulow)) ;
        if(nbin<1) nbin= 1 ;
        Value += RangeIntLog(physicsVector,nbin);
      }
    }
      rangeVector->PutValue(i,Value);
  }
}    

///////////////////////////////////////////////////////////////////
//
//  num. integration, linear binning
//

G4double G4PAIenergyLoss::RangeIntLin(G4PhysicsVector* physicsVector,
                                    G4int nbin)
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

//////////////////////////////////////////////////////////////////
//
//  num. integration, logarithmic binning
//

G4double G4PAIenergyLoss::RangeIntLog(G4PhysicsVector* physicsVector,
                                    G4int nbin)
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

//////////////////////////////////////////////////////////////////////
//
// returns the range of a particle in a given material
//

G4double G4PAIenergyLoss::OldGetRange(const G4DynamicParticle *aParticle,
                                G4Material *aMaterial)
{
  G4int materialIndex;
  G4double KineticEnergy,Range;
  G4double BigRange = DBL_MAX;
  G4bool isOut;
  G4double Charge = aParticle->GetDefinition()->GetPDGCharge() ;
  G4double Chargesquare = Charge*Charge ;
  G4double MassRatio = proton_mass_c2/aParticle->GetDefinition()->GetPDGMass() ;

  G4PhysicsTable* theRangeTable;

 if(Charge>0.) theRangeTable=theRangepTable ;

 else          theRangeTable=theRangepbarTable ;

  materialIndex = aMaterial->GetIndex();
  KineticEnergy = aParticle->GetKineticEnergy()*MassRatio;

  if(KineticEnergy<LowestKineticEnergy)
    Range = BigRange ;
 
  else
  {
      if(KineticEnergy>HighestKineticEnergy) Range = BigRange ;

      else Range = (*theRangeTable)[materialIndex]->GetValue(KineticEnergy,isOut) ;

  }
  return Range/Chargesquare ;

} 

///////////////////////////////////////////////////////////////////////
//
// Build tables of coefficients for the energy loss calculation
//

void G4PAIenergyLoss::
BuildRangeCoeffATable(const G4ParticleDefinition& aParticleType)
{
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

//  create table for coefficients "A"

  G4int numOfMaterials = theMaterialTable->length();
  G4double Charge = aParticleType.GetPDGCharge() ;

  G4PhysicsTable* theRangeTable;
  G4PhysicsTable* theRangeCoeffATable;
  if(Charge>0.)
  {
    if(thepRangeCoeffATable)
    { 
      thepRangeCoeffATable->clearAndDestroy() ;
      delete thepRangeCoeffATable ; 
    }
    thepRangeCoeffATable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffATable = thepRangeCoeffATable ;
    theRangeTable = theRangepTable ;
  }
  else  
  {
    if(thepbarRangeCoeffATable)
    { 
      thepbarRangeCoeffATable->clearAndDestroy() ;
      delete thepbarRangeCoeffATable ; 
    }
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


  for (G4int J=0; J<numOfMaterials; J++)    //  loop for materials
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
     
      if ( i==0 ) Rim = 0. ;
      else
      {
        Tim = Ti/RTable ;
        Rim = rangeVector->GetValue(Tim,isOut);
      }
      if ( i==(TotBin-1)) Rip = Ri ;
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

////////////////////////////////////////////////////////////////////////
//
// Build tables of coefficients for the energy loss calculation
//

void G4PAIenergyLoss::
BuildRangeCoeffBTable(const G4ParticleDefinition& aParticleType)
{
   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

//  create table for coefficients "B"


  G4int numOfMaterials = theMaterialTable->length();
  G4double Charge = aParticleType.GetPDGCharge() ;

  G4PhysicsTable* theRangeTable;
  G4PhysicsTable* theRangeCoeffBTable;
  if(Charge>0.)
  {
    if(thepRangeCoeffBTable)
    { 
      thepRangeCoeffBTable->clearAndDestroy();
      delete thepRangeCoeffBTable ; 
    }
    thepRangeCoeffBTable = new G4PhysicsTable(numOfMaterials);
    theRangeCoeffBTable = thepRangeCoeffBTable ;
    theRangeTable = theRangepTable ;
  }
  else
  {
    if(thepbarRangeCoeffBTable)
    { 
      thepbarRangeCoeffBTable->clearAndDestroy();
      delete thepbarRangeCoeffBTable; 
    }
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
   // create  vector

    G4int binmax=TotBin ;
    G4PhysicsLinearVector* aVector =  new G4PhysicsLinearVector(0.,binmax, TotBin);

   //  loop for kinetic energy

    Ti = LowestKineticEnergy ;

    G4PhysicsVector* rangeVector= (*theRangeTable)[J];
   
    for ( G4int i=0; i<TotBin; i++)
    {
      Ri = rangeVector->GetValue(Ti,isOut) ;
     
      if ( i==0 ) Rim = 0. ;
      else
      {
        Tim = Ti/RTable ;
        Rim = rangeVector->GetValue(Tim,isOut);
      }
      if ( i==(TotBin-1)) Rip = Ri ;
      else
      {
        Tip = Ti*RTable ;
        Rip = rangeVector->GetValue(Tip,isOut);
      }
      Value = (w1*Rip + w2*Ri + w3*Rim)/Ti;  

      aVector->PutValue(i,Value);

      Ti = RTable*Ti ;
    }  
    theRangeCoeffBTable->insert(aVector) ;
  } 
}

///////////////////////////////////////////////////////////////////////////
//
// Build tables of coefficients for the energy loss calculation
//

void G4PAIenergyLoss::BuildRangeCoeffCTable(
                            const G4ParticleDefinition& aParticleType)
{
   const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();

//  create table for coefficients "C"

  G4int numOfMaterials = theMaterialTable->length();
  G4double Charge = aParticleType.GetPDGCharge() ;

  G4PhysicsTable* theRangeTable;
  G4PhysicsTable* theRangeCoeffCTable;
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
     
      if ( i==0 ) Rim = 0. ;
      else
      {
        Tim = Ti/RTable ;
        Rim = rangeVector->GetValue(Tim,isOut);
      }
      if ( i==(TotBin-1)) Rip = Ri ;
      else
      {
        Tip = Ti*RTable ;
        Rip = rangeVector->GetValue(Tip,isOut);
      }
      Value = w1*Rip + w2*Ri + w3*Rim ;

      aVector->PutValue(i,Value);

      Ti = RTable*Ti ;
    }
    theRangeCoeffCTable->insert(aVector) ;
  } 
}

//
//
/////////////////////////////////////////////////////////////////////////
   

























