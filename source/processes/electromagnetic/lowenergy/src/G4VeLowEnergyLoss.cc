//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VeLowEnergyLoss.cc,v 1.16 2001-11-07 20:47:30 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// --------------------------------------------------------------
//
// Modifications:
// 20/09/00 update fluctuations V.Ivanchenko
// 22/11/00 minor fix in fluctuations V.Ivanchenko
// 10/05/01  V.Ivanchenko Clean up againist Linux compilation with -Wall
// 22/05/01  V.Ivanchenko Update range calculation
//
// --------------------------------------------------------------

#include "G4VeLowEnergyLoss.hh"

G4double     G4VeLowEnergyLoss::ParticleMass ;                
G4double     G4VeLowEnergyLoss::taulow       ;                
G4double     G4VeLowEnergyLoss::tauhigh      ;                
G4double     G4VeLowEnergyLoss::ltaulow       ;                
G4double     G4VeLowEnergyLoss::ltauhigh      ;                


G4bool       G4VeLowEnergyLoss::rndmStepFlag   = false;
G4bool       G4VeLowEnergyLoss::EnlossFlucFlag = true;
G4double     G4VeLowEnergyLoss::dRoverRange    = 20*perCent;
G4double     G4VeLowEnergyLoss::finalRange     = 200*micrometer;
G4double     G4VeLowEnergyLoss::c1lim = dRoverRange ;
G4double     G4VeLowEnergyLoss::c2lim = 2.*(1.-dRoverRange)*finalRange ;
G4double     G4VeLowEnergyLoss::c3lim = -(1.-dRoverRange)*finalRange*finalRange;


//    

G4VeLowEnergyLoss::G4VeLowEnergyLoss()
                   :G4VContinuousDiscreteProcess("No Name Loss Process"),
     lastMaterial(0),
     nmaxCont1(4),
     nmaxCont2(16)
{
  G4Exception("G4VeLowEnergyLoss:: default constructor is called");
}

//    

G4VeLowEnergyLoss::G4VeLowEnergyLoss(const G4String& aName , G4ProcessType aType)
                  : G4VContinuousDiscreteProcess(aName, aType),
     lastMaterial(0),
     nmaxCont1(4),
     nmaxCont2(16)
{
}

//    

G4VeLowEnergyLoss::~G4VeLowEnergyLoss()
{
}

//    

G4VeLowEnergyLoss::G4VeLowEnergyLoss(G4VeLowEnergyLoss& right)
                  : G4VContinuousDiscreteProcess(right),
     lastMaterial(0),
     nmaxCont1(4),
     nmaxCont2(16)
{
}

//    

G4PhysicsTable* G4VeLowEnergyLoss::BuildRangeTable(
        G4PhysicsTable* theDEDXTable,G4PhysicsTable* theRangeTable,            
        G4double lowestKineticEnergy,G4double highestKineticEnergy,G4int TotBin)
// Build range table from the energy loss table
{
   
   G4int numOfMaterials = G4Material::GetNumberOfMaterials();

   if(theRangeTable)
   { theRangeTable->clearAndDestroy();
     delete theRangeTable; }
   theRangeTable = new G4PhysicsTable(numOfMaterials);

   // loop for materials

   for (G4int J=0;  J<numOfMaterials; J++)
   {
     G4PhysicsLogVector* aVector;
     aVector = new G4PhysicsLogVector(lowestKineticEnergy,
                              highestKineticEnergy,TotBin);
     BuildRangeVector(theDEDXTable,lowestKineticEnergy,highestKineticEnergy,
                      TotBin,J,aVector);
     theRangeTable->insert(aVector);
   }
   return theRangeTable ;
}   

//    

void G4VeLowEnergyLoss::BuildRangeVector(G4PhysicsTable* theDEDXTable,
                                         G4double lowestKineticEnergy,
                                         G4double,
                                         G4int TotBin,
                                         G4int materialIndex,
                                         G4PhysicsLogVector* rangeVector)

//  create range vector for a material
{
  G4bool isOut;
  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];
  G4double energy1 = lowestKineticEnergy;
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

//    

G4double G4VeLowEnergyLoss::RangeIntLin(G4PhysicsVector* physicsVector,
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

//    

G4double G4VeLowEnergyLoss::RangeIntLog(G4PhysicsVector* physicsVector,
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


//    

G4PhysicsTable* G4VeLowEnergyLoss::BuildLabTimeTable(G4PhysicsTable* theDEDXTable,
                                     G4PhysicsTable* theLabTimeTable,
                                     G4double lowestKineticEnergy,
                                     G4double highestKineticEnergy,G4int TotBin)
                            
{

  G4int numOfMaterials = G4Material::GetNumberOfMaterials();
 
  if(theLabTimeTable)
  { theLabTimeTable->clearAndDestroy();
    delete theLabTimeTable; }
  theLabTimeTable = new G4PhysicsTable(numOfMaterials);


  for (G4int J=0;  J<numOfMaterials; J++)
  {
    G4PhysicsLogVector* aVector;

    aVector = new G4PhysicsLogVector(lowestKineticEnergy,
                            highestKineticEnergy,TotBin);

    BuildLabTimeVector(theDEDXTable,
              lowestKineticEnergy,highestKineticEnergy,TotBin,J,aVector);
    theLabTimeTable->insert(aVector);


  }
  return theLabTimeTable ;
}

//    

G4PhysicsTable* G4VeLowEnergyLoss::BuildProperTimeTable(G4PhysicsTable* theDEDXTable,
                                     G4PhysicsTable* theProperTimeTable,
                                     G4double lowestKineticEnergy,
                                     G4double highestKineticEnergy,G4int TotBin)
                            
{

  G4int numOfMaterials = G4Material::GetNumberOfMaterials();
 
  if(theProperTimeTable)
  { theProperTimeTable->clearAndDestroy();
    delete theProperTimeTable; }
  theProperTimeTable = new G4PhysicsTable(numOfMaterials);


  for (G4int J=0;  J<numOfMaterials; J++)
  {
    G4PhysicsLogVector* aVector;

    aVector = new G4PhysicsLogVector(lowestKineticEnergy,
                            highestKineticEnergy,TotBin);

    BuildProperTimeVector(theDEDXTable,
              lowestKineticEnergy,highestKineticEnergy,TotBin,J,aVector);
    theProperTimeTable->insert(aVector);


  }
  return theProperTimeTable ;
}

//    

void G4VeLowEnergyLoss::BuildLabTimeVector(G4PhysicsTable* theDEDXTable,
                                    G4double lowestKineticEnergy,
                                    G4double highestKineticEnergy,G4int TotBin,
                                    G4int materialIndex, G4PhysicsLogVector* timeVector)
//  create lab time vector for a material
{

  G4int nbin=100;
  G4bool isOut;
  G4double tlim=5.*keV,parlowen=0.4,ppar=0.5-parlowen ;
  G4double losslim,clim,taulim,timelim,ltaulim,ltaumax,
           LowEdgeEnergy,tau,Value ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];

  // low energy part first...
  losslim = physicsVector->GetValue(tlim,isOut);
  taulim=tlim/ParticleMass ;
  clim=sqrt(ParticleMass*tlim/2.)/(c_light*losslim*ppar) ;
  ltaulim = log(taulim);
  ltaumax = log(highestKineticEnergy/ParticleMass) ;

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

//    

void G4VeLowEnergyLoss::BuildProperTimeVector(G4PhysicsTable* theDEDXTable,
                                    G4double lowestKineticEnergy,
                                    G4double highestKineticEnergy,G4int TotBin,
                                    G4int materialIndex, G4PhysicsLogVector* timeVector)
//  create proper time vector for a material
{
  G4int nbin=100;
  G4bool isOut;
  G4double tlim=5.*keV,parlowen=0.4,ppar=0.5-parlowen ;
  G4double losslim,clim,taulim,timelim,ltaulim,ltaumax,
           LowEdgeEnergy,tau,Value ;

  G4PhysicsVector* physicsVector= (*theDEDXTable)[materialIndex];
  //const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // low energy part first...
  losslim = physicsVector->GetValue(tlim,isOut);
  taulim=tlim/ParticleMass ;
  clim=sqrt(ParticleMass*tlim/2.)/(c_light*losslim*ppar) ;
  ltaulim = log(taulim);
  ltaumax = log(highestKineticEnergy/ParticleMass) ;

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

//    

G4double G4VeLowEnergyLoss::LabTimeIntLog(G4PhysicsVector* physicsVector,
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

//    

G4double G4VeLowEnergyLoss::ProperTimeIntLog(G4PhysicsVector* physicsVector,
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

//    

G4PhysicsTable* G4VeLowEnergyLoss::BuildInverseRangeTable(G4PhysicsTable* theRangeTable,
                                   G4PhysicsTable* theRangeCoeffATable,
                                   G4PhysicsTable* theRangeCoeffBTable,
                                   G4PhysicsTable* theRangeCoeffCTable,
                                   G4PhysicsTable* theInverseRangeTable,
                                   G4double lowestKineticEnergy,
                                   G4double highestKineticEnergy,G4int TotBin)
// Build inverse table of the range table
{
  G4double SmallestRange,BiggestRange ;
  G4bool isOut ;

  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

    if(theInverseRangeTable)
    { theInverseRangeTable->clearAndDestroy();
      delete theInverseRangeTable; }
    theInverseRangeTable = new G4PhysicsTable(numOfMaterials);

  // loop for materials
  for (G4int J=0;  J<numOfMaterials; J++)
  {
    SmallestRange = (*theRangeTable)(J)->
                       GetValue(lowestKineticEnergy,isOut) ;
    BiggestRange = (*theRangeTable)(J)->
                       GetValue(highestKineticEnergy,isOut) ;
    G4PhysicsLogVector* aVector;
    aVector = new G4PhysicsLogVector(SmallestRange,
                            BiggestRange,TotBin);

    InvertRangeVector(theRangeTable,
                      theRangeCoeffATable,
                      theRangeCoeffBTable,
                      theRangeCoeffCTable,
         lowestKineticEnergy,highestKineticEnergy,TotBin,J, aVector);

    theInverseRangeTable->insert(aVector);
  }
  return theInverseRangeTable ;
}

//    

void G4VeLowEnergyLoss::InvertRangeVector(G4PhysicsTable* theRangeTable,
                              G4PhysicsTable* theRangeCoeffATable,
                              G4PhysicsTable* theRangeCoeffBTable,
                              G4PhysicsTable* theRangeCoeffCTable,
                              G4double lowestKineticEnergy,
                              G4double highestKineticEnergy,G4int TotBin,
                              G4int  materialIndex,G4PhysicsLogVector* aVector)
//  invert range vector for a material
{
  G4double LowEdgeRange,A,B,C,discr,KineticEnergy ;
  G4double RTable = exp(log(highestKineticEnergy/lowestKineticEnergy)/TotBin) ;
  G4double Tbin = lowestKineticEnergy/RTable ;
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
      KineticEnergy = lowestKineticEnergy ;
    else if(binnumber == TotBin-1)
      KineticEnergy = highestKineticEnergy ;
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

//    

G4PhysicsTable* G4VeLowEnergyLoss::BuildRangeCoeffATable(G4PhysicsTable* theRangeTable,
                                      G4PhysicsTable* theRangeCoeffATable,
                                      G4double lowestKineticEnergy,
                             G4double highestKineticEnergy,G4int TotBin)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "A"
{

  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  if(theRangeCoeffATable)
  { theRangeCoeffATable->clearAndDestroy();
    delete theRangeCoeffATable; }
  theRangeCoeffATable = new G4PhysicsTable(numOfMaterials);

  G4double RTable = exp(log(highestKineticEnergy/lowestKineticEnergy)/TotBin) ;
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
    Ti = lowestKineticEnergy ;
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
  return theRangeCoeffATable ;
}

//    
  
G4PhysicsTable* G4VeLowEnergyLoss::BuildRangeCoeffBTable(G4PhysicsTable* theRangeTable,
                                      G4PhysicsTable* theRangeCoeffBTable,
                                      G4double lowestKineticEnergy,
                             G4double highestKineticEnergy,G4int TotBin)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "B"
{

  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  if(theRangeCoeffBTable)
  { theRangeCoeffBTable->clearAndDestroy();
    delete theRangeCoeffBTable; }
  theRangeCoeffBTable = new G4PhysicsTable(numOfMaterials);

  G4double RTable = exp(log(highestKineticEnergy/lowestKineticEnergy)/TotBin) ;
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
    Ti = lowestKineticEnergy ;
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
  return theRangeCoeffBTable ;
}

//    

G4PhysicsTable* G4VeLowEnergyLoss::BuildRangeCoeffCTable(G4PhysicsTable* theRangeTable,
                                      G4PhysicsTable* theRangeCoeffCTable,
                                      G4double lowestKineticEnergy,
                             G4double highestKineticEnergy,G4int TotBin)
// Build tables of coefficients for the energy loss calculation
//  create table for coefficients "C"
{

  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  if(theRangeCoeffCTable)
  { theRangeCoeffCTable->clearAndDestroy();
    delete theRangeCoeffCTable; }
  theRangeCoeffCTable = new G4PhysicsTable(numOfMaterials);

  G4double RTable = exp(log(highestKineticEnergy/lowestKineticEnergy)/TotBin) ;
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
    Ti = lowestKineticEnergy ;
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
  return theRangeCoeffCTable ;
}

//    

G4double G4VeLowEnergyLoss::GetLossWithFluct(const G4DynamicParticle* aParticle,
                                               G4Material* aMaterial,
                                               G4double MeanLoss,
                                               G4double step)
//  calculate actual loss from the mean loss
//  The model used to get the fluctuation is essentially the same as in Glandz in Geant3.
{
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
  G4double threshold,w1,w2,C,
           beta2,suma,e0,loss,lossc,w;
  G4double a1,a2,a3;
  G4int p1,p2,p3;
  G4int nb;
  G4double Corrfac, na,alfa,rfac,namean,sa,alfa1,ea,sea;
  //  G4double dp1;
  G4double dp3;
  G4double siga ;

  // shortcut for very very small loss 
  if(MeanLoss < minLoss) return MeanLoss ;

  // get particle data
  G4double Tkin   = aParticle->GetKineticEnergy();

  //  G4cout << "MGP -- Fluc Tkin " << Tkin/keV << " keV "  << " MeanLoss = " << MeanLoss/keV << G4endl;

  threshold = G4Electron::Electron()->GetEnergyThreshold(aMaterial);
  G4double rmass = electron_mass_c2/ParticleMass;
  G4double tau   = Tkin/ParticleMass, tau1 = tau+1., tau2 = tau*(tau+2.);
  G4double Tm    = 2.*electron_mass_c2*tau2/(1.+2.*tau1*rmass+rmass*rmass);

  // G4cout << "MGP Particle mass " << ParticleMass/MeV << " Tm " << Tm << G4endl;

  if (Tm <= ipotFluct) Tm = ipotFluct ;
  
  if(Tm > threshold) Tm = threshold;
  beta2 = tau2/(tau1*tau1);

  // Gaussian fluctuation ?
  if(MeanLoss >= kappa*Tm)
  {
    G4double electronDensity = aMaterial->GetElectronDensity() ;
    siga = sqrt(Tm*(1.0-0.5*beta2)*step*
                factor*electronDensity/beta2) ;
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
      a3  = 0. ;
    } 
  
  suma = a1+a2+a3;
  
  loss = 0. ;
  
  if(suma < sumaLim)             // very small Step
    {
      e0 = aMaterial->GetIonisation()->GetEnergy0fluct();
      // G4cout << "MGP e0 = " << e0/keV << G4endl;
      
      if(Tm == ipotFluct)
	{
	  a3 = MeanLoss/e0;
	  
	  if(a3>alim)
	    {
	      siga=sqrt(a3) ;
	      p3 = G4std::max(0,G4int(G4RandGauss::shoot(a3,siga)+0.5));
	    }
	  else p3 = G4Poisson(a3);
	  
	  loss = p3*e0 ;
	  
	  if(p3 > 0) loss += (1.-2.*G4UniformRand())*e0 ;
	  // G4cout << "MGP very small step " << loss/keV << G4endl;
	}
      else
	{
	  //	  G4cout << "MGP old Tm = " << Tm << " " << ipotFluct << " " << e0 << G4endl;
	  Tm = Tm-ipotFluct+e0 ;

	  // MGP ---- workaround to avoid log argument<0, TO BE CHECKED
	  if (Tm <= 0.) 
	    { 
	      loss = MeanLoss;
	      p3 = 0;
	      // G4cout << "MGP correction loss = MeanLoss " << loss/keV << G4endl;
	    }
	  else
	    {
	      a3 = MeanLoss*(Tm-e0)/(Tm*e0*log(Tm/e0));

	      // G4cout << "MGP new Tm = " << Tm << " " << ipotFluct << " " << e0 << " a3= " << a3 << G4endl;
	      
	      if(a3>alim)
		{
		  siga=sqrt(a3) ;
		  p3 = G4std::max(0,G4int(G4RandGauss::shoot(a3,siga)+0.5));
		}
	      else
		p3 = G4Poisson(a3);
	      //G4cout << "MGP p3 " << p3 << G4endl;

	    }
	      
	  if(p3 > 0)
	    {
	      w = (Tm-e0)/Tm ;
	      if(p3 > nmaxCont2)
		{
		  // G4cout << "MGP dp3 " << dp3 << " p3 " << p3 << " " << nmaxCont2 << G4endl;
		  dp3 = G4double(p3) ;
		  Corrfac = dp3/G4double(nmaxCont2) ;
		  p3 = nmaxCont2 ;
		}
	      else
		Corrfac = 1. ;
	      
	      for(G4int i=0; i<p3; i++) loss += 1./(1.-w*G4UniformRand()) ;
	      loss *= e0*Corrfac ; 
	      // G4cout << "MGP Corrfac = " << Corrfac << " e0 = " << e0/keV << " loss = " << loss/keV << G4endl;
	    }        
	}
    }
    
  else                              // not so small Step
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
		  dp3        = G4double(p3);
		  rfac       = dp3/(G4double(nmaxCont2)+dp3);
		  namean     = G4double(p3)*rfac;
		  sa         = G4double(nmaxCont1)*rfac;
		  na         = G4RandGauss::shoot(namean,sa);
		  if (na > 0.)
		    {
		      alfa   = w1*G4double(nmaxCont2+p3)/(w1*G4double(nmaxCont2)+G4double(p3));
		      alfa1  = alfa*log(alfa)/(alfa-1.);
		      ea     = na*ipotFluct*alfa1;
		      sea    = ipotFluct*sqrt(na*(alfa-alfa1*alfa1));
		      lossc += G4RandGauss::shoot(ea,sea);
		    }
		}
	      
	      nb = G4int(G4double(p3)-na);
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

//    
   
