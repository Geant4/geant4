// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IeplusAnnihilation.cc,v 1.4 1999-05-03 11:04:14 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4IeplusAnnihilation process --------
//                   by Michel Maire, 7 July 1996
// **************************************************************
// --------------------------------------------------------------
// ************************************************************
// It is the first implementation of the 
//          eplusANNIHILATION  PROCESS
//   using an INTEGRAL APPROACH instead of the differential
//   one used in the standard implementation .
// ************************************************************
//                by Laszlo Urban, 23 June 1998
// -----------------------------------------------------------
// 28/10/28: some cleanup , L.Urban

#include "G4IeplusAnnihilation.hh"
#include "G4UnitsTable.hh"   
 
// constructor
G4IeplusAnnihilation::G4IeplusAnnihilation(const G4String& processName)
  : G4IVRestDiscreteProcess (processName),
    LowestEnergyLimit ( 10*keV),      // initialization
    HighestEnergyLimit( 10*TeV),
    NumberOfBuildPhysicsTableCalls(0),
    NumbBinTable(100)
 {
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
     G4cout << "LowestEnergy: " << LowestEnergyLimit/keV << "keV ";
     G4cout << "HighestEnergy: " << HighestEnergyLimit/TeV << "TeV " << endl;
   }
     theCrossSectionTable = NULL;
     theMeanFreePathTable = NULL; 
     theMeanFreePathTable = NULL ;
     theNlambdaTable = NULL;
     theInverseNlambdaTable = NULL;
     theCoeffATable = NULL ;
     theCoeffBTable = NULL ;
     theCoeffCTable = NULL ;

     LowestKineticEnergy = LowestEnergyLimit ;
     HighestKineticEnergy= HighestEnergyLimit;
     TotBin = NumbBinTable ;
     RTable = exp(log(HighestKineticEnergy/LowestKineticEnergy)/TotBin) ;
 }
 
// destructor
 
G4IeplusAnnihilation::~G4IeplusAnnihilation()
{
   if (theCrossSectionTable) {
      theCrossSectionTable->clearAndDestroy();
      delete theCrossSectionTable;
   }

   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }
     if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
     }
     if (theNlambdaTable) {
        theNlambdaTable->clearAndDestroy();
        delete theNlambdaTable;
     }
     if (theInverseNlambdaTable) {
        theInverseNlambdaTable->clearAndDestroy();
        delete theInverseNlambdaTable;
     }
     if (theCoeffATable) {
        theCoeffATable->clearAndDestroy();
        delete theCoeffATable;
     }
     if (theCoeffBTable) {
        theCoeffBTable->clearAndDestroy();
        delete theCoeffBTable;
     }
     if (theCoeffCTable) {
        theCoeffCTable->clearAndDestroy();
        delete theCoeffCTable;
     }
}

void G4IeplusAnnihilation::SetPhysicsTableBining(G4double lowE,
                                        G4double highE, G4int nBins)
{
  LowestEnergyLimit = lowE; HighestEnergyLimit = highE; NumbBinTable = nBins;
}

void G4IeplusAnnihilation::BuildPhysicsTable(const G4ParticleDefinition& PositronType)

// Build microscopic total cross section tables and mean free path table
{
 NumberOfBuildPhysicsTableCalls += 1 ;
 if(NumberOfBuildPhysicsTableCalls == 1)
 { ; }
 else
 {

   G4double LowEdgeEnergy, Value;
   G4PhysicsLogVector* ptrVector;

// Build microscopic cross section tables for the e+e- annihilation

   if (theCrossSectionTable) {
          theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable; }

   theCrossSectionTable = new G4PhysicsTable( G4Element::GetNumberOfElements()) ;
   const G4ElementTable* theElementTable = G4Element::GetElementTable() ;
   G4double AtomicNumber;
   G4int J;

   for ( J=0 ; J < G4Element::GetNumberOfElements(); J++ )  
       { 
        //create physics vector then fill it ....
        ptrVector = new G4PhysicsLogVector(LowestEnergyLimit, HighestEnergyLimit,
                                           NumbBinTable ) ;
        AtomicNumber = (*theElementTable)(J)->GetZ();
 
        for ( G4int i = 0 ; i < NumbBinTable ; i++ )      
           {
             LowEdgeEnergy = ptrVector->GetLowEdgeEnergy( i ) ;
             Value = ComputeMicroscopicCrossSection( LowEdgeEnergy, AtomicNumber);  
             ptrVector->PutValue( i , Value ) ;
           }

        theCrossSectionTable->insertAt( J , ptrVector ) ;

      }

// Build mean free path table for the e+e- annihilation

   if (theMeanFreePathTable) {
          theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable; }

   theMeanFreePathTable = new G4PhysicsTable( G4Material::GetNumberOfMaterials() );

   //*******************!!!!!!!!!!!!!!!!!!!!!!********************
   theMeanFreePathTable = theMeanFreePathTable ;
   //*******************!!!!!!!!!!!!!!!!!!!!!!********************

   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;
   G4Material* material;

   for ( J=0 ; J < G4Material::GetNumberOfMaterials(); J++ )  
       { 
        //create physics vector then fill it ....
        ptrVector = new G4PhysicsLogVector(LowestEnergyLimit, HighestEnergyLimit,
                                           NumbBinTable ) ;
        material = (*theMaterialTable)(J);
 
        for ( G4int i = 0 ; i < NumbBinTable ; i++ )      
           {
             LowEdgeEnergy = ptrVector->GetLowEdgeEnergy( i ) ;
             Value = ComputeMeanFreePath( LowEdgeEnergy, material);  
             ptrVector->PutValue( i , Value ) ;

        theMeanFreePathTable->insertAt( J , ptrVector ) ;

           }

      }

    const G4ParticleDefinition& aParticleType = PositronType ;

    BuildNlambdaTable(aParticleType) ;

    BuildCoeffATable(aParticleType) ;
    BuildCoeffBTable(aParticleType) ;
    BuildCoeffCTable(aParticleType) ;

    BuildInverseNlambdaTable(aParticleType) ;

    G4int printflag = 0 ;
    if(printflag>0)
      TestOfInversion(aParticleType,printflag) ;

    NumberOfBuildPhysicsTableCalls = 0 ;

    PrintInfoDefinition() ;

 }
}

void G4IeplusAnnihilation::TestOfInversion(
                             const G4ParticleDefinition& aParticleType,
                                   G4int printflag)
{
  G4double T,Nlambda,Tprime,delta,del,sum,delmean,Tdelta ;
  G4bool isOut ;

  const G4MaterialTable* theMaterialTable=
                          G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length() ;

  G4cout.setf(ios::scientific, ios::floatfield) ;
  if(printflag>1)
  {
  G4cout << endl;
  G4cout << "  particle=" << aParticleType.GetParticleName() << endl;
  G4cout << "----------------------" << endl;
  }
  for (G4int J=0; J<numOfMaterials; J++)
  {
   if(printflag>1)
   {
    G4cout << endl;
    G4cout << " material = " << (*theMaterialTable)[J]->GetName() << endl;
    G4cout << " mat.ind.=" << J << "  T           Nlambda         Tprime"
         << "         (Tprime-T)/T(%)" << endl ;
   }

    G4PhysicsLogVector* aVector ;

    aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                                    HighestKineticEnergy,TotBin) ;
    delta = 0. ;
    delmean = 0.;
    sum = 0.;
    Tdelta = 0. ;
    for (G4int i=0; i<TotBin-1; i++)
    {
      T = aVector->GetLowEdgeEnergy(i) ;

      Nlambda = (*theNlambdaTable)[J]->GetValue(T,isOut) ;

     if(Nlambda>0.)
     {
      Tprime = (*theInverseNlambdaTable)[J]->GetValue(Nlambda,isOut) ;

      if((Nlambda>0.)&&(i<(TotBin-1)))
      {
        del = 100.*(Tprime-T)/T ;
        sum += 1.;
        delmean += abs(del);
        if(abs(del)>abs(delta))
        {
           delta = del ;
           Tdelta = T ;
        }
      }
    if(printflag>1)
    {
      G4cout << setw(18) << setprecision(6) << T << "  " <<
              setw(14)<< setprecision(6) << Nlambda << "  " <<
      setw(14) << setprecision(6) << Tprime << "      " <<
              setw(12) << setprecision(3) << del << endl;
    }
     }
    }
    if(printflag>0)
    {
    G4cout << endl;
    G4cout << "G4IeplusAnnihilation::TestOfInversion (T->Nlambda->Tprime) " << endl
;
    G4cout << "particle= " << aParticleType.GetParticleName() <<
            "   material= " << (*theMaterialTable)[J]->GetName() << endl ;
    G4cout << "max (Tprime-T)/T in % =" << setw(10) << setprecision(3) << delta ;
    G4cout << "  at a kinetic energy " << setw(10) << setprecision(3) <<
            Tdelta/MeV << "  MeV" << endl;
    delmean /= sum ;
    G4cout << "mean rel.diff. (Tprime-T)/T=" << setw(10) <<
            setprecision(3) << delmean <<
            " % (mean is calculated in abs. value)" << endl;
    G4cout << endl;
    }
  }
}

void G4IeplusAnnihilation::BuildNlambdaTable(
                             const G4ParticleDefinition& aParticleType)
{
  const G4MaterialTable* theMaterialTable=
                          G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length() ;

  if(theNlambdaTable)
  {  theNlambdaTable->clearAndDestroy();
     delete theNlambdaTable ;           }

  theNlambdaTable = new G4PhysicsTable(numOfMaterials) ;

  // loop for materials
  for (G4int J=0; J<numOfMaterials; J++)
  {
    G4PhysicsLogVector* aVector ;

    aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                                    HighestKineticEnergy,TotBin) ;

    BuildNlambdaVector(aParticleType,J,aVector) ;

    theNlambdaTable->insert(aVector) ;
  }

}

void G4IeplusAnnihilation::BuildNlambdaVector(
                               const G4ParticleDefinition& aParticleType,
                                       G4int materialIndex,
                                       G4PhysicsLogVector* nlambdaVector)
{
  G4double LowEdgeEnergy,T,Tlast,dEdx,Value,Vlast,u,du,coeff,l ;
  const G4int nbin = 20 ;
  G4bool isOut ;
  const G4double small = 1.e-100;
  const G4double plowloss = 0.5 ;  //this should be a data member of en.loss!
  const G4double lmin=1.e-100,lmax=1.e100;

  const G4MaterialTable* theMaterialTable=
                          G4Material::GetMaterialTable();

  //first integral from 0. to LowestKineticEnergy
  // here assumed that the threshold energy for the process >=
  //                LowestKineticEnergy
  dEdx = G4EnergyLossTables::GetPreciseDEDX(&aParticleType,
                                            LowestKineticEnergy,
                                       (*theMaterialTable)[materialIndex]);
  Value = LowestKineticEnergy/(dEdx*BIGSTEP*(1.-plowloss)) ;
  if(Value<small)
     Value = 0. ;

  nlambdaVector->PutValue(0,Value) ;

  Tlast = LowestKineticEnergy ;
  Vlast = Value ;

  // loop for kinetic energy

  for (G4int i=1; i<TotBin; i++)
  {
    LowEdgeEnergy = nlambdaVector->GetLowEdgeEnergy(i) ;

    u = log(LowEdgeEnergy/Tlast) ;
    du = u/nbin ;

    u = -du ;
    Value = 0. ;
    for(G4int n=0; n<=nbin; n++)
    {
      u += du ;
      T = Tlast*exp(u) ;

      if((n==0)||(n==nbin))
       coeff=0.5 ;
      else
       coeff=1.0 ;

      l = (*theMeanFreePathTable)[materialIndex]->GetValue(T,isOut);
      if((l>lmin) && (l<lmax))
        Value += coeff*T/(G4EnergyLossTables::GetPreciseDEDX(&aParticleType,
                                   T,(*theMaterialTable)[materialIndex])*l);
    }

    Value *= du ;
    Value += Vlast ;
    if(Value<small)
       Value = 0. ;
    nlambdaVector->PutValue(i,Value) ;

    Tlast = LowEdgeEnergy ;
    Vlast = Value ;

  }

}

void G4IeplusAnnihilation::BuildCoeffATable(
                             const G4ParticleDefinition& aParticleType)
{
   const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();

//  create table for coefficients "A"

  G4int numOfMaterials = theMaterialTable->length();

  if(theCoeffATable)

    { theCoeffATable->clearAndDestroy();
      delete theCoeffATable; }
  theCoeffATable = new G4PhysicsTable(numOfMaterials);

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
   G4PhysicsLinearVector* aVector =
                            new G4PhysicsLinearVector(0.,binmax, TotBin);

   //  loop for kinetic energy

   Ti = LowestKineticEnergy ;

    G4PhysicsVector* lVector= (*theNlambdaTable)[J];

   for ( G4int i=0; i<TotBin; i++)
   {
     Ri = lVector->GetValue(Ti,isOut) ;


     if ( i==0 )
        Rim = Ri/sqrt(RTable) ;
     else
     {
        Tim = Ti/RTable ;
        Rim = lVector->GetValue(Tim,isOut);
     }

     Tip = Ti*RTable ;
     Rip = lVector->GetValue(Tip,isOut);

     if( i < (TotBin-1))
        Value = (w1*Rip + w2*Ri + w3*Rim)/(Ti*Ti) ;
     else
     {
        Rip = Rip*Rip/Rim ;
        Value = (w1*Rip + w2*Ri + w3*Rim)/(Ti*Ti) ;
     }

     aVector->PutValue(i,Value);


     Ti = RTable*Ti ;

   }

   theCoeffATable->insert(aVector);

  }

}

void G4IeplusAnnihilation::BuildCoeffBTable(
                             const G4ParticleDefinition& aParticleType)
{
   const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();

//  create table for coefficients "B"


  G4int numOfMaterials = theMaterialTable->length();

  if(theCoeffBTable)
    { theCoeffBTable->clearAndDestroy();
      delete theCoeffBTable; }
  theCoeffBTable = new G4PhysicsTable(numOfMaterials);

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
   G4PhysicsLinearVector* aVector =
                           new G4PhysicsLinearVector(0.,binmax, TotBin);

   //  loop for kinetic energy

   Ti = LowestKineticEnergy ;

   G4PhysicsVector* lVector= (*theNlambdaTable)[J];

   for ( G4int i=0; i<TotBin; i++)
   {
     Ri = lVector->GetValue(Ti,isOut) ;

     if ( i==0 )
        Rim = Ri/sqrt(RTable) ;
     else
     {
        Tim = Ti/RTable ;
        Rim = lVector->GetValue(Tim,isOut);
     }

     Tip = Ti*RTable ;
     Rip = lVector->GetValue(Tip,isOut);

     if(i < (TotBin-1))
        Value = (w1*Rip + w2*Ri + w3*Rim)/Ti;
     else
     {
        Rip = Rip*Rip/Rim ;
        Value = (w1*Rip + w2*Ri + w3*Rim)/Ti ;
     }

     aVector->PutValue(i,Value);

     Ti = RTable*Ti ;

   }

   theCoeffBTable->insert(aVector);

  }

}

void G4IeplusAnnihilation::BuildCoeffCTable(
                             const G4ParticleDefinition& aParticleType)
{
   const G4MaterialTable* theMaterialTable=
                                G4Material::GetMaterialTable();

//  create table for coefficients "C"

  G4int numOfMaterials = theMaterialTable->length();

  if(theCoeffCTable)
    { theCoeffCTable->clearAndDestroy();
      delete theCoeffCTable; }
  theCoeffCTable = new G4PhysicsTable(numOfMaterials);

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
   G4PhysicsLinearVector* aVector =
                          new G4PhysicsLinearVector(0.,binmax, TotBin);

   //  loop for kinetic energy

   Ti = LowestKineticEnergy ;

   G4PhysicsVector* lVector= (*theNlambdaTable)[J];

   for ( G4int i=0; i<TotBin; i++)
   {
     Ri = lVector->GetValue(Ti,isOut) ;

     if ( i==0 )
        Rim = Ri/sqrt(RTable) ;
     else
     {
        Tim = Ti/RTable ;
        Rim = lVector->GetValue(Tim,isOut);
     }

     Tip = Ti*RTable ;
     Rip = lVector->GetValue(Tip,isOut);

     if(i < (TotBin-1))
        Value = w1*Rip + w2*Ri + w3*Rim ;
     else
     {
        Rip = Rip*Rip/Rim ;
        Value = w1*Rip + w2*Ri + w3*Rim ;
     }


     aVector->PutValue(i,Value);

     Ti = RTable*Ti ;

   }

   theCoeffCTable->insert(aVector);

  }

}
void G4IeplusAnnihilation::BuildInverseNlambdaTable(
                             const G4ParticleDefinition& aParticleType)
{
    G4double T,Smallest,Biggest,TT ;
    const G4double small = 1.e-10;
    G4bool isOut ;

//  create table

    const G4MaterialTable* theMaterialTable=
                                  G4Material::GetMaterialTable();

    G4int numOfMaterials = theMaterialTable->length();

    if(theInverseNlambdaTable)
    { theInverseNlambdaTable->clearAndDestroy();
      delete theInverseNlambdaTable; }
    theInverseNlambdaTable = new G4PhysicsTable(numOfMaterials);

// loop for materials

    for (G4int J=0;  J<numOfMaterials; J++)
    {
      T = LowestKineticEnergy ;
      do
      {
        Smallest = (*theNlambdaTable)[J]->
                           GetValue(T,isOut) ;
        T *= RTable ;
      } while (Smallest <= small) ;

      Biggest = (*theNlambdaTable)[J]->
                       GetValue(HighestKineticEnergy,isOut) ;
     //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Biggest *= 1.25 ;   

      // create vector
      G4PhysicsLogVector* aVector;

      aVector = new G4PhysicsLogVector(Smallest,
                            Biggest,TotBin);

      // fill the vector

      InvertNlambdaVector(aParticleType,J, aVector);

      // insert vector to the table

      theInverseNlambdaTable->insert(aVector);

    }
}

void G4IeplusAnnihilation::InvertNlambdaVector(
                           const G4ParticleDefinition& aParticleType,
                                       G4int materialIndex,
                                       G4PhysicsLogVector* nlambdaVector)
{
 G4double LowEdge,A,B,C,discr,KineticEnergy ;
 G4double Tbin = LowestKineticEnergy/RTable ;
 G4double bin = 0.0 ;
 G4int binnumber = -1 ;
 G4bool isOut ;

 //loop for Nlambda values
 for( G4int i=0; i<TotBin; i++)
 {
  LowEdge = nlambdaVector->GetLowEdgeEnergy(i) ;  //i.e. GetLowEdgeValue(i)

  if( bin < LowEdge )
  {
   do
   {
    binnumber += 1 ;
    Tbin *= RTable ;
    bin = (*theNlambdaTable)[materialIndex]->GetValue(Tbin,isOut) ;
   }
   while ((bin < LowEdge) && (binnumber < TotBin-2 )) ;
  }

  if(binnumber == 0)
    KineticEnergy = LowestKineticEnergy ;
  else if(binnumber == TotBin-1)
    KineticEnergy = HighestKineticEnergy ;
  else
  {
    A = (*(*theCoeffATable)(materialIndex))(binnumber-1) ;
    B = (*(*theCoeffBTable)(materialIndex))(binnumber-1) ;
    C = (*(*theCoeffCTable)(materialIndex))(binnumber-1) ;
    if(A==0.)
      KineticEnergy = (LowEdge -C )/B ;
    else
    {
      discr = B*B - 4.*A*(C-LowEdge);
      discr = discr>0. ? sqrt(discr) : 0.;
      KineticEnergy = 0.5*(discr-B)/A ;
    }

  }

  nlambdaVector->PutValue(i,KineticEnergy) ;

 }
}

G4double G4IeplusAnnihilation::ComputeMicroscopicCrossSection
                              (G4double PositKinEnergy, G4double AtomicNumber)
 
// Calculates the microscopic cross section of annihilation into two photons
// from the Heilter formula.
// GEANT4 internal units.

{
 static const G4double pi_rcl2 = pi*classic_electr_radius*classic_electr_radius;

 G4double gama = 1. + PositKinEnergy/electron_mass_c2;
 G4double gama2 = gama*gama, sqgama2 = sqrt(gama2-1.);

 return pi_rcl2*AtomicNumber*((gama2+4*gama+1.)*log(gama+sqgama2) - (gama+3.)*sqgama2) 
                            /((gama2-1.)*(gama+1.));
}
 
 
G4VParticleChange* G4IeplusAnnihilation::PostStepDoIt(const G4Track& aTrack,
                                                    const G4Step&  aStep)
//
// The secondaries Gamma energies are sampled using the Heitler cross section.
//  
// A modified version of the random number techniques of Butcher & Messel is used 
//    (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : The initial electron is assumed free and at rest.
//          
// Note 2 : The annihilation processes producing one or more than two photons are
//          ignored, as negligible compared to the two photons process.         
 
{
   aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();
     

   const G4DynamicParticle* aDynamicPositron = aTrack.GetDynamicParticle();
   G4double PositKinEnergy = aDynamicPositron->GetKineticEnergy();
   G4ParticleMomentum PositDirection = aDynamicPositron->GetMomentumDirection();

   aParticleChange.Initialize(aTrack);

  // Do not make anything if PositKinEnergy=0. , the annihilation then
  // should be performed by the AtRestDoIt!
   if(PositKinEnergy == 0.)
     return &aParticleChange ;

   G4double gama = 1. + PositKinEnergy/electron_mass_c2;
   G4double gamap1 = gama+1. , gamam1 = gama-1. , sqgrate = sqrt(gamam1/gamap1)/2. ,
                                                  sqg2m1  = sqrt(gamam1*gamap1);

   // limits of the energy sampling
   G4double epsil1 = 0.5 - sqgrate , epsil2 = 0.5 + sqgrate;
   G4double epsilqot = epsil2/epsil1;

   //
   // sample the energy rate of the created gammas 
   //
   G4double epsil, greject ;

   do {
        epsil = epsil1*pow(epsilqot,G4UniformRand());
        greject = 1. - epsil + (2*gama*epsil-1.)/(epsil*gamap1*gamap1);
   } while( greject < G4UniformRand() );

   //
   // scattered Gamma angles. ( Z - axis along the parent positron)
   //
   
   G4double cost = (epsil*gamap1-1.)/(epsil*sqg2m1) , sint = sqrt((1.+cost)*(1.-cost));


   G4double phi  = twopi * G4UniformRand() ;
   G4double dirx = sint*cos(phi) , diry = sint*sin(phi) , dirz = cost;
 
   //
   // kinematic of the created pair
   //

   G4double LocalEnerDeposit = 0. ;
   aParticleChange.SetNumberOfSecondaries(2) ; 

   G4double TotalAvailableEnergy = PositKinEnergy + 2*electron_mass_c2;
   G4double Phot1Energy = epsil*TotalAvailableEnergy;

   G4double GammaCut= (G4Gamma::GetCutsInEnergy())[aMaterial->GetIndex()];

   if (Phot1Energy > GammaCut)
      {
        G4ThreeVector Phot1Direction ( dirx, diry, dirz );
        Phot1Direction.rotateUz(PositDirection);   
 
        // create G4DynamicParticle object for the particle1  
        G4DynamicParticle* aParticle1= new G4DynamicParticle (G4Gamma::Gamma(),
                                                 Phot1Direction, Phot1Energy);
        aParticleChange.AddSecondary( aParticle1 ) ; 
       }
    else
       { LocalEnerDeposit += Phot1Energy; }

    G4double Phot2Energy =(1.-epsil)*TotalAvailableEnergy; 

   if (Phot2Energy > GammaCut)
      {
        G4double Eratio = Phot1Energy/Phot2Energy;
        G4double PositP = sqrt(PositKinEnergy*(PositKinEnergy+2.*electron_mass_c2));
        G4ThreeVector Phot2Direction (-dirx*Eratio, -diry*Eratio,
                                         (PositP-dirz*Phot1Energy)/Phot2Energy); 
        Phot2Direction.rotateUz(PositDirection);
 
        // create G4DynamicParticle object for the particle2 
        G4DynamicParticle* aParticle2= new G4DynamicParticle (G4Gamma::Gamma(),
                                                 Phot2Direction, Phot2Energy);
        aParticleChange.AddSecondary( aParticle2 ) ; 
       }
    else
       { LocalEnerDeposit += Phot2Energy; }

   aParticleChange.SetLocalEnergyDeposit( LocalEnerDeposit ) ;

   //
   // Kill the incident positron 
   //

   aParticleChange.SetMomentumChange( 0., 0., 0. ) ;
   aParticleChange.SetEnergyChange( 0. ) ; 
   aParticleChange.SetStatusChange( fStopAndKill ) ;

   return &aParticleChange;
}

 
G4VParticleChange* G4IeplusAnnihilation::AtRestDoIt(const G4Track& aTrack,
                                                  const G4Step&  aStep)
//
// Performs the e+ e- annihilation when both particles are assumed at rest.
// It generates two back to back photons with energy = electron_mass.
// The angular distribution is isotropic. 
// GEANT4 internal units
//
// Note : Effects due to binding of atomic electrons are negliged.
 
{
   aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();

   aParticleChange.SetNumberOfSecondaries(2) ; 

   if (electron_mass_c2 > (G4Gamma::GetCutsInEnergy())[aMaterial->GetIndex()])
      {
        G4double cosTeta = 2*G4UniformRand()-1. , sinTeta = sqrt(1.-cosTeta*cosTeta);
        G4double Phi     = twopi * G4UniformRand() ;
        G4ThreeVector Direction (sinTeta*cos(Phi), sinTeta*sin(Phi), cosTeta);   
 
        aParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                                 Direction, electron_mass_c2) );
        aParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                                -Direction, electron_mass_c2) ); 

        aParticleChange.SetLocalEnergyDeposit(0.);
       }
   else
      { aParticleChange.SetLocalEnergyDeposit( 2*electron_mass_c2 ); }

   // Kill the incident positron 
   //
   aParticleChange.SetStatusChange( fStopAndKill );
      
   return &aParticleChange;
}

void G4IeplusAnnihilation::PrintInfoDefinition()
{
 G4String comments="Total cross section from Heitler formula (2 photon annihilation).\n";
          comments += "        gamma energies sampled according Heitler";

  G4cout << endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestEnergyLimit,"Energy")
         << " to " << G4BestUnit(HighestEnergyLimit,"Energy")
         << " in " << NumbBinTable << " bins. \n";
}


