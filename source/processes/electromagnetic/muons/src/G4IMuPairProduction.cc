// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IMuPairProduction.cc,v 1.1 1999-01-07 16:11:07 gunter Exp $
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
//      -------- G4IMuPairProduction physics process ---------
//                by Laszlo Urban, May 1998 
// **************************************************************
// 04-06-98, in DoIt, secondary production condition: range>min(threshold,safety)
// --------------------------------------------------------------

#include "G4IMuPairProduction.hh"
#include "G4EnergyLossTables.hh"

// static members ........
G4PhysicsTable* G4IMuPairProduction::themuplusLambdaTable = NULL ;
G4PhysicsTable* G4IMuPairProduction::themuminusLambdaTable = NULL ;

G4int G4IMuPairProduction::nzdat = 5 ;
G4double G4IMuPairProduction::zdat[]={1.,4.,13.,29.,92.};
G4int G4IMuPairProduction::ntdat = 8 ;
G4double G4IMuPairProduction::tdat[]={1.e3,1.e4,1.e5,1.e6,1.e7,1.e8,1.e9,1.e10};
G4int G4IMuPairProduction::NBIN = 100 ; //500 ;
G4double G4IMuPairProduction::ya[1000]={0.};
G4double G4IMuPairProduction::proba[5][8][1000]={0.};
 
// constructor
 
G4IMuPairProduction::G4IMuPairProduction(const G4String& processName)
  : G4IMuEnergyLoss("IMuPairProduction"),      // initialization
    LowestKineticEnergy (1.*keV),
    HighestKineticEnergy (1000000.*TeV),
    TotBin(200),
    MinKineticEnergy (1.*GeV),
    MinCutValue (1.*keV),
     NumberOfBuildPhysicsTableCalls(0),
    theElectron (G4Electron::Electron() ),
    thePositron (G4Positron::Positron() ),
    theMuonMinus ( G4MuonMinus::MuonMinus() ),
    theMuonPlus ( G4MuonPlus::MuonPlus() )
{

    theMeanFreePathTable = NULL ;

  cout << endl ;
  cout << "****************************************************************" << endl;
  cout << "****************************************************************" << endl;
  cout << "**                                                            **" << endl; 
  cout << "**   G4IMuPairProduction :                                     **" << endl;
  cout << "**      cross section/energy loss is calculated using         **" << endl;
  cout << "**      the accurate cross section formula of R.Kokoulin,     **" << endl;
  cout << "**      ( building the tables takes a lot of time,            **" << endl;
  cout << "**                   PLEASE BE PATIENT  ! )                   **" << endl;
  cout << "**      sampling of the energy of the pair is generated       **" << endl;
  cout << "**      according to the sampling tables.                     **" << endl;
  cout << "**                                                            **" << endl;
  cout << "****************************************************************" << endl;
  cout << "****************************************************************" << endl;
}
 
// destructor
 
G4IMuPairProduction::~G4IMuPairProduction()
{
   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }

   if (&PartialSumSigma) {
      PartialSumSigma.clearAndDestroy();
   }
}
 
 
// methods.............................................................................
 
void G4IMuPairProduction::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
//  just call BuildLossTable+BuildLambdaTable
{
 NumberOfBuildPhysicsTableCalls += 1 ;
 if(NumberOfBuildPhysicsTableCalls == 1)
 {
    ParticleMass = aParticleType.GetPDGMass() ;

   // make tables ........
 
    BuildLossTable(aParticleType) ;
 
  if(&aParticleType==theMuonMinus)
  {
    RecorderOfmuminusProcess[CounterOfmuminusProcess] = (*this).theLossTable ;
    CounterOfmuminusProcess++;
  }
  else
  {
    RecorderOfmuplusProcess[CounterOfmuplusProcess] = (*this).theLossTable ;
    CounterOfmuplusProcess++;
  }

    BuildLambdaTable(aParticleType) ;

    MakeSamplingTables(&aParticleType) ;

    G4IMuEnergyLoss::BuildDEDXTable(aParticleType) ;
 }
 else
 {

    BuildNlambdaTable(aParticleType) ;

    BuildCoeffATable(aParticleType) ;
    BuildCoeffBTable(aParticleType) ;
    BuildCoeffCTable(aParticleType) ;

    BuildInverseNlambdaTable(aParticleType) ;

    G4int printflag = 0 ;
    if(printflag>0)
      TestOfInversion(aParticleType,printflag) ;

    NumberOfBuildPhysicsTableCalls = 0 ;
 }

}

void G4IMuPairProduction::TestOfInversion(
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
    for (G4int i=0; i<TotBin-2; i++)
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
    G4cout << "G4IMuPairProduction::TestOfInversion (T->Nlambda->Tprime) " << endl ;
    G4cout << "particle= " << aParticleType.GetParticleName() <<
            "   material= " << (*theMaterialTable)[J]->GetName() << endl ;
    G4cout << "max (Tprime-T)/T in % =" << setw(10) << setprecision(3) << delta
;
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

void G4IMuPairProduction::BuildNlambdaTable(
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
void G4IMuPairProduction::BuildNlambdaVector(
                               const G4ParticleDefinition& aParticleType,
                                       G4int materialIndex,
                                       G4PhysicsLogVector* nlambdaVector)
{
  G4double LowEdgeEnergy,T,Tlast,dEdx,Value,Vlast,u,du,coeff ;
  G4double thresholdEnergy ;
  const G4int nbin = 100 ;
  G4bool isOut ;
  const G4double small = 1.e-10;
  const G4double plowloss = 0.5 ;  //this should be a data member of en.loss!

  const G4MaterialTable* theMaterialTable=
                          G4Material::GetMaterialTable();

  //  it is a lower limit for the threshold energy ....?????
    thresholdEnergy =  MinKineticEnergy ;  

  // here assumed that the threshold energy for the process >=
  //                LowestKineticEnergy (temporarily )
  if(thresholdEnergy >= LowestKineticEnergy)
  {
    Value = 0. ;
  }
  else
  {
  // extrapolation needed ..................
  //first integral from thresholdEnergy to LowestKineticEnergy
     Value = 0. ;
  }

  nlambdaVector->PutValue(0,Value) ;
  Tlast = LowestKineticEnergy ;
  Vlast = Value ;

  // loop for kinetic energy
  for (G4int i=1; i<TotBin; i++)
  {
    LowEdgeEnergy = nlambdaVector->GetLowEdgeEnergy(i) ;

    Value = 0. ;

   if(LowEdgeEnergy > thresholdEnergy)
   {
    u = log(LowEdgeEnergy/Tlast) ;
    du = u/nbin ;

    u = -du ;
    for(G4int n=0; n<=nbin; n++)
    {
      u += du ;
      T = Tlast*exp(u) ;

      if((n==0)||(n==nbin))
       coeff=0.5 ;
      else
       coeff=1.0 ;

      Value += coeff*T/(G4EnergyLossTables::GetPreciseDEDX(&aParticleType,
                                   T,(*theMaterialTable)[materialIndex])*
                      (*theMeanFreePathTable)[materialIndex]->GetValue(T,isOut))
;
    }

    Value *= du ;
    Value += Vlast ;
    if(Value<small)
       Value = 0. ;
   }

    nlambdaVector->PutValue(i,Value) ;

    Tlast = LowEdgeEnergy ;
    Vlast = Value ;

  }

}
void G4IMuPairProduction::BuildCoeffATable(
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
        Value = 0. ;

     aVector->PutValue(i,Value);

     Ti = RTable*Ti ;

   }

   theCoeffATable->insert(aVector);

  }

}
void G4IMuPairProduction::BuildCoeffBTable(
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
        Value = RTable*(Ri-Rim)/((RTable-1.)*Ti) ;

     aVector->PutValue(i,Value);

     Ti = RTable*Ti ;

   }

   theCoeffBTable->insert(aVector);

  }

}
void G4IMuPairProduction::BuildCoeffCTable(
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
        Value = (-Ri+RTable*Rim)/(RTable-1.) ;


     aVector->PutValue(i,Value);

     Ti = RTable*Ti ;

   }

   theCoeffCTable->insert(aVector);

  }

}
void G4IMuPairProduction::BuildInverseNlambdaTable(
                             const G4ParticleDefinition& aParticleType)
{
    G4double T,Smallest,Biggest ;
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
      } while ((Smallest <= small)&&(T<HighestKineticEnergy)) ;

      Biggest = (*theNlambdaTable)[J]->
                       GetValue(HighestKineticEnergy,isOut) ;

     // inverse can be built for "meaningful" cut value only!
      if(Smallest >= Biggest)
      {
         G4cout << endl ;
         G4Exception(
        "Cut value is too big , smaller value should be used !");
      }

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

void G4IMuPairProduction::InvertNlambdaVector(
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
    KineticEnergy = HighestKineticEnergy/RTable ;
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


void G4IMuPairProduction::BuildLossTable(const G4ParticleDefinition& aParticleType)
//  build table for energy loss due to soft pairs
//  tables are built for MATERIALS
//                       *********
{
  G4double KineticEnergy,TotalEnergy,pairloss,Z,
           loss,natom,eCut,pCut ;


   const G4MaterialTable* theMaterialTable =
                                G4Material::GetMaterialTable();

  const G4double SmallLoss = DBL_MIN ;

  ParticleMass = aParticleType.GetPDGMass() ;
  ElectronCutInKineticEnergy = (*theElectron).GetEnergyCuts() ;
  PositronCutInKineticEnergy = (*thePositron).GetEnergyCuts() ;

//  create table

  G4int numOfMaterials = theMaterialTable->length() ;

     if (theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable;
     }
  theLossTable = new G4PhysicsTable(numOfMaterials) ;

//  loop for materials

  for (G4int J=0; J<numOfMaterials; J++)
  {

  // create physics vector and fill it

  G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                               LowestKineticEnergy,HighestKineticEnergy,TotBin);

  // get e-/e+ cut in kin.energy for the material

  ElectronCutInKineticEnergyNow = ElectronCutInKineticEnergy[J] ;
  PositronCutInKineticEnergyNow = PositronCutInKineticEnergy[J] ;

  // get elements in the material

  const G4Material* material = (*theMaterialTable)[J] ;
 
  const G4ElementVector* theElementVector =
          material->GetElementVector() ;
  const G4double* theAtomicNumDensityVector =
          material->GetAtomicNumDensityVector() ;
  const G4int NumberOfElements =
          material->GetNumberOfElements() ;

//  loop for the kinetic energy values

    for (G4int i=0; i<TotBin; i++)
    {
      KineticEnergy = aVector->GetLowEdgeEnergy(i) ;

      TotalEnergy = KineticEnergy+ParticleMass ;

      eCut = ElectronCutInKineticEnergyNow ;
      pCut = PositronCutInKineticEnergyNow ;

      if(eCut>KineticEnergy)
        eCut = KineticEnergy ;
      if(pCut>KineticEnergy)
        pCut = KineticEnergy ;

      pairloss = 0.;

   //  loop for elements in the material

        for (G4int iel=0; iel<NumberOfElements; iel++)
        {

          Z=(*theElementVector)(iel)->GetZ();
          natom = theAtomicNumDensityVector[iel] ;

          loss = ComputePairLoss(&aParticleType,
                                 Z,KineticEnergy,eCut,pCut) ;   

          pairloss += natom*loss ;
        }
  
      if(pairloss<0.)
        pairloss = SmallLoss ;

      aVector->PutValue(i,pairloss);
    }

    theLossTable->insert(aVector);
  }
}


G4double G4IMuPairProduction::ComputePairLoss(
                                     const G4ParticleDefinition* ParticleType,
                                             G4double AtomicNumber,
                                             G4double KineticEnergy,
                                             G4double ElectronEnergyCut, 
                                             G4double PositronEnergyCut) 

// compute loss due to soft pairs  
{
  G4double sqrte = sqrt(exp(1.)) ;
  G4double z13 = exp(log(AtomicNumber)/3.) ;

  G4double loss = 0.0 ;

  if ( AtomicNumber < 1. ) return loss;
  if ( KineticEnergy < MinKineticEnergy ) return loss; 

  G4double CutInPairEnergy = ElectronEnergyCut+PositronEnergyCut
                            +2.*electron_mass_c2 ;
  G4double MinPairEnergy = 4.*electron_mass_c2 ;
  if( CutInPairEnergy <= MinPairEnergy ) return loss ;

  G4double MaxPairEnergy = KineticEnergy+ParticleMass*(1.-0.75*sqrte*z13) ;
  if( CutInPairEnergy >= MaxPairEnergy ) 
      CutInPairEnergy = MaxPairEnergy ;

  G4double aaa,bbb,hhh,x,ep ;
  G4int kkk ;
 // calculate the rectricted loss    
 // numerical integration in log(PairEnergy)
  aaa = log(MinPairEnergy) ;
  bbb = log(CutInPairEnergy) ;
  kkk = int((bbb-aaa)/0.05)+1 ;
  hhh = (bbb-aaa)/kkk ;
 
  for (G4int l=0 ; l<kkk; l++)
  {
    x = aaa+hhh*(l+0.5) ;
    ep = exp(x) ;
    loss += ep*ep*hhh*ComputeDMicroscopicCrossSection(ParticleType,
                                                 KineticEnergy,AtomicNumber,
                                                 ep) ;
  }

  if (loss < 0.) loss = 0.;

  return loss ;
}
 

void G4IMuPairProduction::BuildLambdaTable(const G4ParticleDefinition& ParticleType)

// build  mean free path tables for the e+e- emission by muons.
// tables are build for MATERIALS. 
{
   G4double LowEdgeEnergy , Value;
   G4double FixedEnergy = (LowestKineticEnergy + HighestKineticEnergy)/2. ;

   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

   //create table
   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }
   theMeanFreePathTable = new G4PhysicsTable( G4Material::GetNumberOfMaterials() ) ;
   if(&ParticleType == theMuonPlus )
     themuplusLambdaTable = theMeanFreePathTable ;
   if(&ParticleType == theMuonMinus )
     themuminusLambdaTable = theMeanFreePathTable ;

   G4PhysicsLogVector* ptrVector;
   for ( G4int J=0 ; J < G4Material::GetNumberOfMaterials(); J++ )  
       { 
        //create physics vector then fill it ....
        ptrVector = new G4PhysicsLogVector(LowestKineticEnergy, HighestKineticEnergy,
                                           TotBin ) ;

        const G4Material* material= (*theMaterialTable)[J];

        for ( G4int i = 0 ; i < TotBin ; i++ )      
           {
             LowEdgeEnergy = ptrVector->GetLowEdgeEnergy( i ) ;
             Value = ComputeMeanFreePath( &ParticleType, LowEdgeEnergy,
                                         material );  
             ptrVector->PutValue( i , Value ) ;
           }

        theMeanFreePathTable->insertAt( J , ptrVector );

        // Compute the PartialSumSigma table at a given fixed energy
        ComputePartialSumSigma( &ParticleType, FixedEnergy, material) ;       
   }
}

void G4IMuPairProduction::ComputePartialSumSigma(const G4ParticleDefinition* ParticleType,
                                               G4double KineticEnergy,
                                               const G4Material* aMaterial)

// build the table of cross section per element. The table is built for MATERIALS.
// This table is used by DoIt to select randomly an element in the material. 
{
   G4int Imate = aMaterial->GetIndex();
   G4int NbOfElements = aMaterial->GetNumberOfElements();
   const G4ElementVector* theElementVector = aMaterial->GetElementVector(); 
   const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();
   G4double ElectronEnergyCut = (G4Electron::GetCutsInEnergy())[Imate];
   G4double PositronEnergyCut = (G4Positron::GetCutsInEnergy())[Imate];

   PartialSumSigma(Imate) = new G4ValVector(NbOfElements);

   G4double SIGMA = 0. ;

   for ( G4int Ielem=0 ; Ielem < NbOfElements ; Ielem++ )
      {             
        SIGMA += theAtomNumDensityVector[Ielem] * 
                 ComputeMicroscopicCrossSection( ParticleType, KineticEnergy,
                                                 (*theElementVector)(Ielem)->GetZ(), 
                                                 ElectronEnergyCut,PositronEnergyCut );

        PartialSumSigma(Imate)->insertAt(Ielem, SIGMA);
   }
}

G4double G4IMuPairProduction::ComputeMicroscopicCrossSection(
                                     const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy, G4double AtomicNumber,
                                           G4double ElectronEnergyCut,
                                           G4double PositronEnergyCut)
 
// Calculates the microscopic (TOTAL)cross section in GEANT4 internal units.
//
 
{

  G4double sqrte = sqrt(exp(1.)) ;
  G4double z13 = exp(log(AtomicNumber)/3.) ;

  G4double CrossSection = 0.0 ;

  if ( AtomicNumber < 1. ) return CrossSection;
  if ( KineticEnergy < MinKineticEnergy ) return CrossSection; 

  G4double CutInPairEnergy = ElectronEnergyCut+PositronEnergyCut
                            +2.*electron_mass_c2 ;

  if( CutInPairEnergy < 4.*electron_mass_c2 )
    CutInPairEnergy = 4.*electron_mass_c2 ;

  G4double MaxPairEnergy = KineticEnergy+ParticleMass*(1.-0.75*sqrte*z13) ;
  if( CutInPairEnergy >= MaxPairEnergy ) return CrossSection ;

  G4double aaa,bbb,hhh,x,ep ;
  G4int kkk ;
 // calculate the total cross section
 // numerical integration in log(PairEnergy)
  aaa = log(CutInPairEnergy) ;
  bbb = log(MaxPairEnergy) ;
  kkk = int((bbb-aaa)/0.05)+1 ;
  hhh = (bbb-aaa)/kkk ;
 
  for (G4int l=0 ; l<kkk; l++)
  {
    x = aaa+hhh*(l+0.5) ;
    ep = exp(x) ;
    CrossSection += ep*hhh*ComputeDMicroscopicCrossSection(ParticleType,
                                                 KineticEnergy,AtomicNumber,
                                                 ep) ;
  }

  if (CrossSection < 0.) CrossSection = 0.;

  return CrossSection;
                          
}

void G4IMuPairProduction::MakeSamplingTables(const G4ParticleDefinition* ParticleType)
{
 
  G4double epbin[1000],xbin[1000],prbin[1000] ;
  G4int nbin;
  G4double AtomicNumber,KineticEnergy,MinPairEnergy ;  
  G4double c,y,ymin,ymax,dy,yy,dx,x,ep ;

  MinPairEnergy = 4.*electron_mass_c2 ;

  static G4int probtable=0 ;
               
  if(probtable == 0 )
  {
   probtable = 1 ;
   G4double sqrte = sqrt(exp(1.)) ;

   for (G4int iz=0; iz<nzdat; iz++)
   {
    AtomicNumber = zdat[iz];
    G4double z13 = exp(log(AtomicNumber)/3.) ;

    for (G4int it=0; it<ntdat; it++)
    {
     KineticEnergy = tdat[it];
     G4double MaxPairEnergy = KineticEnergy+ParticleMass*(1.-0.75*sqrte*z13) ;

     G4double CrossSection = 0.0 ;

     G4int NbofIntervals ;
    // calculate the differential cross section
    // numerical integration in    
    //  log(PairEnergy/MinPairEnergy)/log(MaxPairEnergy/MinPairEnergy)
    c = log(MaxPairEnergy/MinPairEnergy) ;
    ymin = -5. ;
    ymax = 0. ;
    dy = (ymax-ymin)/NBIN ; 

    nbin=-1;              

    y = ymin - 0.5*dy ;
    yy = ymin - dy ;
    for (G4int i=0 ; i<NBIN; i++)
    {
      y += dy ;
      x = exp(y) ;
      yy += dy ;
      dx = exp(yy+dy)-exp(yy) ;
      
      ep = MinPairEnergy*exp(c*x) ;

     CrossSection += ep*dx*ComputeDMicroscopicCrossSection(ParticleType,
                                                 KineticEnergy,AtomicNumber,
                                                 ep) ;
     if(nbin<NBIN)
     {
      nbin += 1 ;
      epbin[nbin]=ep;
      xbin[nbin]=x;
      prbin[nbin]=CrossSection ;
      ya[nbin]=y ;
      proba[iz][it][nbin] = CrossSection ;
     }
    }

    if(CrossSection > 0.)
    {
  
     G4double e,x,logx,prob,fx ;
     G4double plast = 0.;
     G4double xlast= 0.;
  
     for(G4int ib=0; ib<=nbin; ib++)
     {
       prbin[ib] /= CrossSection ;
       proba[iz][it][ib] /= CrossSection ;

       prob = prbin[ib] ;
       e=epbin[ib] ;
       x = xbin[ib] ;
       logx = log(x) ;
       fx=(prob-plast)/(x-xlast);
       if(ib == 0)
         fx /= 2. ;
       plast = prob ;
       xlast = x;
     }
    }
   }
  }
 
 }

} 
 

G4double G4IMuPairProduction::ComputeDDMicroscopicCrossSection(
                                 const G4ParticleDefinition* ParticleType,
                                 G4double KineticEnergy, G4double AtomicNumber,
                                 G4double PairEnergy,G4double asymmetry)
 // Calculates the double differential (DD) microscopic cross section 
 //   using the cross section formula of R.P. Kokoulin (18/01/98)
{
  G4double sqrte = sqrt(exp(1.)) ;

  G4double bbbtf= 183. ;
  G4double bbbh = 202.4 ; 
  G4double g1tf = 1.95e-5 ;
  G4double g2tf = 5.3e-5 ;
  G4double g1h  = 4.4e-5 ;
  G4double g2h  = 4.8e-5 ;

  G4double massratio = ParticleMass/electron_mass_c2 ;
  G4double massratio2 = massratio*massratio ;
  G4double TotalEnergy = KineticEnergy + ParticleMass ;
  G4double z13 = exp(log(AtomicNumber)/3.) ;
  G4double z23 = z13*z13 ;
  G4double EnergyLoss = TotalEnergy - PairEnergy ;

  G4double c3 = 3.*sqrte*ParticleMass/4. ;

  G4double DDCrossSection = 0. ;
 
  if(EnergyLoss <= c3*z13)
    return DDCrossSection ;

  G4double c7 = 4.*electron_mass_c2 ;
  G4double c8 = 6.*ParticleMass*ParticleMass ;
  G4double alf = c7/PairEnergy ;
  G4double a3 = 1. - alf ;

  if(a3 <= 0.)
    return DDCrossSection ;

 // zeta calculation
  G4double bbb,g1,g2,zeta1,zeta2,zeta,z2 ;
  if( AtomicNumber < 1.5 )
  {
    bbb = bbbh ;
    g1  = g1h ;
    g2  = g2h ;
  }
  else
  {
    bbb = bbbtf ;
    g1  = g1tf ;
    g2  = g2tf ;
  }
  zeta1 = 0.073 * log(TotalEnergy/(ParticleMass+g1*z23*TotalEnergy))-0.26 ;
  if( zeta1 > 0.)
  {
    zeta2 = 0.058*log(TotalEnergy/(ParticleMass+g2*z13*TotalEnergy))-0.14 ;
    zeta  = zeta1/zeta2 ;
  }
  else
  {
    zeta = 0. ;
  }

  z2 = AtomicNumber*(AtomicNumber+zeta) ;

  G4double screen0 = 2.*electron_mass_c2*sqrte*bbb/(z13*PairEnergy) ;
  G4double a0 = TotalEnergy*EnergyLoss ;
  G4double a1 = PairEnergy*PairEnergy/a0 ;
  G4double bet = 0.5*a1 ;
  G4double xi0 = 0.25*massratio2*a1 ;
  G4double del = c8/a0 ; 
 // G4double tmn = log((alf+2.*del*a3)/(1.+(1.-del)*sqrt(a3))) ;
 
  G4double romin = 0. ;
  G4double romax = (1.-del)*sqrt(1.-c7/PairEnergy) ;

  if((asymmetry < romin) || (asymmetry > romax))
    return DDCrossSection ;

  G4double a4 = 1.-asymmetry ;
  G4double a5 = a4*(2.-a4) ;
  G4double a6 = 1.-a5 ;
  G4double a7 = 1.+a6 ;
  G4double a9 = 3.+a6 ;
  G4double xi = xi0*a5 ;
  G4double xii = 1./xi ;
  G4double xi1 = 1.+xi ;
  G4double screen = screen0*xi1/a5 ;
   
  G4double yeu = 5.-a6+4.*bet*a7 ;
  G4double yed = 2.*(1.+3.*bet)*log(3.+xii)-a6-a1*(2.-a6) ;
  G4double yel = 1.+yeu/yed ;
  G4double ale=log(bbb/z13*sqrt(xi1*yel)/(1.+screen*yel)) ;
  G4double cre = 0.5*log(1.+2.25/(massratio2*z23)*xi1*yel) ;
  G4double be ;
  if(xi <= 1.e3)
    be = ((2.+a6)*(1.+bet)+xi*a9)*log(1.+xii)+(a5-bet)/xi1-a9;
  else
    be = (3.-a6+a1*a7)/(2.+xi) ;
  G4double fe = (ale-cre)*be ;
  if( fe < 0.)
    fe = 0. ;

  G4double ymu = 4.+a6 +3.*bet*a7 ;
  G4double ymd = a7*(1.5+a1)*log(3.+xi)+1.-1.5*a6 ;
  G4double ym1 = 1.+ymu/ymd ;
  G4double alm_crm = log(bbb*massratio/(1.5*z23*(1.+screen*ym1))) ;
  G4double a10,bm ;
  if( xi >= 1.e-3)
  {
    a10 = (1.+a1)*a5 ;
    bm  = (a7*(1.+1.5*bet)-a10*xii)*log(xi1)+xi*(a5-bet)/xi1+a10 ;
  }
  else
    bm = (5.-a6+bet*a9)*(xi/2.) ;
  G4double fm = alm_crm*bm ;
  if( fm < 0.)
    fm = 0. ;

  DDCrossSection = (fe+fm/massratio2) ;

  DDCrossSection *= 4.*fine_structure_const*fine_structure_const
                   *classic_electr_radius*classic_electr_radius/(3.*pi) ;
  
  DDCrossSection *= z2*EnergyLoss/(TotalEnergy*PairEnergy) ;
 
 
  return DDCrossSection ;

}
      
G4double G4IMuPairProduction::ComputeDMicroscopicCrossSection(
                                 const G4ParticleDefinition* ParticleType,
                                 G4double KineticEnergy, G4double AtomicNumber,
                                 G4double PairEnergy)
 // Calculates the double differential (D) microscopic cross section 
 //   using the cross section formula of R.P. Kokoulin (18/01/98)
{

  static const G4double
  xgi[] ={ 0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801 };

  static const G4double
  wgi[] ={ 0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506 };

  G4double TotalEnergy = KineticEnergy + ParticleMass ;
  G4double EnergyLoss = TotalEnergy - PairEnergy ;
  //G4double romax = (1.-6.*ParticleMass*ParticleMass/(TotalEnergy*EnergyLoss))
  //                *sqrt(1.-4.*electron_mass_c2/PairEnergy) ;
  //G4double tmn = log(1.-romax) ;
  G4double a = 6.*ParticleMass*ParticleMass/(TotalEnergy*EnergyLoss) ;
  G4double b = 4.*electron_mass_c2/PairEnergy ;
  G4double tmn=log((b+2.*a*(1.-b))/(1.+(1.-a)*sqrt(1.-b))) ;

  G4double DCrossSection = 0. ;
  G4double ro ;
// Gaussian integration in ln(1-ro) ( with 8 points)
  for (G4int i=0; i<7; i++)
  {
    ro = 1.-exp(tmn*xgi[i]) ;
    
    DCrossSection += (1.-ro)*ComputeDDMicroscopicCrossSection(
                                                 ParticleType,KineticEnergy,
                                                 AtomicNumber,PairEnergy,ro)
                            *wgi[i] ;
  }

  DCrossSection *= -tmn ;

  return DCrossSection ;

}
      



G4VParticleChange* G4IMuPairProduction::PostStepDoIt(const G4Track& trackData,
                                                  const G4Step& stepData)
{
   static const G4double esq = sqrt(exp(1.));

   static const G4double MassRatio=ParticleMass/electron_mass_c2 ;

   aParticleChange.Initialize(trackData);
   G4Material* aMaterial=trackData.GetMaterial() ;

   const G4DynamicParticle* aDynamicParticle=trackData.GetDynamicParticle();  

   G4double           KineticEnergy     = aDynamicParticle->GetKineticEnergy();
   G4ParticleMomentum ParticleDirection = aDynamicParticle->GetMomentumDirection();

   // e-e+ cut in this material
   G4double ElectronEnergyCut = (G4Electron::GetCutsInEnergy())[aMaterial->GetIndex()];
   G4double PositronEnergyCut = (G4Electron::GetCutsInEnergy())[aMaterial->GetIndex()];
   G4double CutInPairEnergy = ElectronEnergyCut + PositronEnergyCut ;
   G4double MinPairEnergy = 4.*electron_mass_c2 ;
   if (CutInPairEnergy < MinPairEnergy) CutInPairEnergy = MinPairEnergy ;

   // check against insufficient energy
    if (KineticEnergy < CutInPairEnergy )
       {
         aParticleChange.SetMomentumChange( ParticleDirection );
         aParticleChange.SetEnergyChange( KineticEnergy );
         aParticleChange.SetLocalEnergyDeposit (0.); 
         aParticleChange.SetNumberOfSecondaries(0);
         return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
       }

   // select randomly one element constituing the material  
   G4Element* anElement = SelectRandomAtom(aMaterial);

   // limits of the energy sampling
   G4double TotalEnergy = KineticEnergy + ParticleMass ;
   G4double TotalMomentum = sqrt(KineticEnergy*(TotalEnergy+ParticleMass)) ;
   G4double Z3 = anElement->GetIonisation()->GetZ3() ;
   G4double MaxPairEnergy = TotalEnergy-0.75*esq*ParticleMass*Z3 ;

   if(MinPairEnergy >= MaxPairEnergy)
       {
         aParticleChange.SetMomentumChange( ParticleDirection );
         aParticleChange.SetEnergyChange( KineticEnergy );
         aParticleChange.SetLocalEnergyDeposit (0.); 
         aParticleChange.SetNumberOfSecondaries(0);
         return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
       }

   // sample e-e+ energy, pair energy first
   //  sampling using tables (it have to be accelerated)
   G4double PairEnergy,xc,x,yc,y ;
   G4int iZ,iT,iy ;

   // select sampling table ;
   G4double lnZ = log(anElement->GetZ()) ;
   G4double delmin = 1.e10 ;
   G4double del ;
   G4int izz,itt,NBINminus1 ;
   NBINminus1 = NBIN-1 ;
   for (G4int iz=0; iz<nzdat; iz++)
   {
     del = abs(lnZ-log(zdat[iz])) ;
     if(del<delmin)
     {
        delmin=del ;
        izz=iz ;
     }
   }
   delmin = 1.e10 ;
   for (G4int it=0; it<ntdat; it++)
   {
     del = abs(log(KineticEnergy)-log(tdat[it])) ;
     if(del<delmin)
     {
       del=delmin;
       itt=it ;
     }
   }

   xc = log(CutInPairEnergy/MinPairEnergy)/log(MaxPairEnergy/MinPairEnergy) ;
   yc = log(xc) ;
   
   iy = -1 ;
   do {
       iy += 1 ;
      } while ((ya[iy] < yc )&&(iy < NBINminus1)) ;
   G4double norm = 1./(1.-proba[izz][itt][iy]) ;

   G4double r = G4UniformRand() ;
 
   iy = -1 ;
   do {
        iy += 1 ;
      } while (((norm*proba[izz][itt][iy]) < r)&&(iy < NBINminus1)) ;

   //sampling is uniformly in y in the bin
   if( iy < NBINminus1 )
     y = ya[iy] + G4UniformRand() * ( ya[iy+1] - ya[iy]) ;
   else
     y = ya[iy] ;

   x = exp(y) ;

   PairEnergy = MinPairEnergy*exp(x*log(MaxPairEnergy/MinPairEnergy)) ;

  // sample r=(E+-E-)/PairEnergy  ( uniformly .....)
   G4double rmax = (1.-6.*ParticleMass*ParticleMass/(TotalEnergy*
                (TotalEnergy-PairEnergy)))
          *sqrt(1.-MinPairEnergy/PairEnergy) ;
   r = rmax * (-1.+2.*G4UniformRand()) ;

  // compute energies from PairEnergy,r
   G4double ElectronEnergy=(1.-r)*PairEnergy/2. ;  
   G4double PositronEnergy=(1.+r)*PairEnergy/2. ;
     
    //  angles of the emitted particles ( Z - axis along the parent particle)
   //  mean theta for the moment *********************************
   G4double Teta = electron_mass_c2/TotalEnergy ;
  //************************************************************

   G4double Phi  = twopi * G4UniformRand() ;
   G4double dirx = sin(Teta)*cos(Phi) , diry = sin(Teta)*sin(Phi) , dirz = cos(Teta) ;

   G4double LocalEnerDeposit = 0. ;
   G4int numberofsecondaries = 1 ;
   G4int flagelectron = 0 ;
   G4int flagpositron = 1 ; 
   G4DynamicParticle *aParticle1,*aParticle2 ;
   G4double ElectronMomentum , PositronMomentum ;
   G4double finalPx,finalPy,finalPz ;

   G4double ElectKineEnergy = ElectronEnergy - electron_mass_c2 ;

   if (G4EnergyLossTables::GetRange(G4Electron::Electron(),ElectKineEnergy,aMaterial)
        >= min(G4Electron::GetCuts(),stepData.GetPostStepPoint()->GetSafety()) )         
      {
        numberofsecondaries += 1 ;
        flagelectron = 1 ;
        ElectronMomentum = sqrt(ElectKineEnergy*(ElectronEnergy+electron_mass_c2));
        G4ThreeVector ElectDirection ( dirx, diry, dirz );
        ElectDirection.rotateUz(ParticleDirection);   
 
        // create G4DynamicParticle object for the particle1  
        aParticle1= new G4DynamicParticle (G4Electron::Electron(),
                                                 ElectDirection, ElectKineEnergy);
       }
    else
       { LocalEnerDeposit += ElectKineEnergy ; }

   // the e+ is always created (even with Ekine=0) for further annihilation.

   G4double PositKineEnergy = PositronEnergy - electron_mass_c2 ;
   PositronMomentum = sqrt(PositKineEnergy*(PositronEnergy+electron_mass_c2));

   if (G4EnergyLossTables::GetRange(G4Positron::Positron(),PositKineEnergy,aMaterial)
        < min(G4Positron::GetCuts(),stepData.GetPostStepPoint()->GetSafety()) )         
      {
        LocalEnerDeposit += PositKineEnergy ;
        PositKineEnergy = 0. ;
      }
       G4ThreeVector PositDirection ( -dirx, -diry, dirz );
       PositDirection.rotateUz(ParticleDirection);   
 
      // create G4DynamicParticle object for the particle2 
      aParticle2= new G4DynamicParticle (G4Positron::Positron(),
                                                 PositDirection, PositKineEnergy);

  // fill particle change and update initial particle
   aParticleChange.SetNumberOfSecondaries(numberofsecondaries) ; 
   if(flagelectron==1)
        aParticleChange.AddSecondary( aParticle1 ) ; 
   if(flagpositron==1)       
        aParticleChange.AddSecondary( aParticle2 ) ; 

     G4double NewKinEnergy = KineticEnergy - ElectronEnergy - PositronEnergy ;
     G4double finalMomentum=sqrt(NewKinEnergy*
                         (NewKinEnergy+2.*ParticleMass)) ;

   aParticleChange.SetMomentumChange( ParticleDirection );

   G4double KinEnergyCut = 
           (aDynamicParticle->GetDefinition()->GetEnergyCuts())[aMaterial->GetIndex()];

   if (NewKinEnergy > KinEnergyCut)
      {
       aParticleChange.SetEnergyChange( NewKinEnergy );
      }
   else
      {
       aParticleChange.SetEnergyChange(0.);
       LocalEnerDeposit += NewKinEnergy ;
       aParticleChange.SetStatusChange(fStopButAlive);
      }

   aParticleChange.SetLocalEnergyDeposit( LocalEnerDeposit ) ;


   return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

G4Element* G4IMuPairProduction::SelectRandomAtom(G4Material* aMaterial) const
{
  // select randomly 1 element within the material

  const G4int Index = aMaterial->GetIndex();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  G4double rval = G4UniformRand()*((*PartialSumSigma(Index))(NumberOfElements-1));

  for ( G4int i=0; i < NumberOfElements; i++ )
  {

    if (rval <= (*PartialSumSigma(Index))(i)) return ((*theElementVector)(i));
  }
  cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements, NULL pointer returned." << endl;
  return NULL;
}
