// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IMuBremsstrahlung.cc,v 1.4 2000-04-25 14:19:00 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------- G4IMuBremsstrahlung physics process ---------
//                by Laszlo Urban, September 1997
//
// 08-04-98: remove 'tracking cut' of muon in DoIt, MMa
// --------------------------------------------------------------
#include "G4IMuBremsstrahlung.hh"
 
// constructor
 
G4IMuBremsstrahlung::G4IMuBremsstrahlung(const G4String& processName)
  : G4VIMuEnergyLoss("IMuBremsstrahlung"),      // initialization
    LowestKineticEnergy (1.*keV),
    HighestKineticEnergy (1000000.*TeV),
    TotBin(200),
    MinKineticEnergy (1.*GeV),
    MinCutValue (1.*keV),
    NuclearFormFactor (2),
     NumberOfBuildPhysicsTableCalls(0),
    theGamma (G4Gamma::Gamma() ),
    theMuonMinus ( G4MuonMinus::MuonMinus() ),
    theMuonPlus ( G4MuonPlus::MuonPlus() )
{

    theMeanFreePathTable = NULL ;
}
 
// destructor
 
G4IMuBremsstrahlung::~G4IMuBremsstrahlung()
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
 
void G4IMuBremsstrahlung::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
//  just call BuildLossTable+BuildLambdaTable
{
 NumberOfBuildPhysicsTableCalls += 1 ;
 if(NumberOfBuildPhysicsTableCalls == 1)
 {

    ParticleMass = aParticleType.GetPDGMass() ;

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
  
    G4VIMuEnergyLoss::BuildDEDXTable(aParticleType) ;

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

void G4IMuBremsstrahlung::TestOfInversion(
                             const G4ParticleDefinition& aParticleType,
                                   G4int printflag)
{
  G4double T,Nlambda,Tprime,delta,del,sum,delmean,Tdelta ;
  G4bool isOut ;

  const G4MaterialTable* theMaterialTable=
                          G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length() ;

  G4cout.setf(G4std::ios::scientific, G4std::ios::floatfield) ;
  if(printflag>1)
  {
  G4cout << G4endl;
  G4cout << "  particle=" << aParticleType.GetParticleName() << G4endl;
  G4cout << "----------------------" << G4endl;
  }
  for (G4int J=0; J<numOfMaterials; J++)
  {
   if(printflag>1)
   {
    G4cout << G4endl;
    G4cout << " material = " << (*theMaterialTable)[J]->GetName() << G4endl;
    G4cout << " mat.ind.=" << J << "  T           Nlambda         Tprime"
         << "         (Tprime-T)/T(%)" << G4endl ;
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
      G4cout << G4std::setw(18) << G4std::setprecision(6) << T << "  " <<
              G4std::setw(14)<< G4std::setprecision(6) << Nlambda << "  " <<
      G4std::setw(14) << G4std::setprecision(6) << Tprime << "      " <<
              G4std::setw(12) << G4std::setprecision(3) << del << G4endl;
    }
     }
    }
    if(printflag>0)
    {
    G4cout << G4endl;
    G4cout << "G4IMuBremsstrahlung::TestOfInversion (T->Nlambda->Tprime) " << G4endl ;
    G4cout << "particle= " << aParticleType.GetParticleName() <<
            "   material= " << (*theMaterialTable)[J]->GetName() << G4endl ;
    G4cout << "max (Tprime-T)/T in % =" << G4std::setw(10) << G4std::setprecision(3) << delta
;
    G4cout << "  at a kinetic energy " << G4std::setw(10) << G4std::setprecision(3) <<
            Tdelta/MeV << "  MeV" << G4endl;
    delmean /= sum ;
    G4cout << "mean rel.diff. (Tprime-T)/T=" << G4std::setw(10) <<
            G4std::setprecision(3) << delmean <<
            " % (mean is calculated in abs. value)" << G4endl;
    G4cout << G4endl;
    }
  }

}

void G4IMuBremsstrahlung::BuildNlambdaTable(
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
void G4IMuBremsstrahlung::BuildNlambdaVector(
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

  //   threshold energy ....?????
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
void G4IMuBremsstrahlung::BuildCoeffATable(
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

void G4IMuBremsstrahlung::BuildCoeffBTable(
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
void G4IMuBremsstrahlung::BuildCoeffCTable(
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

void G4IMuBremsstrahlung::BuildInverseNlambdaTable(
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
         G4cout << G4endl ;
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
void G4IMuBremsstrahlung::InvertNlambdaVector(
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


void G4IMuBremsstrahlung::BuildLossTable(const G4ParticleDefinition& aParticleType)
//  Build table for energy loss due to soft brems
//  tables are built for MATERIALS
//                       *********
{
  G4double KineticEnergy,TotalEnergy,bremloss,Z,
           loss,natom,Cut ;


   const G4MaterialTable* theMaterialTable =
                                G4Material::GetMaterialTable();

  const G4double SmallLoss = DBL_MIN ;

  ParticleMass = aParticleType.GetPDGMass();
  GammaCutInKineticEnergy = (*theGamma).GetEnergyCuts() ;

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

  // get photon cut in kin.energy for the material

  GammaCutInKineticEnergyNow = GammaCutInKineticEnergy[J] ;

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
      Cut = GammaCutInKineticEnergyNow ;
      if(Cut>KineticEnergy) Cut = KineticEnergy ;


      //loop for elements in the material
      bremloss = 0.;      
      for (G4int iel=0; iel<NumberOfElements; iel++)
        {
          Z=(*theElementVector)(iel)->GetZ();
          natom = theAtomicNumDensityVector[iel] ;
          loss = ComputeBremLoss(Z,KineticEnergy,Cut) ;   
          bremloss += natom*loss ;
        } 
      if(bremloss<0.) bremloss = SmallLoss;
      
      aVector->PutValue(i,bremloss);
    }
    theLossTable->insert(aVector);
  }
}


G4double G4IMuBremsstrahlung::ComputeBremLoss(G4double AtomicNumber,
                                                        G4double T,G4double Cut)

// compute loss due to soft brems  
{

 static const G4double sigmafactor = 1.*cm2/Avogadro ;
 static const G4double apar1[]={
                       0.15328e+0,-0.28426e-1, 0.91526e-1, 0.50875e-2,-0.56839e-2,
                      -0.28592e-2,-0.22830e-3, 0.10264e-3, 0.83753e-3, 0.11703e-2,
                      -0.21296e-1, 0.65650e-2, 0.14409e-2,-0.52079e-3, 0.95797e-3,
                       0.64648e-3, 0.12752e-4,-0.10587e-3,-0.47799e-4, 0.10457e-3};
 static const G4double apar2[]={
                       0.19564e+0,-0.40379e-1, 0.86428e-1, 0.67609e-2,-0.58381e-2,
                      -0.16173e-2,-0.28879e-3, 0.37766e-4, 0.10055e-2, 0.10146e-2,
                      -0.81069e-2, 0.28804e-2, 0.51924e-2,-0.10308e-3, 0.10440e-3,
                      -0.45960e-4,-0.38156e-5,-0.76242e-4, 0.92117e-4, 0.12112e-3};
 static const G4double cpar1[]={
                       0.44198e-6,-0.23098e-8,-0.17542e-8,-0.61447e-7 };
 static const G4double cpar2[]={
                       0.45488e-6, 0.12797e-7,-0.29005e-8,-0.90486e-7 };
 G4double apar[20],cpar[4];

 if(NuclearFormFactor==1)
 {
   G4int ipar1;
   for (ipar1=0; ipar1<20; ipar1++)
   {  apar[ipar1]=apar1[ipar1] ; }
   for (ipar1=0; ipar1<4; ipar1++)
   {  cpar[ipar1]=cpar1[ipar1] ; }
 }
 else
 {
   G4int ipar2;
   for (ipar2=0; ipar2<20; ipar2++)
   {  apar[ipar2]=apar2[ipar2] ; }
   for (ipar2=0; ipar2<4; ipar2++)
   {  cpar[ipar2]=cpar2[ipar2] ; }
 }

 static const G4double CutLim = 1.*MeV ;

 static const G4double T100TeV = 100.*TeV ;
   
 G4double loss = 0.0 ;

 if ( AtomicNumber < 1. ) return loss ;
 if( T < MinKineticEnergy ) return loss ;

 if( Cut > T )
   Cut = T ;

 if( Cut < MinCutValue )
     Cut = MinCutValue ;

 // extrapolation for Cut < Cutlim , part 1.
 G4double CutSave = 0. ;
 if( Cut < CutLim )
 {   CutSave = Cut ;
     Cut     = CutLim ;
 }            

 //extrapolation for  T > T100TeV , part 1.
 G4double Textr,Cutextr ;
 Textr = 0. ;
 if( T > T100TeV )
 {
   Textr = T ;
   Cutextr = Cut ;
   T = T100TeV ;
   if(Cut > T)
     Cut = T ; 
 }

 G4double u = log(AtomicNumber) ;
 G4double v = log(1.+T/ParticleMass) ;
 G4double w = log(T/Cut) ;

 G4double v2,v3,w2,w3,a,c ;
 v2=v*v ;
 v3=v2*v ;
 w2=w*w ;
 w3=w2*w ;

 a=   apar[0]+apar[1]*v+apar[2]*w+apar[3]*v2+apar[4]*v*w+apar[5]*w2+
      apar[6]*v3+apar[7]*v2*w+apar[8]*v*w2+apar[9]*w3+
   u*(apar[10]+apar[11]*v+apar[12]*w+apar[13]*v2+apar[14]*v*w+apar[15]*w2+
      apar[16]*v3+apar[17]*v2*w+apar[18]*v*w2+apar[19]*w3) ;

 c=cpar[0]+u*(cpar[1]+cpar[2]*u)+cpar[3]*Cut/T ;

 loss = c*(1.-exp(-a*v)) ;

 loss *= sigmafactor*Cut*AtomicNumber*(AtomicNumber+1.)/(1.+0.06*u) ;

 // extrapolation for T > T100TeV , part 2.
 if( Textr > 0. )
 {
   loss *= Cutextr/Cut ;
 }

 // extrapolation for Cut < CutLim , part 2.
 if( CutSave > 0.)
   loss *= CutSave/CutLim ;


 if (loss < 0.) loss = 0.;

  return loss ;
}


void G4IMuBremsstrahlung::BuildLambdaTable(const G4ParticleDefinition& ParticleType)

// Build  mean free path tables for the gamma emission by muons.
// tables are Build for MATERIALS. 
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
   PartialSumSigma.resize(G4Material::GetNumberOfMaterials());

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

void G4IMuBremsstrahlung::ComputePartialSumSigma(const G4ParticleDefinition* ParticleType,
                                               G4double KineticEnergy,
                                               const G4Material* aMaterial)

// Build the table of cross section per element. The table is built for MATERIALS.
// This table is used by DoIt to select randomly an element in the material. 
{
   G4int Imate = aMaterial->GetIndex();
   G4int NbOfElements = aMaterial->GetNumberOfElements();
   const G4ElementVector* theElementVector = aMaterial->GetElementVector(); 
   const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();
   G4double GammaEnergyCut = (G4Gamma::GetCutsInEnergy())[Imate];

   PartialSumSigma(Imate) = new G4ValVector(NbOfElements);

   G4double SIGMA = 0. ;

   for ( G4int Ielem=0 ; Ielem < NbOfElements ; Ielem++ )
      {             
        SIGMA += theAtomNumDensityVector[Ielem] * 
                 ComputeMicroscopicCrossSection( ParticleType, KineticEnergy,
                                                 (*theElementVector)(Ielem)->GetZ(), 
                                                 GammaEnergyCut );
        PartialSumSigma(Imate)->insertAt(Ielem, SIGMA);
   }
}

G4double G4IMuBremsstrahlung::ComputeMicroscopicCrossSection(
                                     const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy, G4double AtomicNumber,
                                           G4double GammaEnergyCut)
 
// Calculates the microscopic cross section in GEANT4 internal units.
// A parametrized formula from L. Urban is used to estimate the total cross section.
//
 
{

 static const G4double sigmafactor = 1.*cm2/Avogadro ;

 static const G4double apar1[]={
                       0.83102e-1,-0.24945e-1, 0.32233e-1, 0.47013e-2,-0.15698e-2,
                       0.18918e-2,-0.21786e-3, 0.50011e-4,-0.35375e-4,-0.84328e-4,
                      -0.35881e-1, 0.13027e-1, 0.29771e-2,-0.14009e-2,-0.74999e-3,
                       0.15088e-3, 0.47484e-4, 0.33027e-4, 0.10146e-4,-0.14978e-4};

 static const G4double apar2[]={
                       0.68931e-1,-0.15069e-1, 0.39220e-1, 0.32425e-2,-0.34125e-2,
                       0.16335e-2,-0.15817e-3, 0.14446e-3, 0.89327e-5,-0.10687e-3,
                      -0.11443e-1, 0.36703e-2, 0.19308e-2,-0.27223e-3,-0.19030e-3,
                       0.10404e-3, 0.51231e-5,-0.22346e-5, 0.13321e-5,-0.36106e-5};
 static const G4double cpar1[]={
                       0.41842e-6,-0.63439e-9,-0.18174e-8};

 static const G4double cpar2[]={
                       0.43455e-6, 0.12501e-7,-0.28730e-8};

 static const G4double bpar1[]={
                       0.15328e+0,-0.28426e-1, 0.91526e-1, 0.50875e-2,-0.56839e-2,
                      -0.28592e-2,-0.22830e-3, 0.10264e-3, 0.83753e-3, 0.11703e-2,
                      -0.21296e-1, 0.65650e-2, 0.14409e-2,-0.52079e-3, 0.95797e-3,
                       0.64648e-3, 0.12752e-4,-0.10587e-3,-0.47799e-4, 0.10457e-3};

 static const G4double bpar2[]={
                       0.19564e+0,-0.40379e-1, 0.86428e-1, 0.67609e-2,-0.58381e-2,
                      -0.16173e-2,-0.28879e-3, 0.37766e-4, 0.10055e-2, 0.10146e-2,
                      -0.81069e-2, 0.28804e-2, 0.51924e-2,-0.10308e-3, 0.10440e-3,
                      -0.45960e-4,-0.38156e-5,-0.76242e-4, 0.92117e-4, 0.12112e-3};

 static const G4double dpar1[]={
                       0.44198e-6,-0.23098e-8,-0.17542e-8,-0.61447e-7 };

 static const G4double dpar2[]={
                       0.45488e-6, 0.12797e-7,-0.29005e-8,-0.90486e-7 };

 G4double apar[20],cpar[3],bpar[20],dpar[4];

 if(NuclearFormFactor==1)
 {
   G4int ipar1;
   for (ipar1=0; ipar1<20; ipar1++)
   {  apar[ipar1]=apar1[ipar1] ; 
      bpar[ipar1]=bpar1[ipar1] ; }
   for (ipar1=0; ipar1<3; ipar1++)
   {  cpar[ipar1]=cpar1[ipar1] ; 
      dpar[ipar1]=dpar1[ipar1] ; } 
      dpar[3]=dpar1[3];  
 }
 else
 {
   G4int ipar2;
   for (ipar2=0; ipar2<20; ipar2++)
   {  apar[ipar2]=apar2[ipar2] ; 
      bpar[ipar2]=bpar2[ipar2] ; }
   for (ipar2=0; ipar2<3; ipar2++)
   {  cpar[ipar2]=cpar2[ipar2] ; 
      dpar[ipar2]=dpar2[ipar2] ; }
      dpar[3]=dpar2[3];  
 }

 static const G4double CutLim = 1.*MeV ;
 static const G4double T100TeV = 100.*TeV ;
  
 G4double CrossSection = 0.0 ;

 if ( AtomicNumber < 1. ) return CrossSection;
 if ( KineticEnergy < MinKineticEnergy ) return CrossSection; 
 if ( KineticEnergy <= GammaEnergyCut ) return CrossSection;

 G4double CutSav = 0. ;
 
 if ( GammaEnergyCut < MinCutValue )
      GammaEnergyCut = MinCutValue ;

 // extrapolation for GammaEnergyCut < CutLim , part 1.
 if( GammaEnergyCut < CutLim)
 {
     CutSav = GammaEnergyCut ;
     GammaEnergyCut = CutLim ;
 }  

 // extrapolation for KineticEnergy > T100TeV , part 1.
 G4double Textr,Cutextr ;
 Textr = 0. ;
 if( KineticEnergy > T100TeV )
 {
   Textr = KineticEnergy ;
   Cutextr = GammaEnergyCut ;
   KineticEnergy = T100TeV ;
 } 

 G4double u = log(AtomicNumber) ;
 G4double v = log(1.+(KineticEnergy-GammaEnergyCut)/ParticleMass) ;
 G4double w = log(KineticEnergy/GammaEnergyCut) ;

 G4double v2,v3,w2,w3,a,c ;
 v2=v*v ;
 v3=v2*v ;
 w2=w*w ;
 w3=w2*w ;

 a=   apar[0]+apar[1]*v+apar[2]*w+apar[3]*v2+apar[4]*v*w+apar[5]*w2+
      apar[6]*v3+apar[7]*v2*w+apar[8]*v*w2+apar[9]*w3+
   u*(apar[10]+apar[11]*v+apar[12]*w+apar[13]*v2+apar[14]*v*w+apar[15]*w2+
      apar[16]*v3+apar[17]*v2*w+apar[18]*v*w2+apar[19]*w3) ;

 c=cpar[0]+u*(cpar[1]+cpar[2]*u) ;

 CrossSection = c*(1.-exp(-a*v)) ;

 CrossSection *= sigmafactor*AtomicNumber*(AtomicNumber+1.)*w/(1.+0.06*u) ;

 //extrapolation for KineticEnergy > T100TeV , part 2.
 if( Textr > 0.)
 {
   CrossSection *= log(Textr/Cutextr)/w ;
 }

 // extrapolation for GammaEnergyCut < CutLim , part 2.
 if( CutSav > 0.)
 {
   G4double loss ;

   v=log(1.+KineticEnergy/ParticleMass) ;
   v2=v*v ;
   v3=v2*v ;

   a=   bpar[0]+bpar[1]*v+bpar[2]*w+bpar[3]*v2+bpar[4]*v*w+bpar[5]*w2+
        bpar[6]*v3+bpar[7]*v2*w+bpar[8]*v*w2+bpar[9]*w3+
     u*(bpar[10]+bpar[11]*v+bpar[12]*w+bpar[13]*v2+bpar[14]*v*w+bpar[15]*w2+
        bpar[16]*v3+bpar[17]*v2*w+bpar[18]*v*w2+bpar[19]*w3) ;

   c=dpar[0]+u*(dpar[1]+dpar[2]*u)+dpar[3]*GammaEnergyCut/KineticEnergy ;

   loss = c*(1.-exp(-a*v)) ;

   loss *= sigmafactor*AtomicNumber*(AtomicNumber+1.)/(1.+0.06*u) ;

   if( GammaEnergyCut < KineticEnergy)
      loss *= GammaEnergyCut ;
   else
      loss *= KineticEnergy ;  


   CrossSection += loss*log(CutLim/CutSav)/CutLim ;

 }

 if (CrossSection < 0.) CrossSection = 0.;


 return CrossSection;
}
 

G4VParticleChange* G4IMuBremsstrahlung::PostStepDoIt(const G4Track& trackData,
                                                  const G4Step& stepData)
//
// 
{
 //nuclear form factor correction fn=2./(3.*Z**(1/3)).............
 //           or
 //nuclear form factor correction fn=exp(-fn1*(fn2*A**(1/3)-fn3).............
   static const G4double fn11=2./3.;
   static const G4double fn1=0.128 ;
   static const G4double fn2=1.18  ;
   static const G4double fn3=0.48  ; 
   
   static const G4double R=189. ;
   static const G4double esq=sqrt(exp(1.)) ;

   static const G4double MassRatio=ParticleMass/electron_mass_c2 ;

   aParticleChange.Initialize(trackData);
   G4Material* aMaterial=trackData.GetMaterial() ;

   const G4DynamicParticle* aDynamicParticle=trackData.GetDynamicParticle();  

   G4double           KineticEnergy     = aDynamicParticle->GetKineticEnergy();
   G4ParticleMomentum ParticleDirection = aDynamicParticle->GetMomentumDirection();

   // Gamma cut in this material
   G4double GammaEnergyCut = (G4Gamma::GetCutsInEnergy())[aMaterial->GetIndex()];

   // check against insufficient energy
    if (KineticEnergy < GammaEnergyCut)
       {
         aParticleChange.SetMomentumChange( ParticleDirection );
         aParticleChange.SetEnergyChange( KineticEnergy );
         aParticleChange.SetLocalEnergyDeposit (0.); 
         aParticleChange.SetNumberOfSecondaries(0);
         return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
       }

   // select randomly one element constituing the material  
   G4Element* anElement = SelectRandomAtom(aMaterial);

   // Extract data for this Element
   G4double Z3 = anElement->GetIonisation()->GetZ3();
   G4double A = (anElement->GetA())*mole/g ;
   G4double A3 = exp(log(A)/3.) ;

   // nuclear form factor correction
   G4double fn ;
   if(NuclearFormFactor==1)
    fn=fn11/Z3 ;
   else
    fn=exp(-fn1*(fn2*A3-fn3)) ;

   // limits of the energy sampling
   G4double TotalEnergy = KineticEnergy + ParticleMass ;
   G4double vmin = GammaEnergyCut/TotalEnergy ; 
   G4double vmax = 1.-0.5*esq*ParticleMass/(fn*TotalEnergy) ;

   if(vmin >= vmax)
       {
         aParticleChange.SetMomentumChange( ParticleDirection );
         aParticleChange.SetEnergyChange( KineticEnergy );
         aParticleChange.SetLocalEnergyDeposit (0.); 
         aParticleChange.SetNumberOfSecondaries(0);
       //  ResetNumberOfInteractionLengthLeft();
       //  return &aParticleChange ;
        return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
       }

   //compute maximum of the rejection function
   G4double qmin=0.5*ParticleMass*vmin/(TotalEnergy*(1.-vmin)) ;
   G4double w1=R/Z3 ;
   G4double w2=MassRatio*esq*w1 ;
   G4double w3=fn*MassRatio ;
   G4double w4=w1/(1.+qmin*w2) ;
   G4double rejmax=log(w3*w4)*(4.*(1.-vmin)/3.+vmin*vmin) ; 

   // sample photon energy
   G4double v,reject ;
       do {
             v = vmin*pow((vmax/vmin), G4UniformRand());
             qmin=0.5*ParticleMass*v/(TotalEnergy*(1.-v)) ;
             w4=w1/(1.+qmin*w2) ;
             reject=log(w3*w4)*(4.*(1.-v)/3.+v*v) ;  
          }  while( reject < G4UniformRand()*rejmax );

   if( v <= 0.)
     return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);


   //  angles of the emitted gamma. ( Z - axis along the parent particle)
   //  Teta = electron_mass_c2/TotalEnergy for the moment .....

   G4double Teta = electron_mass_c2/TotalEnergy ;
   G4double Phi  = twopi * G4UniformRand() ;
   G4double dirx = sin(Teta)*cos(Phi) , diry = sin(Teta)*sin(Phi) , dirz = cos(Teta) ;

   G4ThreeVector GammaDirection ( dirx, diry, dirz);
   GammaDirection.rotateUz(ParticleDirection);   
 
   // create G4DynamicParticle object for the Gamma 
   G4double GammaEnergy = v*TotalEnergy; 
   G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
                                                  GammaDirection, GammaEnergy);

   aParticleChange.SetNumberOfSecondaries(1);
   aParticleChange.AddSecondary(aGamma); 

   //
   // Update the incident particle 
   //

   G4double NewKinEnergy = KineticEnergy - GammaEnergy;
   if (NewKinEnergy > 0.)      
      {
       aParticleChange.SetMomentumChange(ParticleDirection);
       aParticleChange.SetEnergyChange(NewKinEnergy);
       aParticleChange.SetLocalEnergyDeposit (0.); 
      }
   else
      {
       aParticleChange.SetEnergyChange(0.);
       aParticleChange.SetLocalEnergyDeposit (0.);
       aParticleChange.SetStatusChange(fStopButAlive);
      }

   return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
                    
}

G4Element* G4IMuBremsstrahlung::SelectRandomAtom(G4Material* aMaterial) const
{
  // select randomly 1 element within the material

  const G4int Index = aMaterial->GetIndex();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  G4double rval = G4UniformRand()*((*PartialSumSigma(Index))(NumberOfElements-1));
  for ( G4int i=0; i < NumberOfElements; i++ )
    if (rval <= (*PartialSumSigma(Index))(i)) return ((*theElementVector)(i));
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements, NULL pointer returned." << G4endl;
  return NULL;
}
