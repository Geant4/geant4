// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IeBremsstrahlung.cc,v 1.1 1999-01-07 16:11:20 gunter Exp $
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
//      ------------ G4IeBremsstrahlung physics process --------
//                     by Michel Maire, 24 July 1996
// **************************************************************
// ************************************************************
// It is the first implementation of the BREMSSTRAHLUNG
// PROCESS. (  photons   + continuous energy loss)
//   using an INTEGRAL APPROACH instead of the differential
//   one used in the standard implementation .
// ************************************************************
//                by Laszlo Urban, 23 June 1998
// --------------------------------------------------------------
// 28/10/98:LPM effect implemented+ small changes, cleanup  L.Urban


#include "G4IeBremsstrahlung.hh"
#include "G4UnitsTable.hh"
 
 
G4IeBremsstrahlung::G4IeBremsstrahlung(const G4String& processName)
  : G4IeEnergyLoss(processName),     
    theMeanFreePathTable(NULL),
    theNlambdaTable(NULL),
    theInverseNlambdaTable(NULL),
    theCoeffATable(NULL),
    theCoeffBTable(NULL),
    theCoeffCTable(NULL),
    LowestKineticEnergy (10.*keV),
    HighestKineticEnergy (100.*TeV),
    TotBin(100),
    NumberOfBuildPhysicsTableCalls(0),
    theGamma (G4Gamma::Gamma() ),
    theElectron ( G4Electron::Electron() ),
    thePositron ( G4Positron::Positron() )
{  }
 
G4IeBremsstrahlung::~G4IeBremsstrahlung()
{
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
   if (&PartialSumSigma) {
      PartialSumSigma.clearAndDestroy();
   }
}
 
void G4IeBremsstrahlung::SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins)
{
  LowestKineticEnergy = lowE;  HighestKineticEnergy = highE; TotBin = nBins;
}

 
void G4IeBremsstrahlung::BuildPhysicsTable(
                                 const G4ParticleDefinition& aParticleType)
{
  RTable = exp(log(HighestKineticEnergy/LowestKineticEnergy)/TotBin) ;
  NumberOfBuildPhysicsTableCalls += 1 ;
  if(NumberOfBuildPhysicsTableCalls == 1)
  { 
    BuildLossTable(aParticleType) ;
    if(&aParticleType==theElectron)
    {
      RecorderOfElectronProcess[CounterOfElectronProcess] =
                                                      (*this).theLossTable ;
      CounterOfElectronProcess++;
    }
    else
    {
      RecorderOfPositronProcess[CounterOfPositronProcess] =
                                                      (*this).theLossTable ;
      CounterOfPositronProcess++;
    }
    BuildLambdaTable(aParticleType) ;
    BuildDEDXTable(aParticleType) ;
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

    if(&aParticleType==theElectron)
      PrintInfoDefinition();
 }
}

void G4IeBremsstrahlung::TestOfInversion(
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
     G4cout << "G4IeBremsstrahlung::TestOfInversion (T->Nlambda->Tprime) " 
            << endl ;
     G4cout << "particle= " << aParticleType.GetParticleName() <<
            "   material= " << (*theMaterialTable)[J]->GetName() << endl ;
     G4cout << "max (Tprime-T)/T in % =" << setw(10) << setprecision(3)
            << delta ;
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

void G4IeBremsstrahlung::BuildNlambdaTable(
                             const G4ParticleDefinition& aParticleType)
{
  const G4MaterialTable* theMaterialTable=
                          G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length() ;

  if(theNlambdaTable)
  {  theNlambdaTable->clearAndDestroy();
     delete theNlambdaTable ;           }

  theNlambdaTable = new G4PhysicsTable(numOfMaterials) ;

  for (G4int J=0; J<numOfMaterials; J++)
  {
    G4PhysicsLogVector* aVector ;
    aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                                    HighestKineticEnergy,TotBin) ;
    BuildNlambdaVector(aParticleType,J,aVector) ;
    theNlambdaTable->insert(aVector) ;
  }
}

void G4IeBremsstrahlung::BuildNlambdaVector(
                               const G4ParticleDefinition& aParticleType,
                                       G4int materialIndex,
                                       G4PhysicsLogVector* nlambdaVector)
{
  G4double LowEdgeEnergy,T,Tlast,dEdx,Value,Vlast,u,du,coeff ;
  G4double thresholdEnergy ;
  const G4int nbin = 100 ;
  G4bool isOut ;
  const G4double small = 1.e-100;
  const G4double plowloss = 0.5 ;

  const G4MaterialTable* theMaterialTable=
                          G4Material::GetMaterialTable();

  thresholdEnergy = GammaCutInKineticEnergy[materialIndex] ;

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
    if(LowEdgeEnergy >= thresholdEnergy)
    {
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

        Value += coeff*T/(G4EnergyLossTables::GetPreciseDEDX(&aParticleType,
                                   T,(*theMaterialTable)[materialIndex])*
                    (*theMeanFreePathTable)[materialIndex]->GetValue(T,isOut));
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

void G4IeBremsstrahlung::BuildCoeffATable(
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

void G4IeBremsstrahlung::BuildCoeffBTable(
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

void G4IeBremsstrahlung::BuildCoeffCTable(
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

void G4IeBremsstrahlung::BuildInverseNlambdaTable(
                             const G4ParticleDefinition& aParticleType)
{
    G4double T,Smallest,Biggest ;
    const G4double small = 1.e-10;
    G4bool isOut ;

    const G4MaterialTable* theMaterialTable=
                                  G4Material::GetMaterialTable();
    G4int numOfMaterials = theMaterialTable->length();

    if(theInverseNlambdaTable)
    { theInverseNlambdaTable->clearAndDestroy();
      delete theInverseNlambdaTable; }
    theInverseNlambdaTable = new G4PhysicsTable(numOfMaterials);

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

      InvertNlambdaVector(aParticleType,J, aVector);

      theInverseNlambdaTable->insert(aVector);
    }
}

void G4IeBremsstrahlung::InvertNlambdaVector(
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


void G4IeBremsstrahlung::BuildLossTable(
                                 const G4ParticleDefinition& aParticleType)
//  Build table for energy loss due to soft brems
//  tables are built for MATERIALS
//                       *********
{
  G4double KineticEnergy,TotalEnergy,bremloss,Z,x,
           losslim,loss,rate,natom,Cut;

  const G4double MinKineticEnergy = 1.*keV;
  const G4double MinCut = 1.*keV;
  const G4double Thigh = 100.*GeV;
  const G4double Cuthigh = 50.*GeV;
  const G4double Factorhigh = 36./(1450.*GeV);
  const G4double coef1 = -0.5, coef2 = 2./9.;

  ParticleMass = aParticleType.GetPDGMass() ;
  //G4double* GammaCutInKineticEnergy = G4Gamma::Gamma()->GetEnergyCuts();
  GammaCutInKineticEnergy = G4Gamma::Gamma()->GetEnergyCuts();

  //  create table

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length() ;

  if (theLossTable) { theLossTable->clearAndDestroy();
                         delete theLossTable;
                    }

  theLossTable = new G4PhysicsTable(numOfMaterials);

//  loop for materials

  for (G4int J=0; J<numOfMaterials; J++)
    {
     // create physics vector and fill it

     G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                               LowestKineticEnergy,HighestKineticEnergy,TotBin);

     // get elements in the material
     const G4Material* material = (*theMaterialTable)[J];

     const G4ElementVector* theElementVector = material->GetElementVector();
     const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();
     const G4int NumberOfElements = material->GetNumberOfElements();

       //  loop for the kinetic energy values
       for (G4int i=0; i<TotBin; i++)
         {
          KineticEnergy = aVector->GetLowEdgeEnergy(i) ;
          TotalEnergy = KineticEnergy+ParticleMass ;

          Cut = GammaCutInKineticEnergy[J] ;
          if (Cut < MinCut) Cut = MinCut ;
          if (Cut > KineticEnergy) Cut = KineticEnergy ;

          bremloss = 0.;

          if (KineticEnergy>MinKineticEnergy)
            {
             if (Cut > KineticEnergy) Cut = KineticEnergy ;

             //  loop for elements in the material
             for (G4int iel=0; iel<NumberOfElements; iel++)
               {
                Z=(*theElementVector)(iel)->GetZ();
                natom = theAtomicNumDensityVector[iel] ;
                if (KineticEnergy <= Thigh)
                  {
                   //loss for MinKineticEnergy<KineticEnergy<=100 GeV
                   x=log(TotalEnergy/ParticleMass);
                   loss = ComputeBremLoss(Z,natom,KineticEnergy,Cut,x) ;
                   if (&aParticleType==G4Positron::Positron())
                      loss *= ComputePositronCorrFactorLoss(Z,KineticEnergy,Cut) ;
                  }
                else
                  {
                   // extrapolation for KineticEnergy>100 GeV
                   x=log(Thigh/ParticleMass) ;
                   if (Cut<Thigh)
                     {
                      losslim = ComputeBremLoss(Z,natom,Thigh,Cut,x) ;
                      if (&aParticleType==G4Positron::Positron())
                         loss *= ComputePositronCorrFactorLoss(Z,Thigh,Cut) ;
                      rate = Cut/TotalEnergy ;
                      loss = losslim*(1.+coef1*rate+coef2*rate*rate) ;
                      rate = Cut/Thigh ;
                      loss /= (1.+coef1*rate+coef2*rate*rate) ;
                     }
                   else
                     {
                      losslim = ComputeBremLoss(Z,natom,Thigh,Cuthigh,x) ;
                      if (&aParticleType==G4Positron::Positron())
                         loss *= ComputePositronCorrFactorLoss(Z,Thigh,Cuthigh) ;
                      rate = Cut/TotalEnergy ;
                      loss = losslim*(1.+coef1*rate+coef2*rate*rate) ;
                      loss *= Factorhigh*Cut ;
                     }

                  }
                bremloss += natom*loss;
               }

            }

           // now compute the correction due to the LPM effect
           const G4double MigdalConstant = classic_electr_radius*
                                           electron_Compton_length*
                                           electron_Compton_length/pi ;

           const G4double LPMconstant = fine_structure_const*electron_mass_c2*
                                electron_mass_c2/(8.*pi*hbarc) ;
           const G4double kmin = 1.*eV ;
           const G4double klim = 1.*keV ;

           G4double LPMEnergy = LPMconstant*(material->GetRadlen()) ;
           G4double TotalEnergysquare = TotalEnergy*TotalEnergy ;
           G4double LPMGammaEnergyLimit = TotalEnergysquare/LPMEnergy ;

           if(LPMGammaEnergyLimit > klim)
           {
             G4double kmax = min(Cut,LPMGammaEnergyLimit) ;

             G4double floss = 0. ;
             G4int nmax = 1000 ;
             G4int nn ;
             G4double vmin=log(kmin);
             G4double vmax=log(Cut) ;
             nn = int(nmax*(vmax-vmin)/(log(HighestKineticEnergy)-vmin)) ;
             G4double u,uu,s2lpm,sp,fac,c,v,dv,w ;
             dv = (vmax-vmin)/nn ;
             v = vmin-dv ;
             for(G4int n=0; n<=nn; n++)
             {
               v += dv ;
               u = exp(v) ;
               uu = u*u ;
               if(u<=kmax)
               {
                 s2lpm=LPMEnergy*u/TotalEnergysquare ;
                 sp=uu/(uu+MigdalConstant*TotalEnergysquare*
                           (material->GetElectronDensity())) ;
                 w=s2lpm*(1.+1./sp) ;
                 fac=0.5*(sqrt(w*w+4.*s2lpm)-w)/sp;
                 if(fac>1.)
                 fac=1. ;
               }
               else
               {
                 fac=1. ;
               }

               fac *= uu*u ;

               if((n==0)||(n==nn))
                 c=0.5;
               else
                 c=1.;

               fac *= c ;
               floss += fac ;
             }

             floss *=dv*3./(Cut*Cut*Cut-kmin*kmin*kmin) ;
             if(floss > 1.) floss = 1. ;

             // correct the loss
             bremloss *= floss ;
          }
 
          if(bremloss < 0.) bremloss = 0. ;
          aVector->PutValue(i,bremloss);
        }

       theLossTable->insert(aVector);
    }
}

G4double G4IeBremsstrahlung::ComputeXYPolynomial(G4double x,  G4double y,
                                             G4int xSize, G4int ySize,
                                             const G4double coeff[])
{
  // Computes the polynomial (1 y y^2 ...) * matrix * (1 x x^2 ...) .
  // xSize and ySize are the dimensions of the matrix,
  // coeff containts the elements, stored row-wise.    
  G4double* a= new G4double[xSize];
  G4int i, j;

  for (i=0; i<xSize; i++) {
    a[i]= 0.0;
  }  
  G4int index= 0;
  G4double yy= 1.0;
  for (j=0; j<ySize; j++) {
    for (i=0; i<xSize; i++) {
      a[i]+= coeff[index++]*yy;     
    }
    yy*= y;
  }
  
  G4double r= a[0];
  G4double xx= x;
  for (i=1; i<xSize; i++) {
    r+= a[i]*xx;
    xx*= x;
  }  
  
  delete[] a;
  return r;
}                                             

G4double G4IeBremsstrahlung::ComputeBremLoss(G4double Z,G4double natom,
                         G4double T,G4double Cut,G4double x)
{
  const G4double beta=0.99,ksi=2.51,ve=0.00004 ;
  const G4double corrfac = classic_electr_radius*electron_Compton_length*electron_Compton_length/pi  ;

  static const G4double
  CMbarn[]= {
    -0.960613e-1, 0.631029e-1,-0.142819e-1, 0.150437e-2,-0.733286e-4, 0.131404e-5,
     0.859343e-1,-0.529023e-1, 0.131899e-1,-0.159201e-2, 0.926958e-4,-0.208439e-5,
    -0.684096e+1, 0.370364e+1,-0.786752e0,  0.822670e-1,-0.424710e-2, 0.867980e-4,
    -0.200856e+1, 0.129573e+1,-0.306533e0,  0.343682e-1,-0.185931e-2, 0.392432e-4,
     0.127538e+1,-0.515705e0,  0.820644e-1,-0.641997e-2, 0.245913e-3,-0.365789e-5,
     0.115792e0, -0.463143e-1, 0.725442e-2,-0.556266e-3, 0.208049e-4,-0.300895e-6};

  static const G4double
  CPbarn[]= {
    -0.960613e-1, 0.631029e-1,-0.142819e-1, 0.150437e-2,-0.733286e-4, 0.131404e-5,
     0.859343e-1,-0.529023e-1, 0.131899e-1,-0.159201e-2, 0.926958e-4,-0.208439e-5,
    -0.271082e-1, 0.173949e-1,-0.452531e-2, 0.569405e-3,-0.344856e-4, 0.803964e-6,
     0.419855e-2,-0.277188e-2, 0.737658e-3,-0.939463e-4, 0.569748e-5,-0.131737e-6,
    -0.318752e-3, 0.215144e-3,-0.579787e-4, 0.737972e-5,-0.441485e-6, 0.994726e-8,
     0.938233e-5,-0.651642e-5, 0.177303e-5,-0.224680e-6, 0.132080e-7,-0.288593e-9};

  static const G4double
  CCMbarn[]= {
    -0.245667e-3, 0.833406e-4,-0.129217e-4, 0.915099e-6,-0.247179e-7,
     0.147696e-3,-0.498793e-4, 0.402375e-5, 0.989281e-7,-0.133378e-7,
    -0.737702e-2, 0.333057e-2,-0.553141e-3, 0.402464e-4,-0.107977e-5,
    -0.641533e-2, 0.290113e-2,-0.477641e-3, 0.342008e-4,-0.900582e-6,
     0.574303e-5, 0.908521e-4,-0.256900e-4, 0.239921e-5,-0.741271e-7};

  static const G4double
  CCPbarn[]= {
    -0.245667e-3, 0.833406e-4,-0.129217e-4, 0.915099e-6,-0.247179e-7,
     0.147696e-3,-0.498793e-4, 0.402375e-5, 0.989281e-7,-0.133378e-7,
    -0.341260e-4, 0.971711e-5,-0.172031e-6,-0.119455e-6, 0.704166e-8,
     0.341740e-5,-0.775867e-6,-0.653231e-7, 0.225605e-7,-0.114860e-8,
    -0.119391e-6, 0.194885e-7, 0.588959e-8,-0.127589e-8, 0.608247e-10};

  G4double CM[36],CP[36],CCM[25],CCP[25]; //Set the unit: barn

  for (G4int i=0; i<36; i++)    { CM[i] = CMbarn[i]*barn;
                                  CP[i] = CPbarn[i]*barn;
                                }
  for (G4int ii=0; ii<25; ii++) { CCM[ii] = CCMbarn[ii]*barn;
                                  CCP[ii] = CCPbarn[ii]*barn;
                                }
   //  -----------------------------------------------------------

  G4double TotalEnergy = T + electron_mass_c2;
  G4double y=log(Cut/(ve*TotalEnergy));

  G4double loss;

  if (y <= 0.) loss = ComputeXYPolynomial(x, y, 6, 6, CM)
                     + Z * ComputeXYPolynomial(x, y, 5, 5, CCM);
  else         loss = ComputeXYPolynomial(x, y, 6, 6, CP)
                     + Z * ComputeXYPolynomial(x, y, 5, 5, CCP);

  G4double rate = TotalEnergy/Cut ;
  G4double corr = 1./(1.+corrfac*natom*rate*rate) ;

  G4double factor = pow(Cut*corr/T,beta);
  factor *= Z*(Z+ksi)*TotalEnergy*TotalEnergy/(TotalEnergy+electron_mass_c2) ;

  loss   *= factor ;

  return loss ;
}


  G4double G4IeBremsstrahlung::ComputePositronCorrFactorLoss(
                             G4double Z,G4double KineticEnergy,G4double GammaCut)

  // calculates the correction factor for the energy loss due to soft bremsstrahlung for positrons
  //***********************************************************************************************
  //  the same correction is in the (discrete) bremsstrahlung 
  //**********************************************************************************************  
  {
    static const G4double K = 132.9416*eV ;
    static const G4double a1=4.15e-1,a3=2.10e-3,a5=54.0e-5 ;
    G4double x,x2,x3,eta,e0,factor ;

    x = log(KineticEnergy/(K*Z*Z)) ;
    x2 = x*x ;
    x3 = x2*x ;
    eta = 0.5+atan(a1*x+a3*x3+a5*x3*x2)/pi ;
    e0 = GammaCut/KineticEnergy ;

    if (e0==1.0) {
      factor = 0;
    } else {
      factor = log(1.-e0)/eta ;
      factor = exp(factor) ;
    }  
    factor = eta*(1.-factor)/e0 ;

    return factor ;
  }
      
void G4IeBremsstrahlung::BuildLambdaTable(const G4ParticleDefinition& ParticleType)

// Build  mean free path tables for the gamma emission by e- or e+.
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

   PartialSumSigma.resize( G4Material::GetNumberOfMaterials() );
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

void G4IeBremsstrahlung::ComputePartialSumSigma(const G4ParticleDefinition* ParticleType,
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

G4double G4IeBremsstrahlung::ComputeMicroscopicCrossSection(
                                     const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy, G4double AtomicNumber,
                                           G4double GammaEnergyCut)
 
// Calculates the microscopic cross section in GEANT4 internal units.
// A parametrized formula from L. Urban is used to estimate the total cross section.
// This parametrization is derived from :
//        tabulated cross-section values of Seltzer and Berger below 10 GeV,
//        screened Bethe Heilter differential cross section above 10 GeV,
//        Migdal corrections in both case. 
//  Seltzer & Berger: Nim B 12:95 (1985)
//  Nelson, Hirayama & Rogers: Technical report 265 SLAC (1985)
//  Migdal: Phys Rev 103:1811 (1956); Messel & Crawford: Pergamon Press (1970)
//
// Above 100 GeV the Cross section is scaled in log(KineticEnergy).
 
{
 G4double CrossSection = 0.0 ;
 if ( AtomicNumber < 1. ) return CrossSection;
 if ( KineticEnergy < 10*keV ) return CrossSection;
 if ( KineticEnergy <= GammaEnergyCut ) return CrossSection;

 G4double LocalKineticEnergy = KineticEnergy,  LocalGammaEnergyCut = GammaEnergyCut;

 static const G4double KinLimitScale = 100.*GeV,   CutLimitScale = 50.*GeV;

 if ( KineticEnergy > KinLimitScale )
     { LocalKineticEnergy = KinLimitScale;
       if (GammaEnergyCut > KinLimitScale) LocalGammaEnergyCut = CutLimitScale;
     }

 static const G4double
      aay0x0= 0.430748E-02*barn, aay0x1= 0.576058E-02*barn, aay0x2=-0.122564E-02*barn,
      aay0x3= 0.114843E-03*barn, aay0x4=-0.489452E-05*barn, aay0x5= 0.795991E-07*barn;

 static const G4double
      aay1x0= 0.326746E-02*barn, aay1x1=-0.132872E-02*barn, aay1x2= 0.217197E-03*barn,
      aay1x3=-0.179769E-04*barn, aay1x4= 0.766114E-06*barn, aay1x5=-0.125603E-07*barn;

 static const G4double
      amy2x0= 0.326452E-02*barn, amy2x1=-0.175331E-02*barn, amy2x2= 0.415488E-03*barn,
      amy2x3=-0.507652E-04*barn, amy2x4= 0.297569E-05*barn, amy2x5=-0.651741E-07*barn;

 static const G4double
      amy3x0= 0.847189E-03*barn, amy3x1=-0.433923E-03*barn, amy3x2= 0.116672E-03*barn,
      amy3x3=-0.166799E-04*barn, amy3x4= 0.110237E-05*barn, amy3x5=-0.263383E-07*barn;

 static const G4double
      amy4x0= 0.846052E-04*barn, amy4x1=-0.415764E-04*barn, amy4x2= 0.129610E-04*barn,
      amy4x3=-0.212844E-05*barn, amy4x4= 0.152871E-06*barn, amy4x5=-0.384393E-08*barn;

 static const G4double
      amy5x0= 0.300838E-05*barn, amy5x1=-0.136833E-05*barn, amy5x2= 0.507296E-06*barn,
      amy5x3=-0.943623E-07*barn, amy5x4= 0.720305E-08*barn, amy5x5=-0.187210E-09*barn;

 static const G4double
      apy2x0= 0.448230E-01*barn, apy2x1=-0.210048E-01*barn, apy2x2= 0.379434E-02*barn,
      apy2x3=-0.328431E-03*barn, apy2x4= 0.136710E-04*barn, apy2x5=-0.220593E-06*barn;

 static const G4double
      apy3x0=-0.539248E-02*barn, apy3x1= 0.330244E-02*barn, apy3x2=-0.733726E-03*barn,
      apy3x3= 0.732312E-04*barn, apy3x4=-0.336810E-05*barn, apy3x5= 0.583913E-07*barn;

 static const G4double
      apy4x0=-0.106983E-02*barn, apy4x1= 0.378021E-03*barn, apy4x2=-0.384854E-04*barn,
      apy4x3= 0.978156E-06*barn, apy4x4= 0.410622E-07*barn, apy4x5=-0.174250E-08*barn;

 static const G4double
      apy5x0=-0.117501E-04*barn, apy5x1=-0.983887E-05*barn, apy5x2= 0.239644E-05*barn,
      apy5x3=-0.190104E-06*barn, apy5x4= 0.619226E-08*barn, apy5x5=-0.680932E-10*barn;

 static const G4double
      bby0x0= 0.168074E-03*barn, bby0x1=-0.934609E-04*barn, bby0x2= 0.141293E-04*barn,
      bby0x3=-0.854216E-06*barn, bby0x4= 0.183287E-07*barn;

 static const G4double
      bby1x0= 0.932144E-04*barn, bby1x1=-0.234926E-04*barn, bby1x2= 0.136656E-05*barn,
      bby1x3= 0.351109E-07*barn, bby1x4=-0.330189E-08*barn;

 static const G4double
      bmy2x0= 0.174523E-04*barn, bmy2x1= 0.253854E-05*barn, bmy2x2=-0.171643E-05*barn,
      bmy2x3= 0.183074E-06*barn, bmy2x4=-0.566331E-08*barn;

 static const G4double
      bmy3x0= 0.111970E-05*barn, bmy3x1= 0.112776E-05*barn, bmy3x2=-0.386924E-06*barn,
      bmy3x3= 0.367597E-07*barn, bmy3x4=-0.108504E-08*barn;

 static const G4double
      bmy4x0= 0.171604E-07*barn, bmy4x1= 0.738801E-07*barn, bmy4x2=-0.218761E-07*barn,
      bmy4x3= 0.199032E-08*barn, bmy4x4=-0.576173E-10*barn;

 static const G4double
      bpy2x0=-0.105531E-03*barn, bpy2x1= 0.362995E-04*barn, bpy2x2=-0.433334E-05*barn,
      bpy2x3= 0.207664E-06*barn, bpy2x4=-0.330250E-08*barn;

 static const G4double
      bpy3x0=-0.168293E-05*barn, bpy3x1=-0.773204E-06*barn, bpy3x2= 0.227974E-06*barn,
      bpy3x3=-0.159385E-07*barn, bpy3x4= 0.321958E-09*barn;

 static const G4double
      bpy4x0= 0.167046E-05*barn, bpy4x1=-0.440761E-06*barn, bpy4x2= 0.396377E-07*barn,
      bpy4x3=-0.151053E-08*barn, bpy4x4= 0.215624E-10*barn;

 static const G4double ksi=1.8, alfa=0.98, vs= 1.E-4;

 G4double TotalEnergy = LocalKineticEnergy + electron_mass_c2;
 G4double X = log(TotalEnergy/electron_mass_c2),  X2=X*X, X3=X2*X, X4=X3*X, X5=X4*X;
 G4double Y = log(vs*TotalEnergy/LocalGammaEnergyCut), Y2=Y*Y, Y3=Y2*Y, Y4=Y3*Y, Y5=Y4*Y;

 G4double ay0, ay1, ay2, ay3, ay4, ay5, by0, by1, by2, by3, by4;
 if (Y < 0.) {
    ay0 = aay0x0 + aay0x1*X + aay0x2*X2 + aay0x3*X3 + aay0x4*X4 + aay0x5*X5;
    ay1 = aay1x0 + aay1x1*X + aay1x2*X2 + aay1x3*X3 + aay1x4*X4 + aay1x5*X5;
    ay2 = amy2x0 + amy2x1*X + amy2x2*X2 + amy2x3*X3 + amy2x4*X4 + amy2x5*X5;
    ay3 = amy3x0 + amy3x1*X + amy3x2*X2 + amy3x3*X3 + amy3x4*X4 + amy3x5*X5;
    ay4 = amy4x0 + amy4x1*X + amy4x2*X2 + amy4x3*X3 + amy4x4*X4 + amy4x5*X5;
    ay5 = amy5x0 + amy5x1*X + amy5x2*X2 + amy5x3*X3 + amy5x4*X4 + amy5x5*X5;

    by0 = bby0x0 + bby0x1*X + bby0x2*X2 + bby0x3*X3 + bby0x4*X4;
    by1 = bby1x0 + bby1x1*X + bby1x2*X2 + bby1x3*X3 + bby1x4*X4;
    by2 = bmy2x0 + bmy2x1*X + bmy2x2*X2 + bmy2x3*X3 + bmy2x4*X4;
    by3 = bmy3x0 + bmy3x1*X + bmy3x2*X2 + bmy3x3*X3 + bmy3x4*X4;
    by4 = bmy4x0 + bmy4x1*X + bmy4x2*X2 + bmy4x3*X3 + bmy4x4*X4;
  }
 else {
    ay0 = aay0x0 + aay0x1*X + aay0x2*X2 + aay0x3*X3 + aay0x4*X4 + aay0x5*X5;
    ay1 = aay1x0 + aay1x1*X + aay1x2*X2 + aay1x3*X3 + aay1x4*X4 + aay1x5*X5;
    ay2 = apy2x0 + apy2x1*X + apy2x2*X2 + apy2x3*X3 + apy2x4*X4 + apy2x5*X5;
    ay3 = apy3x0 + apy3x1*X + apy3x2*X2 + apy3x3*X3 + apy3x4*X4 + apy3x5*X5;
    ay4 = apy4x0 + apy4x1*X + apy4x2*X2 + apy4x3*X3 + apy4x4*X4 + apy4x5*X5;
    ay5 = apy5x0 + apy5x1*X + apy5x2*X2 + apy5x3*X3 + apy5x4*X4 + apy5x5*X5;

    by0 = bby0x0 + bby0x1*X + bby0x2*X2 + bby0x3*X3 + bby0x4*X4;
    by1 = bby1x0 + bby1x1*X + bby1x2*X2 + bby1x3*X3 + bby1x4*X4;
    by2 = bpy2x0 + bpy2x1*X + bpy2x2*X2 + bpy2x3*X3 + bpy2x4*X4;
    by3 = bpy3x0 + bpy3x1*X + bpy3x2*X2 + bpy3x3*X3 + bpy3x4*X4;
    by4 = bpy4x0 + bpy4x1*X + bpy4x2*X2 + bpy4x3*X3 + bpy4x4*X4;
  }

 G4double F0 = ay0 + ay1*Y + ay2*Y2 + ay3*Y3 + ay4*Y4 + ay5*Y5,
          F1 = by0 + by1*Y + by2*Y2 + by3*Y3 + by4*Y4;

 CrossSection = AtomicNumber*(AtomicNumber+ksi)*TotalEnergy*TotalEnergy
               * pow(log(LocalKineticEnergy/LocalGammaEnergyCut),alfa)
               * (F0 + F1*AtomicNumber)
               / (LocalKineticEnergy*(LocalKineticEnergy+2*electron_mass_c2));

 if (ParticleType == G4Positron::Positron())
     CrossSection *= ComputePositronCorrFactorSigma(AtomicNumber, LocalKineticEnergy,
                                                             LocalGammaEnergyCut);

 // now comes the scaling above 100GeV
 if (KineticEnergy > KinLimitScale)
     { G4double X1 = GammaEnergyCut/KineticEnergy, 
                X2 = LocalGammaEnergyCut/LocalKineticEnergy;
       CrossSection *= (-log(X1) -2./3. + X1 - X1*X1/3.)/(-log(X2) -2./3. + X2 - X2*X2/3.);
     }

 if (CrossSection < 0.) CrossSection = 0.;
 return CrossSection;
}
 
G4double G4IeBremsstrahlung::ComputePositronCorrFactorSigma( G4double AtomicNumber,
                                           G4double KineticEnergy, G4double GammaEnergyCut)
 
// Calculates the correction factor for the total cross section of the positron bremsstrahl.
// Eta is the ratio of positron to electron energy loss by bremstrahlung. 
// A parametrized formula from L. Urban is used to estimate eta. It is a fit to the results
// of L. Kim & al: Phys Rev. A33,3002 (1986)
 
{
 static const G4double K = 132.9416*eV;
 static const G4double a1 = 4.15e-1, a3 = 2.10e-3, a5 = 54.0e-5;

 G4double x = log(KineticEnergy/(K*AtomicNumber*AtomicNumber));
 G4double eta = 0.5 + atan(a1*x + a3*x*x*x + a5*x*x*x*x*x)/pi ;
 G4double alfa = (1. - eta)/eta;
 return eta*pow((1. - GammaEnergyCut/KineticEnergy) , alfa);
}


G4VParticleChange* G4IeBremsstrahlung::PostStepDoIt(const G4Track& trackData,
                                                  const G4Step& stepData)
//
// The emitted gamma energy is sampled using a parametrized formula from L. Urban.
// This parametrization is derived from :
//    cross-section values of Seltzer and Berger for electron energies 1 keV - 10 GeV,
//    screened Bethe Heilter differential cross section above 10 GeV,
//    Migdal corrections in both case.
//  Seltzer & Berger: Nim B 12:95 (1985)
//  Nelson, Hirayama & Rogers: Technical report 265 SLAC (1985)
//  Migdal: Phys Rev 103:1811 (1956); Messel & Crawford: Pergamon Press (1970)
//
// A modified version of the random number techniques of Butcher & Messel is used
//    (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
{

  static const G4double
     ah10 = 4.67733E+00, ah11 =-6.19012E-01, ah12 = 2.02225E-02,
     ah20 =-7.34101E+00, ah21 = 1.00462E+00, ah22 =-3.20985E-02,
     ah30 = 2.93119E+00, ah31 =-4.03761E-01, ah32 = 1.25153E-02;

  static const G4double
     bh10 = 4.23071E+00, bh11 =-6.10995E-01, bh12 = 1.95531E-02,
     bh20 =-7.12527E+00, bh21 = 9.69160E-01, bh22 =-2.74255E-02,
     bh30 = 2.69925E+00, bh31 =-3.63283E-01, bh32 = 9.55316E-03;

  static const G4double
     al00 =-2.05398E+00, al01 = 2.38815E-02, al02 = 5.25483E-04,
     al10 =-7.69748E-02, al11 =-6.91499E-02, al12 = 2.22453E-03,
     al20 = 4.06463E-02, al21 =-1.01281E-02, al22 = 3.40919E-04;

  static const G4double
     bl00 = 1.04133E+00, bl01 =-9.43291E-03, bl02 =-4.54758E-04,
     bl10 = 1.19253E-01, bl11 = 4.07467E-02, bl12 =-1.30718E-03,
     bl20 =-1.59391E-02, bl21 = 7.27752E-03, bl22 =-1.94405E-04;

  static const G4double MigdalConstant = classic_electr_radius
                                        *electron_Compton_length
                                        *electron_Compton_length/pi;
  const G4double LPMconstant = fine_structure_const*electron_mass_c2*
                                electron_mass_c2/(8.*pi*hbarc) ;
   aParticleChange.Initialize(trackData);
   G4Material* aMaterial=trackData.GetMaterial() ;

   G4double LPMEnergy = LPMconstant*(aMaterial->GetRadlen()) ;

   const G4DynamicParticle* aDynamicParticle=trackData.GetDynamicParticle();
   G4double charge = aDynamicParticle->GetDefinition()->GetPDGCharge();

   G4double           KineticEnergy     = aDynamicParticle->GetKineticEnergy();
   G4ParticleMomentum ParticleDirection = aDynamicParticle->GetMomentumDirection();

   // Gamma production cut in this material
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

   // Extract Z factors for this Element
   G4double lnZ = 3.*(anElement->GetIonisation()->GetlogZ3());
   G4double FZ = lnZ* (4.- 0.55*lnZ);
   G4double ZZ = anElement->GetIonisation()->GetZZ3();

   // limits of the energy sampling
   G4double TotalEnergy = KineticEnergy + electron_mass_c2;
   G4double TotalEnergysquare = TotalEnergy*TotalEnergy ;
   G4double LPMGammaEnergyLimit = TotalEnergysquare/LPMEnergy ;

   G4double xmin = GammaEnergyCut/KineticEnergy, epsilmin = GammaEnergyCut/TotalEnergy;
   G4double epsilmax = KineticEnergy/TotalEnergy;

   // Migdal factor
   G4double MigdalFactor = (aMaterial->GetElectronDensity())*MigdalConstant
                          /(epsilmax*epsilmax);

   //
   G4double x, epsil, greject, migdal, grejmax;
   G4double U = log(KineticEnergy/electron_mass_c2), U2 = U*U;

   //
   //  sample the energy rate of the emitted gamma for electron kinetic energy > 1 MeV
   //

   if (KineticEnergy > 1.*MeV)
     {
       // parameters
       G4double ah1 = ah10 + ZZ* (ah11 + ZZ* ah12),
                ah2 = ah20 + ZZ* (ah21 + ZZ* ah22),
                ah3 = ah30 + ZZ* (ah31 + ZZ* ah32);

       G4double bh1 = bh10 + ZZ* (bh11 + ZZ* bh12),
                bh2 = bh20 + ZZ* (bh21 + ZZ* bh22),
                bh3 = bh30 + ZZ* (bh31 + ZZ* bh32);

       G4double ah = 1.   + (ah1*U2 + ah2*U + ah3) / (U2*U);
       G4double bh = 0.75 + (bh1*U2 + bh2*U + bh3) / (U2*U);

       // limit of the screening variable
       G4double screenfac =
       136.*electron_mass_c2/((anElement->GetIonisation()->GetZ3())*TotalEnergy);
       G4double screenmin = screenfac*epsilmin/(1.-epsilmin);

       // Compute the maximum of the rejection function
       G4double F1 = max(ScreenFunction1(screenmin) - FZ ,0.);
       G4double F2 = max(ScreenFunction2(screenmin) - FZ ,0.);
       grejmax = (F1 - epsilmin* (F1*ah - bh*epsilmin*F2))/(42.392 - FZ);

       // sample the energy rate of the emitted Gamma
       G4double screenvar;

       do {

             x = pow(xmin, G4UniformRand()); 
             epsil = x*KineticEnergy/TotalEnergy;
             screenvar = screenfac*epsil/(1-epsil);
             F1 = max(ScreenFunction1(screenvar) - FZ ,0.);
             F2 = max(ScreenFunction2(screenvar) - FZ ,0.);
             migdal = (1. + MigdalFactor)/(1. + MigdalFactor/(x*x));
             greject = migdal*(F1 - epsil* (ah*F1 - bh*epsil*F2))/(42.392 - FZ);
        }  while( greject < G4UniformRand()*grejmax );

    }
   else
     { 
       // sample the energy rate of the emitted gamma for electron kinetic energy < 1 MeV
       //
       // parameters
       G4double al0 = al00 + ZZ* (al01 + ZZ* al02),
                al1 = al10 + ZZ* (al11 + ZZ* al12),
                al2 = al20 + ZZ* (al21 + ZZ* al22);

       G4double bl0 = bl00 + ZZ* (bl01 + ZZ* bl02),
                bl1 = bl10 + ZZ* (bl11 + ZZ* bl12),
                bl2 = bl20 + ZZ* (bl21 + ZZ* bl22);

       G4double al = al0 + al1*U + al2*U2;
       G4double bl = bl0 + bl1*U + bl2*U2;

       // Compute the maximum of the rejection function
       grejmax = max(1. + xmin* (al + bl*xmin), 1.+al+bl);
       G4double xm = -al/(2.*bl);
       if ((xmin < xm)&&(xm < 1.)) grejmax = max(grejmax, 1.+ xm* (al + bl*xm));

       // sample the energy rate of the emitted Gamma

       do {  x = pow(xmin, G4UniformRand());
             migdal = (1. + MigdalFactor)/(1. + MigdalFactor/(x*x));
             greject = migdal*(1. + x* (al + bl*x));
        }  while( greject < G4UniformRand()*grejmax );
   }

   G4double GammaEnergy = x*KineticEnergy;

   // now comes the supression due to the LPM effect
   if(GammaEnergy < LPMGammaEnergyLimit)
   {
     G4double S2LPM = LPMEnergy*GammaEnergy/TotalEnergysquare ;
     G4double Spol  = GammaEnergy*GammaEnergy/(GammaEnergy*GammaEnergy +
                      MigdalConstant*(aMaterial->GetElectronDensity())*
                      TotalEnergysquare) ;
     G4double w = S2LPM*(1.+1./Spol) ;
     G4double Supr = 0.5*(sqrt(w*w+4.*S2LPM)-w)/Spol ;

     //
     if (G4UniformRand() > Supr )
       GammaEnergy = 0. ;
   }

   //protection: DO NOT PRODUCE a gamma with energy 0. !
   if (GammaEnergy <= 0.)
       return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

   //
   //  angles of the emitted gamma. ( Z - axis along the parent particle)
   //
   //  universal distribution suggested by L. Urban (Geant3 manual (1993) Phys211),
   //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

   G4double u;
   const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

   if (9./(9.+d) > G4UniformRand()) u = - log(G4UniformRand()*G4UniformRand())/a1 ;
      else                          u = - log(G4UniformRand()*G4UniformRand())/a2 ;

   G4double Teta = u*electron_mass_c2/TotalEnergy ;
   G4double Phi  = twopi * G4UniformRand() ;
   G4double dirx = sin(Teta)*cos(Phi) , diry = sin(Teta)*sin(Phi) , dirz = cos(Teta) ;

   G4ThreeVector GammaDirection ( dirx, diry, dirz);
   GammaDirection.rotateUz(ParticleDirection);

   // create G4DynamicParticle object for the Gamma
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
      aParticleChange.SetMomentumChange( ParticleDirection );
      aParticleChange.SetEnergyChange( NewKinEnergy );
      aParticleChange.SetLocalEnergyDeposit (0.);
     }
   else
     {
      aParticleChange.SetEnergyChange( 0. );
      aParticleChange.SetLocalEnergyDeposit (0.);
      if (charge<0.) aParticleChange.SetStatusChange(fStopAndKill);
          else       aParticleChange.SetStatusChange(fStopButAlive);
     }   

   return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}
G4Element* G4IeBremsstrahlung::SelectRandomAtom(G4Material* aMaterial) const
{
  // select randomly 1 element within the material

  const G4int Index = aMaterial->GetIndex();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  G4double rval = G4UniformRand()*((*PartialSumSigma(Index))(NumberOfElements-1));
  for ( G4int i=0; i < NumberOfElements; i++ )
    if (rval <= (*PartialSumSigma(Index))(i)) return ((*theElementVector)(i));
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements, NULL pointer returned." << endl;
  return NULL;
}

void G4IeBremsstrahlung::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from a parametrisation(L.Urban). ";
           comments += "Good description from 1 KeV to 100 GeV.\n";
           comments += "        log scale extrapolation above 100 GeV \n";
           comments += "        Gamma energy sampled from a parametrised formula.";

  G4cout << endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,"Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins. \n";
}        


