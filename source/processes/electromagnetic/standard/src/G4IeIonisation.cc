// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IeIonisation.cc,v 1.1 1999-01-07 16:11:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4IeIonisation physics process -----------
//                by Laszlo Urban, 23 June 1998
// ************************************************************
// It is the first implementation of the IONISATION
// PROCESS. ( delta rays + continuous energy loss)
//   using an INTEGRAL APPROACH instead of the differential
//   one used in the standard implementation .
// ************************************************************
// 27/10/98: minor changes + cleanup, L.Urban
// --------------------------------------------------------------
 

#include "G4IeIonisation.hh"
#include "G4UnitsTable.hh"

 
G4IeIonisation::G4IeIonisation(const G4String& processName)
   : G4IeEnergyLoss(processName),
     theMeanFreePathTable(NULL),
     theNlambdaTable(NULL),
     theInverseNlambdaTable(NULL),
     theCoeffATable(NULL),
     theCoeffBTable(NULL),
     theCoeffCTable(NULL),
     LowestKineticEnergy(1.00*keV),
     HighestKineticEnergy(100.*TeV),
     TotBin(100),
     NumberOfBuildPhysicsTableCalls(0),
     theElectron ( G4Electron::Electron() ),
     thePositron ( G4Positron::Positron() )
{  }

G4IeIonisation::~G4IeIonisation() 
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
     if(theCoeffATable) {
        theCoeffATable->clearAndDestroy();
        delete theCoeffATable;
     }
     if(theCoeffBTable) {
        theCoeffBTable->clearAndDestroy();
        delete theCoeffBTable;
     }
     if(theCoeffCTable) {
        theCoeffCTable->clearAndDestroy();
        delete theCoeffCTable;
     }
}
 
void G4IeIonisation::SetPhysicsTableBining(
                                  G4double lowE,G4double highE,G4int nBins)
{
  LowestKineticEnergy=lowE; HighestKineticEnergy=highE; TotBin=nBins;  
}

void G4IeIonisation::BuildPhysicsTable(
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

    if(&aParticleType==G4Electron::Electron())
      PrintInfoDefinition();
  }
}

void G4IeIonisation::TestOfInversion(
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
    G4cout << "G4IeIonisation::TestOfInversion (T->Nlambda->Tprime) " << endl ; 
    G4cout << "particle= " << aParticleType.GetParticleName() <<
            "   material= " << (*theMaterialTable)[J]->GetName() << endl ;  
    G4cout << "max (Tprime-T)/T in % =" << setw(10) << setprecision(3) <<
                                                                    delta ; 
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

void G4IeIonisation::BuildNlambdaTable(
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

void G4IeIonisation::BuildNlambdaVector(
                               const G4ParticleDefinition& aParticleType,
                                       G4int materialIndex,
                                       G4PhysicsLogVector* nlambdaVector)
{
  G4double LowEdgeEnergy,T,Tlast,dEdx,Value,Vlast,u,du,coeff ;
  G4double thresholdEnergy ;
  const G4int nbin = 100 ;
  G4bool isOut ;
  const G4double small = 1.e-100;
  const G4double plowloss = 0.5 ;  //this should be a data member of en.loss!

  const G4MaterialTable* theMaterialTable=
                          G4Material::GetMaterialTable();

  if(&aParticleType == theElectron)
    thresholdEnergy = 2.*DeltaCutInKineticEnergy[materialIndex];
  else
    thresholdEnergy =    DeltaCutInKineticEnergy[materialIndex];

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

void G4IeIonisation::BuildCoeffATable(
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

  for (G4int J=0; J<numOfMaterials; J++)
  {
    G4int binmax=TotBin ;
    G4PhysicsLinearVector* aVector =
                            new G4PhysicsLinearVector(0.,binmax, TotBin);

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

void G4IeIonisation::BuildCoeffBTable(
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

  for (G4int J=0; J<numOfMaterials; J++)
  {
    G4int binmax=TotBin ;
    G4PhysicsLinearVector* aVector =
                           new G4PhysicsLinearVector(0.,binmax, TotBin);

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

void G4IeIonisation::BuildCoeffCTable(
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

  for (G4int J=0; J<numOfMaterials; J++)
  {

    G4int binmax=TotBin ;
    G4PhysicsLinearVector* aVector =
                          new G4PhysicsLinearVector(0.,binmax, TotBin);


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

void G4IeIonisation::BuildInverseNlambdaTable(
                             const G4ParticleDefinition& aParticleType)
{
    G4double T,Smallest,Biggest ;
    const G4double small = 1.e-100;
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
  G4cout << "material=" << (*theMaterialTable)[J]->GetName() <<
            "  smallest=" << Smallest << "  biggest=" << Biggest << endl;
         G4Exception(
        "Cut value is too big , smaller value should be used !");  
      }

      G4PhysicsLogVector* aVector;

      aVector = new G4PhysicsLogVector(Smallest,
                            Biggest,TotBin);

      InvertNlambdaVector(aParticleType,J, aVector);

      theInverseNlambdaTable->insert(aVector);
    }
}

void G4IeIonisation::InvertNlambdaVector(
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

void G4IeIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
    G4double LowEdgeEnergy , ionloss ;
    G4bool isOutRange ;
    const G4MaterialTable* theMaterialTable=
                                     G4Material::GetMaterialTable();
    const G4double twoln10 = 2.*log(10.) ;
    const G4double Factor = twopi_mc2_rcl2 ;

    ParticleMass = aParticleType.GetPDGMass();

    ParticleCutInKineticEnergy = aParticleType.GetEnergyCuts() ;

    G4int numOfMaterials = theMaterialTable->length();

     if (theLossTable) {
        theLossTable->clearAndDestroy();
        delete theLossTable;
     }

    theLossTable = new G4PhysicsTable(numOfMaterials);

    for (G4int J=0; J<numOfMaterials; J++)
    {
      G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);
  
      G4double ElectronDensity,Eexc,Eexcm2,Cden,Mden,Aden,X0den,X1den ;
      const G4Material* material= (*theMaterialTable)[J];
      ElectronDensity = material->GetElectronDensity();
      Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
      Eexcm2 = Eexc/ParticleMass ;
      Eexcm2 *= Eexcm2 ;
      Cden = material->GetIonisation()->GetCdensity();
      Mden = material->GetIonisation()->GetMdensity();
      Aden = material->GetIonisation()->GetAdensity();
      X0den = material->GetIonisation()->GetX0density();
      X1den = material->GetIonisation()->GetX1density();

      ParticleCutInKineticEnergyNow = ParticleCutInKineticEnergy[J] ;

      // some local variables -------------------
      G4double tau,Tmax,gamma,gamma2,bg2,beta2,d,d2,d3,d4,delta,x,y ;

      for (G4int i = 0 ; i < TotBin ; i++)
      {
        LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
        tau = LowEdgeEnergy/ParticleMass ;

          // Seltzer-Berger formula 
          gamma = tau +1. ;
          bg2 = tau*(tau+2.) ;
          gamma2 = gamma*gamma ;
          beta2 = bg2/gamma2 ;

          // electron .................................
          if(&aParticleType==theElectron)
          {
            Tmax = LowEdgeEnergy/2. ;
  
            d = min(ParticleCutInKineticEnergyNow, Tmax)/ParticleMass;

            ionloss = log(2.*(tau+2.)/Eexcm2)-1.-beta2 ;
            ionloss += log((tau-d)*d)+tau/(tau-d) ;
            ionloss += (0.5*d*d+(2.*tau+1.)*log(1.-d/tau))/gamma2 ;
          }
          //positron ...............................
          else
          {
            Tmax = LowEdgeEnergy ;
  
            d = min(ParticleCutInKineticEnergyNow, Tmax)/ParticleMass;

            d2=d*d/2. ;
            d3=d*d*d/3. ;
            d4=d*d*d*d/4. ;

            y=1./(1.+gamma) ;

            ionloss = log(2.*(tau+2.)/Eexcm2)+log(tau*d) ;
            ionloss-= beta2*(tau+2.*d-y*(3.*d2+y*(d-d3+y*(d2-tau*d3+d4))))/tau;

          } 

          // density correction   ................................
          x = log(bg2)/twoln10 ;
          if ( x < X0den )
             delta = 0. ;
          else 
          {
             delta = twoln10*x - Cden ;
             if ( x < X1den )
               delta += Aden*pow((X1den-x),Mden) ;
          } 

          // now you can compute the total ionization loss
          ionloss -= delta ;
          ionloss *= Factor*ElectronDensity/beta2 ;
          if ( ionloss <= 0.)
             ionloss = 0.   ;
   
          aVector->PutValue(i,ionloss) ;
      }  
      theLossTable->insert(aVector);
    }
}

void G4IeIonisation::BuildLambdaTable(const G4ParticleDefinition& aParticleType)
{
      G4double LowEdgeEnergy , Value ,sigma ;
      G4bool isOutRange ;
      const G4MaterialTable* theMaterialTable=
                                         G4Material::GetMaterialTable();

      G4int numOfMaterials = theMaterialTable->length();

      if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
      }

      theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);

      DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
      ParticleCutInKineticEnergy = aParticleType.GetEnergyCuts() ;

      for (G4int J=0 ; J < numOfMaterials; J++)
      { 
        G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
               LowestKineticEnergy, HighestKineticEnergy, TotBin);

        const G4Material* material= (*theMaterialTable)[J];
        
        const G4ElementVector* theElementVector=
                         material->GetElementVector() ;
        const G4double* theAtomicNumDensityVector =
                         material->GetAtomicNumDensityVector();
        const G4int NumberOfElements=
                         material->GetNumberOfElements() ;
 
        DeltaKineticEnergyCutNow = DeltaCutInKineticEnergy[J] ;

        for ( G4int i = 0 ; i < TotBin ; i++ )
        {
           LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
           sigma = 0. ;
           for (G4int iel=0; iel<NumberOfElements; iel++ )
           {
                sigma +=  theAtomicNumDensityVector[iel]*
                        ComputeMicroscopicCrossSection(aParticleType,
                          LowEdgeEnergy,
                          (*theElementVector)(iel)->GetZ() ) ;
           }
           Value = sigma<=0 ? DBL_MAX : 1./sigma ;     
           aVector->PutValue(i, Value) ;
        }
        theMeanFreePathTable->insert(aVector);
      }
}


G4double G4IeIonisation::ComputeMicroscopicCrossSection(
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber)
{
  // calculates the microscopic cross section in GEANT4 internal units
  //    ( it is called for elements , AtomicNumber = Z )
    G4double TotalEnergy,MaxKineticEnergyTransfer,
             betasquare,gamma,gamma2,x,x2,y,y2,y12,b1,b2,b3,b4,
             TotalCrossSection;

    TotalEnergy=KineticEnergy + ParticleMass;
    betasquare = KineticEnergy*(TotalEnergy+ParticleMass)
                /(TotalEnergy*TotalEnergy);
    gamma = TotalEnergy/ParticleMass ;
    x=DeltaKineticEnergyCutNow/KineticEnergy ;

    if(&aParticleType==theElectron)
      MaxKineticEnergyTransfer = 0.5*KineticEnergy ;
    else
      MaxKineticEnergyTransfer =     KineticEnergy ;

    if( MaxKineticEnergyTransfer > DeltaKineticEnergyCutNow )
    {
      //Moller (e-e-) scattering
      if(&aParticleType==theElectron)
      {
        gamma2 = gamma*gamma ;
        TotalCrossSection = (gamma-1.)*(gamma-1.)*(0.5-x)/gamma2 + 1./x -
                          1./(1.-x)-(2.*gamma-1.)*log((1.-x)/x)/gamma2 ; 
        TotalCrossSection /= betasquare ; 
      }
      //Bhabha (e+e-) scattering
      else
      {
        x2=x*x ;
        y=1./(1.+gamma) ;
        y2=y*y ;
        y12=1.-2.*y ;
        b1=2.-y2 ;
        b2=y12*(3.+y2) ;
        b4=y12*y12*y12 ;
        b3=b4+y12*y12 ;
        TotalCrossSection = (1./x-1.)/betasquare+b1*log(x)+b2*(1.-x)-
                          b3*(1.-x2)/2.+b4*(1.-x2*x)/3. ;
      }
    
      TotalCrossSection = twopi_mc2_rcl2 * AtomicNumber
                           *TotalCrossSection/KineticEnergy ;
    }
    else
       TotalCrossSection = 0. ;
 
    return TotalCrossSection ;
}
 
 
G4VParticleChange* G4IeIonisation::PostStepDoIt(
                                              const G4Track& trackData,   
                                              const G4Step& stepData)         
{
  aParticleChange.Initialize(trackData) ;
 
  G4Material*               aMaterial = trackData.GetMaterial() ;
  const G4DynamicParticle*  aParticle = trackData.GetDynamicParticle() ;

  G4double Charge = aParticle->GetDefinition()->GetPDGCharge();
  ParticleMass = aParticle->GetDefinition()->GetPDGMass();
  G4double KineticEnergy = aParticle->GetKineticEnergy();
  G4double TotalEnergy = KineticEnergy + ParticleMass;
  G4double Psquare = KineticEnergy*(TotalEnergy+ParticleMass);
  G4double TotalMomentum = sqrt(Psquare);
  G4double Esquare=TotalEnergy*TotalEnergy;
  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection();

  //  get kinetic energy cut for the electron
  DeltaCutInKineticEnergy = G4Electron::Electron()->GetCutsInEnergy();
  G4double DeltaThreshold = DeltaCutInKineticEnergy[aMaterial->GetIndex()];

  // some kinematics
  G4double MaxKineticEnergyTransfer;
  if (Charge < 0.) MaxKineticEnergyTransfer = 0.5*KineticEnergy;
  else             MaxKineticEnergyTransfer =     KineticEnergy;

  // sampling kinetic energy of the delta ray

  if (MaxKineticEnergyTransfer <= DeltaThreshold) 
     return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  // normal case
  G4double cc,y,y2,c2,b0,b1,b2,b3,b4,x,x1,grej,grejc;

  G4double tau = KineticEnergy/ParticleMass;
  G4double gamma = tau+1., gamma2=gamma*gamma;
  G4double xc = DeltaThreshold/KineticEnergy, xc1=1.-xc;

  if (Charge < 0.)  // Moller (e-e-) scattering
    {
      b1=4./(9.*gamma2-10.*gamma+5.);
      b2=tau*tau*b1; b3=(2.*gamma2+2.*gamma-1.)*b1;
      cc=1.-2.*xc;
      do {
           x    = xc/(1.-cc*G4UniformRand()); x1 = 1.-x;
           grej = b2*x*x-b3*x/x1+b1*gamma2/(x1*x1);
         } while (G4UniformRand()>grej) ;
    }
  else             // Bhabha (e+e-) scattering
    {
      y=1./(gamma+1.); y2=y*y; cc=1.-2.*y;
      b1=2.-y2; b2=cc*(3.+y2);
      c2=cc*cc; b4=c2*cc; b3=c2+b4;
      b0=gamma2/(gamma2-1.);
      grejc=(((b4*xc-b3)*xc+b2)*xc-b1)*xc+b0;
      do {
           x    = xc/(1.-xc1*G4UniformRand());
           grej = ((((b4*x-b3)*x+b2)*x-b1)*x+b0)/grejc;
         } while (G4UniformRand()>grej);
    }

  G4double DeltaKineticEnergy = x * KineticEnergy;

  // protection :do not produce a secondary with 0. kinetic energy !
  if (DeltaKineticEnergy <= 0.)
      return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  G4double DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
                                                       2. * electron_mass_c2 ));
  
  G4double costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2)
                      /(DeltaTotalMomentum * TotalMomentum);

  if (costheta < -1.) costheta = -1.;
  if (costheta > +1.) costheta = +1.;

  //  direction of the delta electron

  G4double phi = twopi * G4UniformRand();
  G4double sintheta = sqrt((1.+costheta)*(1.-costheta));
  G4double dirx = sintheta * cos(phi), diry = sintheta * sin(phi),
                                                       dirz = costheta;

  G4ThreeVector DeltaDirection(dirx,diry,dirz);
  DeltaDirection.rotateUz(ParticleDirection);

  // create G4DynamicParticle object for delta ray

  G4DynamicParticle* theDeltaRay = new G4DynamicParticle;
  theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
  theDeltaRay->SetMomentumDirection(
                   DeltaDirection.x(),DeltaDirection.y(),DeltaDirection.z());
  theDeltaRay->SetDefinition(G4Electron::Electron());
  
  // fill aParticleChange
  // changed energy and momentum of the actual particle
  G4double finalKineticEnergy = KineticEnergy - DeltaKineticEnergy;
 
   if (finalKineticEnergy > 0.)
    {
      G4double finalMomentum=sqrt(finalKineticEnergy*
                         (finalKineticEnergy+2.*ParticleMass));

      G4double finalPx = (TotalMomentum*ParticleDirection.x()
                        - DeltaTotalMomentum*DeltaDirection.x())/finalMomentum;
      G4double finalPy = (TotalMomentum*ParticleDirection.y()
                        - DeltaTotalMomentum*DeltaDirection.y())/finalMomentum;
      G4double finalPz = (TotalMomentum*ParticleDirection.z()
                        - DeltaTotalMomentum*DeltaDirection.z())/finalMomentum;

      aParticleChange.SetMomentumChange( finalPx,finalPy,finalPz );
    }
  else
    {
      finalKineticEnergy = 0.;
      if (Charge < 0.) aParticleChange.SetStatusChange(fStopAndKill);
      else             aParticleChange.SetStatusChange(fStopButAlive);
    }
     
  aParticleChange.SetEnergyChange( finalKineticEnergy );
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary( theDeltaRay );
  aParticleChange.SetLocalEnergyDeposit (0.);

  return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

void G4IeIonisation::PrintInfoDefinition()
{
  G4String comments = "delta cross sections from Moller+Bhabha. ";
           comments += "Good description from 1 KeV to 100 GeV.\n";
           comments += "        delta ray energy sampled from  differential Xsection.";

  G4cout << endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,"Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins. \n";
}

