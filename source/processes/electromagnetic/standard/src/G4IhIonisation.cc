// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IhIonisation.cc,v 1.8 2000-07-31 08:17:49 gracia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4IhIonisation physics process -----------
//                by Laszlo Urban, 30 May 1997 
// **************************************************************
// It is the first implementation of the NEW IONISATION PROCESS.
// It calculates the ionisation of charged hadrons.
// **************************************************************
// corrected by L.Urban on 24/09/97
// several bugs corrected by L.Urban on 13/01/98
// 07-04-98: remove 'tracking cut' of the ionizing particle, MMa
// 29-10-98: small changes , some cleanup, L.Urban
// --------------------------------------------------------------
 

#include "G4IhIonisation.hh"
#include "G4UnitsTable.hh"

// constructor and destructor
 
G4IhIonisation::G4IhIonisation(const G4String& processName)
   : G4VIhEnergyLoss(processName),
     theMeanFreePathTable(NULL),
     theNlambdaTable(NULL),
     theInverseNlambdaTable(NULL),
     theCoeffATable(NULL),
     theCoeffBTable(NULL),
     theCoeffCTable(NULL),
     NumberOfBuildPhysicsTableCalls(0),
     theProton (G4Proton::Proton()),
     theAntiProton (G4AntiProton::AntiProton()),
     theElectron ( G4Electron::Electron() )
{  }

G4IhIonisation::~G4IhIonisation() 
{
     if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
     }

}
 
void G4IhIonisation::SetPhysicsTableBining(G4double lowE, G4double highE,
                                 G4int nBins)
{
  LowestKineticEnergy = lowE;  HighestKineticEnergy = highE;
  TotBin = nBins ;
}

void G4IhIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
//  just call BuildLossTable+BuildLambdaTable
{
 NumberOfBuildPhysicsTableCalls += 1 ;
 if(NumberOfBuildPhysicsTableCalls == 1)
 {
  if(&aParticleType == G4Proton::Proton())
  {
                        LowestKineticEnergy= 1.00*keV;
                        HighestKineticEnergy= 100.*TeV;
                        TotBin=100  ;
  }

    ParticleMass = aParticleType.GetPDGMass() ;

    G4double Charge = aParticleType.GetPDGCharge();

    G4double newCutInRange = aParticleType.GetLengthCuts(); 

  if(Charge>0.)
  {
   if( (CutInRange != newCutInRange) || (theDEDXpTable == NULL))
   {
    BuildLossTable(aParticleType) ;
    RecorderOfpProcess[CounterOfpProcess] = theLossTable ;
    CounterOfpProcess++;
   }
  }
  else
  {
   if( (CutInRange != newCutInRange) || (theDEDXpbarTable == NULL))
   {
    BuildLossTable(aParticleType) ;
    RecorderOfpbarProcess[CounterOfpbarProcess] = theLossTable ;
    CounterOfpbarProcess++;
   }
  }
 
     G4double saveCutInRange = CutInRange ;
     CutInRange = newCutInRange ;
     BuildLambdaTable(aParticleType) ;
     CutInRange = saveCutInRange ;

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

    if(&aParticleType == G4Proton::Proton())
      PrintInfoDefinition();
 }
}

void G4IhIonisation::TestOfInversion(
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
    G4cout << "G4IhIonisation::TestOfInversion (T->Nlambda->Tprime) " << G4endl ;
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

void G4IhIonisation::BuildNlambdaTable(
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
void G4IhIonisation::BuildNlambdaVector(
                               const G4ParticleDefinition& aParticleType,
                                       G4int materialIndex,
                                       G4PhysicsLogVector* nlambdaVector)
{
  G4double LowEdgeEnergy,T,Tlast,dEdx,Value,Vlast,u,du,coeff ;
  G4double Tcut,thresholdEnergy,w1,w2,l ;
  const G4int nbin = 20  ;
  G4bool isOut ;
  const G4double small = 1.e-100;
  const G4double lmin=1.e-100,lmax=1.e100;

  const G4MaterialTable* theMaterialTable=
                          G4Material::GetMaterialTable();

  //  the threshold energy for delta production
    Tcut =  DeltaCutInKineticEnergy[materialIndex];
    w1 = electron_mass_c2*(2.*ParticleMass-Tcut) ;
    w2 = electron_mass_c2+ParticleMass ;
    thresholdEnergy = 0.5*(sqrt(w1*w1+2.*electron_mass_c2*Tcut*w2*w2)-w1)
                      /electron_mass_c2 ;
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

      l = (*theMeanFreePathTable)[materialIndex]->GetValue(T,isOut);
      if((l>lmin) && (l<lmax))
        Value += coeff*T/(G4EnergyLossTables::GetPreciseDEDX(&aParticleType,
                                 T,(*theMaterialTable)[materialIndex])*l) ;
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
void G4IhIonisation::BuildCoeffATable(
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

void G4IhIonisation::BuildCoeffBTable(
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
void G4IhIonisation::BuildCoeffCTable(
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
void G4IhIonisation::BuildInverseNlambdaTable(
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

void G4IhIonisation::InvertNlambdaVector(
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

void G4IhIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
// Build tables for the ionization energy loss
//  the tables are built for MATERIALS
//                           *********
    G4double Charge = aParticleType.GetPDGCharge() ;

 // cuts for p/pbar and electron ....................
   if(Charge>0.)
    ParticleCutInKineticEnergy = theProton->GetCutsInEnergy() ;
   else
    ParticleCutInKineticEnergy = theAntiProton->GetCutsInEnergy() ;

    DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;

  G4double LowEdgeEnergy , ionloss ;
  G4double  RateMass ;
  G4double deltaloss ;

  G4bool isOutRange ;
  static const G4MaterialTable* theMaterialTable=
                                   G4Material::GetMaterialTable();
  const G4double twoln10 = 2.*log(10.) ;
  const G4double Factor = twopi_mc2_rcl2 ;
  const G4double bg2lim = 0.0169 , taulim = 8.4146e-3 ;

  RateMass = electron_mass_c2/proton_mass_c2 ;

  //  create table

  G4int numOfMaterials = theMaterialTable->length();

  if ( theLossTable) {
     theLossTable->clearAndDestroy();
     delete theLossTable;
  }
  theLossTable = new G4PhysicsTable(numOfMaterials);

  //  loop for materials

  for (G4int J=0; J<numOfMaterials; J++)
  {

    // create physics vector and fill it

    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                    LowestKineticEnergy, HighestKineticEnergy, TotBin);
  
    // get material parameters needed for the energy loss calculation

    G4double ElectronDensity,Eexc,Eexc2,Cden,Mden,Aden,X0den,X1den,taul ;
    G4double* ShellCorrectionVector;
   
    const G4Material* material= (*theMaterialTable)[J];

    ElectronDensity = material->GetElectronDensity();
    Eexc = material->GetIonisation()->GetMeanExcitationEnergy();
    Eexc2 = Eexc*Eexc ;
    Cden = material->GetIonisation()->GetCdensity();
    Mden = material->GetIonisation()->GetMdensity();
    Aden = material->GetIonisation()->GetAdensity();
    X0den = material->GetIonisation()->GetX0density();
    X1den = material->GetIonisation()->GetX1density();
    taul = material->GetIonisation()->GetTaul() ;
    ShellCorrectionVector =
                     material->GetIonisation()->GetShellCorrectionVector();

    // get elements in the actual material,
    // they are needed for the low energy part ....

    const G4ElementVector* theElementVector=
                   material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector=
                   material->GetAtomicNumDensityVector() ;
    const G4int NumberOfElements=
                   material->GetNumberOfElements() ;
 
    // get  electron cut in kin. energy for the material

    DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;

    // some local variables -------------------
    G4double tau,tau0,Tmax,gamma,bg2,beta2,rcut,delta,x,sh ;

    for (G4int i = 0 ; i < TotBin ; i++)
    {
      LowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
      tau = LowEdgeEnergy/proton_mass_c2 ;
      gamma = tau +1. ;
      bg2 = tau*(tau+2.) ;
      beta2 = bg2/(gamma*gamma) ;
      Tmax = 2.*electron_mass_c2*bg2
             /(1.+2.*gamma*RateMass+RateMass*RateMass) ;

      if ( tau < taul )
      //  low energy part , parametrized energy loss formulae
      {
        ionloss = 0. ;
        deltaloss = 0. ;

        //  loop for the elements in the material
        for (G4int iel=0; iel<NumberOfElements; iel++)
        {
          const G4Element* element = (*theElementVector)(iel);
          
          if ( tau < element->GetIonisation()->GetTau0())  
            ionloss += theAtomicNumDensityVector[iel]
                       *( element->GetIonisation()->GetAlow()*sqrt(tau)
                       +element->GetIonisation()->GetBlow()*tau) ;
          else
            ionloss += theAtomicNumDensityVector[iel]
                       *  element->GetIonisation()->GetClow()/sqrt(tau) ;
        }
        if ( DeltaCutInKineticEnergyNow < Tmax)
        {
          deltaloss = log(Tmax/DeltaCutInKineticEnergyNow)-
                      beta2*(1.-DeltaCutInKineticEnergyNow/Tmax) ;
          if(aParticleType.GetPDGSpin() == 0.5)
            deltaloss += 0.25*(Tmax-DeltaCutInKineticEnergyNow)*
                              (Tmax-DeltaCutInKineticEnergyNow)/
                        (LowEdgeEnergy*LowEdgeEnergy+proton_mass_c2*proton_mass_c2) ;
            deltaloss *= Factor*ElectronDensity/beta2 ;
        }
        ionloss -= deltaloss ;
      }
      else
      // high energy part , Bethe-Bloch formula 
      {
        if ( DeltaCutInKineticEnergyNow < Tmax)
          rcut = DeltaCutInKineticEnergyNow/Tmax ;
        else
          rcut = 1.;

        ionloss = log(2.*electron_mass_c2*bg2*Tmax/Eexc2)
                  +log(rcut)-(1.+rcut)*beta2 ;


        // density correction 
        x = log(bg2)/twoln10 ;
        if ( x < X0den )
          delta = 0. ;
        else 
        {
          delta = twoln10*x - Cden ;
          if ( x < X1den )
            delta += Aden*pow((X1den-x),Mden) ;
        } 

        // shell correction 
        if ( bg2 > bg2lim ) {
          sh = 0. ;      
          x = 1. ;
          for (G4int k=0; k<=2; k++) {
            x *= bg2 ;
            sh += ShellCorrectionVector[k]/x;
          }
        }
        else {
          sh = 0. ;      
          x = 1. ;
          for (G4int k=0; k<=2; k++) {
             x *= bg2lim ;
             sh += ShellCorrectionVector[k]/x;
          }
          sh *= log(tau/taul)/log(taulim/taul) ;     
        }

        // now you can compute the total ionization loss
        ionloss -= delta + sh ;
        ionloss *= Factor*ElectronDensity/beta2 ;
      }
      if ( ionloss <= 0.)
        ionloss = 0. ;
      aVector->PutValue(i,ionloss) ;

    }
    theLossTable->insert(aVector);
  }

}

void G4IhIonisation::BuildLambdaTable(const G4ParticleDefinition& aParticleType)
{
  // Build mean free path tables for the delta ray production process
  //     tables are built for MATERIALS 

      G4double LowEdgeEnergy , Value ,sigma ;
      G4bool isOutRange ;
      const G4MaterialTable* theMaterialTable=
                                         G4Material::GetMaterialTable();
      //create table

      G4int numOfMaterials = theMaterialTable->length();


      if (theMeanFreePathTable) {
        theMeanFreePathTable->clearAndDestroy();
        delete theMeanFreePathTable;
      }

      theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);

  // get electron and particle cuts in kinetic energy

      DeltaCutInKineticEnergy = theElectron->GetCutsInEnergy() ;
 
      ParticleCutInKineticEnergy = aParticleType.GetEnergyCuts() ;

  // loop for materials 

      for (G4int J=0 ; J < numOfMaterials; J++)
      { 
        //create physics vector then fill it ....

        G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
               LowestKineticEnergy, HighestKineticEnergy, TotBin);

  // compute the (macroscopic) cross section first
 
        const G4Material* material= (*theMaterialTable)[J];
        
        const G4ElementVector* theElementVector=
                         material->GetElementVector() ;
        const G4double* theAtomicNumDensityVector =
                         material->GetAtomicNumDensityVector();
        const G4int NumberOfElements=
                         material->GetNumberOfElements() ;
        DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[J] ;

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


G4double G4IhIonisation::ComputeMicroscopicCrossSection(
                                 const G4ParticleDefinition& aParticleType,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber)
{
  //******************************************************************
  // cross section formula is OK for spin=0 and 1/2 only !
  // *****************************************************************

  // calculates the microscopic cross section in GEANT4 internal units
  //    ( it is called for elements , AtomicNumber = Z )

    G4double TotalEnergy,
             betasquare,
             MaxKineticEnergyTransfer,TotalCrossSection,tempvar;

    // get particle data ...................................

    TotalEnergy=KineticEnergy + ParticleMass;

    // some kinematics......................

    betasquare = KineticEnergy*(TotalEnergy+ParticleMass)
                 /(TotalEnergy*TotalEnergy);
    tempvar = ParticleMass+electron_mass_c2;
    MaxKineticEnergyTransfer = 2.*electron_mass_c2*KineticEnergy
                     *(TotalEnergy+ParticleMass)
                     /(tempvar*tempvar+2.*electron_mass_c2*KineticEnergy);

    // now you can calculate the total cross section ------------------

    if( MaxKineticEnergyTransfer > DeltaCutInKineticEnergyNow )
    {
       tempvar=DeltaCutInKineticEnergyNow/MaxKineticEnergyTransfer;
       TotalCrossSection = (1.-tempvar*(1.-betasquare*log(tempvar)))
                           /DeltaCutInKineticEnergyNow;
  // +term for spin=1/2 particle
     if(aParticleType.GetPDGSpin() == 1)
     {
       TotalCrossSection +=  0.5
                       *(MaxKineticEnergyTransfer-DeltaCutInKineticEnergyNow)
                       /(TotalEnergy*TotalEnergy);
     }
       TotalCrossSection = twopi_mc2_rcl2 * AtomicNumber
                           *TotalCrossSection/betasquare;
    }
    else
       TotalCrossSection= 0. ;
 
    return TotalCrossSection ;
}
 
G4VParticleChange* G4IhIonisation::PostStepDoIt(const G4Track& trackData,
                                                const G4Step& stepData)         
{
  const G4DynamicParticle* aParticle ;
  G4Material* aMaterial;
  G4double KineticEnergy,TotalEnergy,TotalMomentum,
           betasquare,MaxKineticEnergyTransfer,
           DeltaKineticEnergy,DeltaTotalMomentum,costheta,sintheta,phi,
           dirx,diry,dirz,finalKineticEnergy,finalPx,finalPy,finalPz,
           x,xc,te2,grej,Psquare,Esquare,summass,rate,grejc,finalMomentum ;

  aParticleChange.Initialize(trackData) ;
  aMaterial = trackData.GetMaterial() ;

  aParticle = trackData.GetDynamicParticle() ;

  G4double Charge=aParticle->GetDefinition()->GetPDGCharge();
  KineticEnergy=aParticle->GetKineticEnergy();
  TotalEnergy=KineticEnergy + ParticleMass ;
  Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ;
  Esquare=TotalEnergy*TotalEnergy ;
  summass = ParticleMass + electron_mass_c2 ;    
  G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection() ;

  //  get kinetic energy cut for particle and for the electron....
  ParticleCutInKineticEnergyNow =
                       ParticleCutInKineticEnergy[aMaterial->GetIndex()];
  DeltaCutInKineticEnergyNow = DeltaCutInKineticEnergy[aMaterial->GetIndex()];

  // some kinematics......................

  betasquare=Psquare/Esquare ;
  MaxKineticEnergyTransfer = 2.*electron_mass_c2*Psquare
                      /(summass*summass+2.*electron_mass_c2*KineticEnergy);

  // sampling kinetic energy of the delta ray 

  if( MaxKineticEnergyTransfer <= DeltaCutInKineticEnergyNow )
  {
    // pathological case (it should not happen ,
    // there is no change at all).....

    // return &aParticleChange;
    return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
  }
  else
  {
   // normal case ......................................
      xc=DeltaCutInKineticEnergyNow/MaxKineticEnergyTransfer ;
      rate=MaxKineticEnergyTransfer/TotalEnergy ;

     if(aParticle->GetDefinition()->GetPDGSpin() == 1)     
      te2=0.5*rate*rate ;
     else
      te2=0. ;

   // sampling follows ...
        grejc=1.-betasquare*xc+te2*xc*xc ;
      do {
        x=xc/(1.-(1.-xc)*G4UniformRand());
        grej=(1.-x*(betasquare-x*te2))/grejc ;
      } while( G4UniformRand()>grej );
   }
   
   DeltaKineticEnergy = x * MaxKineticEnergyTransfer ;
   if(DeltaKineticEnergy <= 0.)
     return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

   DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
                                               2. * electron_mass_c2 )) ;
   TotalMomentum = sqrt(Psquare) ;
   costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2)
            /(DeltaTotalMomentum * TotalMomentum) ;

   //  protection against costheta > 1 or < -1   ---------------
   if ( costheta < -1. ) 
          costheta = -1. ;
   if ( costheta > +1. ) 
          costheta = +1. ;

   //  direction of the delta electron  ........
   phi = twopi * G4UniformRand() ; 
   sintheta = sqrt((1.+costheta)*(1.-costheta));
   dirx = sintheta * cos(phi) ;
   diry = sintheta * sin(phi) ;
   dirz = costheta ;

   G4ThreeVector DeltaDirection(dirx,diry,dirz) ;
   DeltaDirection.rotateUz(ParticleDirection) ;

   // create G4DynamicParticle object for delta ray
   G4DynamicParticle *theDeltaRay = new G4DynamicParticle;
   theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
   theDeltaRay->SetMomentumDirection(
                   DeltaDirection.x(),DeltaDirection.y(),DeltaDirection.z()); 
   theDeltaRay->SetDefinition(G4Electron::Electron());

   // fill aParticleChange 
   finalKineticEnergy = KineticEnergy - DeltaKineticEnergy ;
   if (finalKineticEnergy > 0.)
     {
      finalPx = TotalMomentum*ParticleDirection.x()
                        - DeltaTotalMomentum*DeltaDirection.x();
      finalPy = TotalMomentum*ParticleDirection.y()
                        - DeltaTotalMomentum*DeltaDirection.y();
      finalPz = TotalMomentum*ParticleDirection.z()
                        - DeltaTotalMomentum*DeltaDirection.z();
      finalMomentum =
                sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz) ;
      finalPx /= finalMomentum ;
      finalPy /= finalMomentum ;
      finalPz /= finalMomentum ;

      aParticleChange.SetMomentumChange( finalPx,finalPy,finalPz );
     }
   else
     {
       finalKineticEnergy = 0. ;
       if (aParticle->GetDefinition()->GetParticleName() == "proton")
             aParticleChange.SetStatusChange(fStopAndKill);
       else  aParticleChange.SetStatusChange(fStopButAlive);
     }

   aParticleChange.SetEnergyChange( finalKineticEnergy );
   aParticleChange.SetNumberOfSecondaries(1);   
   aParticleChange.AddSecondary( theDeltaRay );
   aParticleChange.SetLocalEnergyDeposit (0.);
      
  return G4IVContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}
void G4IhIonisation::PrintInfoDefinition()
{
  G4String comments = "  Knock-on electron cross sections . ";
           comments += "\n         Good description above the mean excitation energy.\n";
           comments += "         delta ray energy sampled from  differential Xsection.";

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
                                                  "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins. \n";
}


