// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuPairProduction.cc,v 1.10 2000-02-10 09:10:39 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------- G4MuPairProduction physics process ---------
//                by Laszlo Urban, May 1998 
// **************************************************************
// 04-06-98, in DoIt,secondary production condition:range>G4std::min(threshold,safety)
// 26/10/98, new stuff from R. Kokoulin + cleanup , L.Urban
// 06/05/99 , bug fixed , L.Urban
// 10/02/00  modifications+bug fix , new e.m. structure, L.Urban
// --------------------------------------------------------------

#include "G4MuPairProduction.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"

// static members ........

G4int G4MuPairProduction::nzdat = 5 ;
G4double G4MuPairProduction::zdat[]={1.,4.,13.,26.,92.};
G4int G4MuPairProduction::ntdat = 8 ;
G4double G4MuPairProduction::tdat[]={1.e3,1.e4,1.e5,1.e6,1.e7,1.e8,1.e9,1.e10};
G4int G4MuPairProduction::NBIN = 100 ; //500 ;
G4double G4MuPairProduction::ya[1000]={0.};
G4double G4MuPairProduction::proba[5][8][1000]={0.};
 
G4MuPairProduction::G4MuPairProduction(const G4String& processName)
  : G4MuEnergyLoss(processName),  
    theMeanFreePathTable(NULL),
    LowerBoundLambda(1.*keV),
    UpperBoundLambda(10000.*TeV),
    NbinLambda(100)
{  }
 
 
G4MuPairProduction::~G4MuPairProduction()
{
   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }

   if (&PartialSumSigma) {
      PartialSumSigma.clearAndDestroy();
   }
}

void G4MuPairProduction::BuildPhysicsTable(
                               const G4ParticleDefinition& aParticleType)
//  just call BuildLossTable+BuildLambdaTable
{
 
    // get bining from EnergyLoss
    LowestKineticEnergy  = GetLowerBoundEloss() ;
    HighestKineticEnergy = GetUpperBoundEloss() ;
    TotBin               = GetNbinEloss() ;

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

  // sampling table should be made only once !
  if(theMeanFreePathTable == NULL)
     MakeSamplingTables(&aParticleType) ;

  G4double electronCutInRange = G4Electron::Electron()->GetCuts();
  if(electronCutInRange != lastelectronCutInRange)
    BuildLambdaTable(aParticleType) ;

  G4MuEnergyLoss::BuildDEDXTable(aParticleType) ;

  if(&aParticleType==theMuonPlus)
    PrintInfoDefinition() ;
}

void G4MuPairProduction::BuildLossTable(
                                   const G4ParticleDefinition& aParticleType)
{
  G4double KineticEnergy,TotalEnergy,pairloss,Z,
           loss,natom,eCut,pCut ;

  const G4MaterialTable* theMaterialTable =
                                G4Material::GetMaterialTable();
  ParticleMass = aParticleType.GetPDGMass() ;
  ElectronCutInKineticEnergy = (*theElectron).GetEnergyCuts() ;
  PositronCutInKineticEnergy = (*thePositron).GetEnergyCuts() ;

  G4int numOfMaterials = theMaterialTable->length() ;

  if (theLossTable) {
      theLossTable->clearAndDestroy();
      delete theLossTable;
    }

  theLossTable = new G4PhysicsTable(numOfMaterials) ;

  for (G4int J=0; J<numOfMaterials; J++)
  {
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                               LowestKineticEnergy,HighestKineticEnergy,TotBin);
    ElectronCutInKineticEnergyNow = ElectronCutInKineticEnergy[J] ;
    PositronCutInKineticEnergyNow = PositronCutInKineticEnergy[J] ;
    const G4Material* material = (*theMaterialTable)[J] ;
    const G4ElementVector* theElementVector =
                                     material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector =
                                     material->GetAtomicNumDensityVector() ;
    const G4int NumberOfElements =
                                     material->GetNumberOfElements() ;

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
      for (G4int iel=0; iel<NumberOfElements; iel++)
      {
        Z=(*theElementVector)(iel)->GetZ();
        natom = theAtomicNumDensityVector[iel] ;
        loss = ComputePairLoss(&aParticleType,
                                 Z,KineticEnergy,eCut,pCut) ;   
        pairloss += natom*loss ;
      }
      if(pairloss<0.)
        pairloss = 0. ;
      aVector->PutValue(i,pairloss);
    }

    theLossTable->insert(aVector);
  }
}

G4double G4MuPairProduction::ComputePairLoss(
                                     const G4ParticleDefinition* ParticleType,
                                             G4double AtomicNumber,
                                             G4double KineticEnergy,
                                             G4double ElectronEnergyCut, 
                                             G4double PositronEnergyCut) 
{
  static const G4double
  xgi[] ={ 0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801 };
  static const G4double
  wgi[] ={ 0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506 };
  static const G4double ak1=6.9 ;
  static const G4double ak2=1.0 ;
  static const G4double sqrte = sqrt(exp(1.)) ;
  G4double z13 = exp(log(AtomicNumber)/3.) ;

  G4double loss = 0.0 ;

  if ( AtomicNumber < 1. ) return loss;

  G4double CutInPairEnergy = ElectronEnergyCut+PositronEnergyCut
                            +2.*electron_mass_c2 ;
  G4double MinPairEnergy = 4.*electron_mass_c2 ;
  if( CutInPairEnergy <= MinPairEnergy ) return loss ;

  G4double MaxPairEnergy = KineticEnergy+ParticleMass*(1.-0.75*sqrte*z13) ;
  if(MaxPairEnergy < MinPairEnergy)
     MaxPairEnergy = MinPairEnergy ;
     
  if( CutInPairEnergy >= MaxPairEnergy ) 
      CutInPairEnergy = MaxPairEnergy ;

  if(MaxPairEnergy <= MinPairEnergy) return loss ;

  G4double aaa,bbb,hhh,x,epln,ep ;
  G4int kkk ;
 // calculate the rectricted loss    
 // numerical integration in log(PairEnergy)
  aaa = log(MinPairEnergy) ;
  bbb = log(CutInPairEnergy) ;
  kkk = int((bbb-aaa)/ak1+ak2) ;
  hhh = (bbb-aaa)/kkk ;
 
  for (G4int l=0 ; l<kkk; l++)
  {
    x = aaa+hhh*l ;
    for (G4int ll=0; ll<8; ll++)
    {
      epln=x+xgi[ll]*hhh ;
      ep = exp(epln) ;
      loss += wgi[ll]*ep*ep*ComputeDMicroscopicCrossSection(ParticleType,
                                             KineticEnergy,AtomicNumber,
                                             ep) ;
    }
  }
  loss *= hhh ;
  if (loss < 0.) loss = 0.;
  return loss ;
}
 

void G4MuPairProduction::BuildLambdaTable(
                                   const G4ParticleDefinition& ParticleType)
{
  G4double LowEdgeEnergy , Value;
  G4double FixedEnergy = (LowestKineticEnergy + HighestKineticEnergy)/2. ;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

  if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }
  theMeanFreePathTable = new 
                      G4PhysicsTable( G4Material::GetNumberOfMaterials() ) ;
   PartialSumSigma.resize(G4Material::GetNumberOfMaterials());

   G4PhysicsLogVector* ptrVector;
   for ( G4int J=0 ; J < G4Material::GetNumberOfMaterials(); J++ )  
   { 
      ptrVector = new G4PhysicsLogVector(
               LowerBoundLambda,UpperBoundLambda,NbinLambda);


     const G4Material* material= (*theMaterialTable)[J];

     for ( G4int i = 0 ; i < NbinLambda ; i++ )      
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

void G4MuPairProduction::ComputePartialSumSigma(
                                    const G4ParticleDefinition* ParticleType,
                                               G4double KineticEnergy,
                                               const G4Material* aMaterial)
{
  G4int Imate = aMaterial->GetIndex();
  G4int NbOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector(); 
  const G4double* theAtomNumDensityVector = aMaterial->
                                                GetAtomicNumDensityVector();
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

G4double G4MuPairProduction::ComputeMicroscopicCrossSection(
                                     const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy,
                                           G4double AtomicNumber,
                                           G4double ElectronEnergyCut,
                                           G4double PositronEnergyCut)
 
{
  static const G4double
  xgi[] ={ 0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801 };
  static const G4double
  wgi[] ={ 0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506 };
  static const G4double ak1=6.9 ;
  static const G4double ak2=1.0 ;

  static const G4double sqrte = sqrt(exp(1.)) ;
  G4double z13 = exp(log(AtomicNumber)/3.) ;

  G4double CrossSection = 0.0 ;

  if ( AtomicNumber < 1. ) return CrossSection;

  G4double CutInPairEnergy = ElectronEnergyCut+PositronEnergyCut
                            +2.*electron_mass_c2 ;

  if( CutInPairEnergy < 4.*electron_mass_c2 )
    CutInPairEnergy = 4.*electron_mass_c2 ;

  G4double MaxPairEnergy = KineticEnergy+ParticleMass*(1.-0.75*sqrte*z13) ;
  if(MaxPairEnergy < CutInPairEnergy)
     MaxPairEnergy = CutInPairEnergy ;
  if( CutInPairEnergy >= MaxPairEnergy ) return CrossSection ;

  G4double aaa,bbb,hhh,x,epln,ep ;
  G4int kkk ;
 // calculate the total cross section
 // numerical integration in log(PairEnergy)
  aaa = log(CutInPairEnergy) ;
  bbb = log(MaxPairEnergy) ;
  kkk = int((bbb-aaa)/ak1+ak2) ;
  hhh = (bbb-aaa)/kkk ;
  for (G4int l=0 ; l<kkk; l++)
  {
    x = aaa+hhh*l ;
    for (G4int ll=0; ll<8; ll++)
    {
      epln = x+xgi[ll]*hhh;
      ep = exp(epln) ;
      CrossSection += wgi[ll]*ep*ComputeDMicroscopicCrossSection(ParticleType,
                                                 KineticEnergy,AtomicNumber,
                                                 ep) ;
    }
  }
  CrossSection *= hhh ;

  if (CrossSection < 0.) CrossSection = 0.;

  return CrossSection;
}

void G4MuPairProduction::MakeSamplingTables(
                                   const G4ParticleDefinition* ParticleType)
{
  G4double epbin[1000],xbin[1000],prbin[1000] ;
  G4int nbin;
  G4double AtomicNumber,KineticEnergy,MinPairEnergy ;  
  G4double c,y,ymin,ymax,dy,yy,dx,x,ep ;

  MinPairEnergy = 4.*electron_mass_c2 ;

  static const G4double sqrte = sqrt(exp(1.)) ;

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
                                          KineticEnergy,AtomicNumber,ep);
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
      ya[NBIN]=0. ;

      if(CrossSection > 0.)
      {
        for(G4int ib=0; ib<=nbin; ib++)
        {
          prbin[ib] /= CrossSection ;
          proba[iz][it][ib] /= CrossSection ;
        }
      }
    }
  }
} 
 

G4double G4MuPairProduction::ComputeDDMicroscopicCrossSection(
                                 const G4ParticleDefinition* ParticleType,
                                 G4double KineticEnergy, G4double AtomicNumber,
                                 G4double PairEnergy,G4double asymmetry)
 // Calculates the double differential (DD) microscopic cross section 
 //   using the cross section formula of R.P. Kokoulin (18/01/98)
{
  static const G4double sqrte = sqrt(exp(1.)) ;

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
      
G4double G4MuPairProduction::ComputeDMicroscopicCrossSection(
                                 const G4ParticleDefinition* ParticleType,
                                 G4double KineticEnergy, G4double AtomicNumber,
                                 G4double PairEnergy)
 // Calculates the  differential (D) microscopic cross section 
 //   using the cross section formula of R.P. Kokoulin (18/01/98)
{

  static const G4double
  xgi[] ={ 0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801 };

  static const G4double
  wgi[] ={ 0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506 };

  G4double TotalEnergy = KineticEnergy + ParticleMass ;
  G4double EnergyLoss = TotalEnergy - PairEnergy ;
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
      
G4VParticleChange* G4MuPairProduction::PostStepDoIt(const G4Track& trackData,
                                                  const G4Step& stepData)
{
   static const G4double esq = sqrt(exp(1.));

   aParticleChange.Initialize(trackData);
   G4Material* aMaterial=trackData.GetMaterial() ;
   const G4DynamicParticle* aDynamicParticle=trackData.GetDynamicParticle();  
   G4double           KineticEnergy     = aDynamicParticle->GetKineticEnergy();
   G4ParticleMomentum ParticleDirection = 
                                      aDynamicParticle->GetMomentumDirection();

   // e-e+ cut in this material
   G4double ElectronEnergyCut = electron_mass_c2+
      ((*G4Electron::Electron()).GetCutsInEnergy())[aMaterial->GetIndex()];
   G4double PositronEnergyCut = electron_mass_c2+
      ((*G4Positron::Positron()).GetCutsInEnergy())[aMaterial->GetIndex()];
   G4double CutInPairEnergy = ElectronEnergyCut + PositronEnergyCut ;

   G4double MinPairEnergy = 4.*electron_mass_c2 ;
   if (CutInPairEnergy < MinPairEnergy) CutInPairEnergy = MinPairEnergy ;

   // check against insufficient energy
   if(KineticEnergy < CutInPairEnergy )
     return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

   // select randomly one element constituing the material  
   G4Element* anElement = SelectRandomAtom(aMaterial);

   // limits of the energy sampling
   G4double TotalEnergy = KineticEnergy + ParticleMass ;
   G4double TotalMomentum = sqrt(KineticEnergy*(TotalEnergy+ParticleMass)) ;
   G4double Z3 = anElement->GetIonisation()->GetZ3() ;
   G4double MaxPairEnergy = TotalEnergy-0.75*esq*ParticleMass*Z3 ;

   if(MinPairEnergy >= MaxPairEnergy)
     return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

   // sample e-e+ energy, pair energy first
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
       delmin=del;
       itt=it ;
     }
   }

   xc = log(CutInPairEnergy/MinPairEnergy)/log(MaxPairEnergy/MinPairEnergy) ;
   yc = log(xc) ;
   
   iy = -1 ;
   do {
       iy += 1 ;
      } while ((ya[iy] < yc )&&(iy < NBINminus1)) ;
   G4double norm = proba[izz][itt][iy] ;

   G4double r = norm+G4UniformRand()*(1.-norm) ;
 
   iy -= 1 ;
   do {
        iy += 1 ;
      } while ((proba[izz][itt][iy] < r)&&(iy < NBINminus1)) ;

   //sampling is uniformly in y in the bin
   if( iy < NBIN )
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
   //      (mean theta for the moment)
   G4double Teta = electron_mass_c2/TotalEnergy ;

   G4double Phi  = twopi * G4UniformRand() ;
   G4double dirx = sin(Teta)*cos(Phi) , diry = sin(Teta)*sin(Phi) ,
            dirz = cos(Teta) ;

   G4double LocalEnerDeposit = 0. ;
   G4int numberofsecondaries = 1 ;
   G4int flagelectron = 0 ;
   G4int flagpositron = 1 ; 
   G4DynamicParticle *aParticle1,*aParticle2 ;
   G4double ElectronMomentum , PositronMomentum ;
   G4double finalPx,finalPy,finalPz ;

   G4double ElectKineEnergy = ElectronEnergy - electron_mass_c2 ;

   if((ElectKineEnergy > ElectronEnergyCut) ||
      (G4EnergyLossTables::GetRange(
             G4Electron::Electron(),ElectKineEnergy,aMaterial) >=
                          stepData.GetPostStepPoint()->GetSafety()))
      {
        numberofsecondaries += 1 ;
        flagelectron = 1 ;
        ElectronMomentum = sqrt(ElectKineEnergy*
                                          (ElectronEnergy+electron_mass_c2));
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

   if((PositKineEnergy < PositronEnergyCut) &&
      (G4EnergyLossTables::GetRange(
             G4Positron::Positron(),PositKineEnergy,aMaterial) <=
                          stepData.GetPostStepPoint()->GetSafety()))
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

   G4double KinEnergyCut = (aDynamicParticle->GetDefinition()->
                            GetEnergyCuts())[aMaterial->GetIndex()];

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

   return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

G4Element* G4MuPairProduction::SelectRandomAtom(G4Material* aMaterial) const
{
  // select randomly 1 element within the material

  const G4int Index = aMaterial->GetIndex();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  G4double rval = G4UniformRand()*((*PartialSumSigma(Index))
                                                      (NumberOfElements-1));

  for ( G4int i=0; i < NumberOfElements; i++ )
  {
    if (rval <= (*PartialSumSigma(Index))(i)) return ((*theElementVector)(i));
  }
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements, NULL pointer returned." << G4endl;
  return NULL;
}
void G4MuPairProduction::PrintInfoDefinition()
{
  G4String comments = "theoretical cross sections \n ";
           comments += "         Good description up to 1000 TeV.";

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n    PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
                                                     "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins. \n";
}


