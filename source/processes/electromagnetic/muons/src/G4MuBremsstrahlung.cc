// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuBremsstrahlung.cc,v 1.14 2001-02-05 17:50:29 gcosmo Exp $
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
//      -------- G4MuBremsstrahlung physics process ---------
//                by Laszlo Urban, September 1997
//
// 08-04-98: remove 'tracking cut' of muon in oIt, MMa
// 26/10/98: new cross section of R.Kokoulin,cleanup , L.Urban
// 10/02/00  modifications , new e.m. structure, L.Urban
// --------------------------------------------------------------
#include "G4MuBremsstrahlung.hh"
#include "G4UnitsTable.hh"
 
// static members ........
G4int G4MuBremsstrahlung::nzdat =  5 ;
G4double G4MuBremsstrahlung::zdat[]={1.,4.,13.,29.,92.};
G4double G4MuBremsstrahlung::adat[]={1.01,9.01,26.98,63.55,238.03};
G4int G4MuBremsstrahlung::ntdat = 8 ;
G4double G4MuBremsstrahlung::tdat[]={1.e3,1.e4,1.e5,1.e6,1.e7,1.e8,1.e9,1.e10};
G4int G4MuBremsstrahlung::NBIN = 1000;    // 100 ;
G4double G4MuBremsstrahlung::ya[1001]={0.};
G4double G4MuBremsstrahlung::proba[5][8][1001]={0.};
G4double G4MuBremsstrahlung::CutFixed=0.98*keV ;


G4double G4MuBremsstrahlung::LowerBoundLambda = 1.*keV ;
G4double G4MuBremsstrahlung::UpperBoundLambda = 1000000.*TeV ;
G4int	 G4MuBremsstrahlung::NbinLambda = 150 ;


// constructor
 
G4MuBremsstrahlung::G4MuBremsstrahlung(const G4String& processName)
  : G4VMuEnergyLoss(processName),  
    theMeanFreePathTable(NULL)
{  }
 
G4MuBremsstrahlung::~G4MuBremsstrahlung()
{
   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();

      delete theMeanFreePathTable;
   }
   if (&PartialSumSigma) {
      PartialSumSigma.clearAndDestroy();
   }
}
 
void G4MuBremsstrahlung::BuildPhysicsTable(
                                 const G4ParticleDefinition& aParticleType)
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

   if(theMeanFreePathTable == NULL)
     MakeSamplingTables(&aParticleType) ;

  G4double gammaCutInRange = G4Gamma::Gamma()->GetCuts();
  if(gammaCutInRange != lastgammaCutInRange)
    BuildLambdaTable(aParticleType) ;
  
  G4VMuEnergyLoss::BuildDEDXTable(aParticleType) ;

  if(&aParticleType == theMuonPlus)
    PrintInfoDefinition() ;
}


void G4MuBremsstrahlung::BuildLossTable(
                              const G4ParticleDefinition& aParticleType)
{
  G4double KineticEnergy,TotalEnergy,bremloss,Z,
           loss,natom,Cut ;

  const G4MaterialTable* theMaterialTable =
                                G4Material::GetMaterialTable();
  ParticleMass = aParticleType.GetPDGMass();
  GammaCutInKineticEnergy = G4Gamma::Gamma()->GetEnergyCuts() ;

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

    GammaCutInKineticEnergyNow = GammaCutInKineticEnergy[J] ;
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
      Cut = GammaCutInKineticEnergyNow ;
      if(Cut>KineticEnergy) Cut = KineticEnergy ;
      bremloss = 0.;      

      for (G4int iel=0; iel<NumberOfElements; iel++)
      {
        Z=(*theElementVector)(iel)->GetZ();
        natom = theAtomicNumDensityVector[iel] ;
        loss = ComputeBremLoss((&aParticleType),Z,
                                        (*theElementVector)(iel)->GetA(),
                                         KineticEnergy,Cut) ;   
        bremloss += natom*loss ;
      } 
      if(bremloss<0.) bremloss = 0. ; 
      
      aVector->PutValue(i,bremloss);
    }
    theLossTable->insert(aVector);
  }
}

G4double G4MuBremsstrahlung::ComputeBremLoss(
                              const G4ParticleDefinition* aParticleType,
                              G4double AtomicNumber,G4double AtomicMass,
                              G4double KineticEnergy,G4double GammaEnergyCut)
{
  G4double TotalEnergy,vcut,vmax,aaa,bbb,hhh,aa,x,ep ;
  G4int kkk ;
  G4double ak1=0.05 ;
  G4int k2=5 ;
  G4double xgi[]={0.03377,0.16940,0.38069,0.61931,0.83060,0.96623};
  G4double wgi[]={0.08566,0.18038,0.23396,0.23396,0.18038,0.08566};
  G4double loss = 0. ;

  TotalEnergy=KineticEnergy+ParticleMass ;
  vcut = GammaEnergyCut/TotalEnergy ;
  vmax = KineticEnergy/TotalEnergy ;

  aaa=0.;
  bbb=vcut ;
  if(vcut>vmax) bbb=vmax ; 
  kkk=int((bbb-aaa)/ak1)+k2 ;
  hhh=(bbb-aaa)/float(kkk) ;

  for(G4int l=0; l<kkk; l++)
  {
    aa=aaa+hhh*float(l) ;
    for(G4int i=0; i<6; i++)
    {
      x=aa+xgi[i]*hhh ;
      ep=x*TotalEnergy ;
      loss += ep*wgi[i]*ComputeDMicroscopicCrossSection(
                                       aParticleType,KineticEnergy,
                                       AtomicNumber,AtomicMass,ep) ;
    }
  }

  loss *=hhh*TotalEnergy ;

  return loss ;
}

void G4MuBremsstrahlung::BuildLambdaTable(
                                  const G4ParticleDefinition& ParticleType)
{

  G4double LowEdgeEnergy , Value;
  G4double FixedEnergy = (LowestKineticEnergy + HighestKineticEnergy)/2. ;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

  if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
  }
  theMeanFreePathTable = new G4PhysicsTable(G4Material::GetNumberOfMaterials());
  PartialSumSigma.clearAndDestroy();
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

void G4MuBremsstrahlung::ComputePartialSumSigma(
                                    const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy,
                                     const G4Material* aMaterial)
// Build the table of cross section per element.The table is built for MATERIALS.
// This table is used by DoIt to select randomly an element in the material. 
{
   G4int Imate = aMaterial->GetIndex();
   G4int NbOfElements = aMaterial->GetNumberOfElements();
   const G4ElementVector* theElementVector = aMaterial->GetElementVector(); 
   const G4double* theAtomNumDensityVector = 
                                    aMaterial->GetAtomicNumDensityVector();
   G4double GammaEnergyCut = (G4Gamma::Gamma()->GetCutsInEnergy())[Imate];

   PartialSumSigma[Imate] = new G4DataVector();

   G4double SIGMA = 0. ;

   for ( G4int Ielem=0 ; Ielem < NbOfElements ; Ielem++ )
      {             
        SIGMA += theAtomNumDensityVector[Ielem] * 
                 ComputeMicroscopicCrossSection( ParticleType, KineticEnergy,
                                            (*theElementVector)(Ielem)->GetZ(), 
                                            (*theElementVector)(Ielem)->GetA(), 
                                                 GammaEnergyCut );
        PartialSumSigma[Imate]->push_back(SIGMA);
   }
}

G4double G4MuBremsstrahlung::ComputeMicroscopicCrossSection(
                                     const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy,
                                           G4double AtomicNumber,
                                           G4double AtomicMass,
                                           G4double GammaEnergyCut)
// Cross section is calculated according to a formula of R.Kokoulin.
{
  G4double TotalEnergy,vcut,vmax,aaa,bbb,hhh,aa,x,ep ;
  G4int kkk ;
  G4double ak1=2.3 ;
  G4int k2=4 ;
  G4double xgi[]={0.03377,0.16940,0.38069,0.61931,0.83060,0.96623};
  G4double wgi[]={0.08566,0.18038,0.23396,0.23396,0.18038,0.08566};
  G4double CrossSection = 0. ;

  TotalEnergy=KineticEnergy+ParticleMass ;
  vcut = GammaEnergyCut/TotalEnergy ;
  vmax = KineticEnergy/TotalEnergy ;
  if(vmax <= vcut) return CrossSection;

  // numerical integration
  aaa=log(vcut) ;
  bbb=log(vmax);
  kkk=int((bbb-aaa)/ak1)+k2 ;
  hhh=(bbb-aaa)/float(kkk) ;

  for(G4int l=0; l<kkk; l++)
  {
    aa=aaa+hhh*float(l) ;
    for(G4int i=0; i<6; i++)
    {
      x=aa+xgi[i]*hhh ;
      ep=exp(x)*TotalEnergy ;
      CrossSection += ep*wgi[i]*ComputeDMicroscopicCrossSection(
                                      ParticleType,KineticEnergy,
                                       AtomicNumber,AtomicMass,ep) ;
    }
  }

  CrossSection *= hhh ; 

  return CrossSection;
}
 
G4double G4MuBremsstrahlung::GetDMicroscopicCrossSection(
                                     const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy,
                                           G4double AtomicNumber,
                                           G4double AtomicMass,   
                                           G4double GammaEnergy)
// get differential cross section
{
   return ComputeDMicroscopicCrossSection(ParticleType,KineticEnergy,
                                          AtomicNumber,AtomicMass,GammaEnergy) ;
} 

G4double G4MuBremsstrahlung::ComputeDMicroscopicCrossSection(
                                     const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy,
                                           G4double AtomicNumber,
                                           G4double AtomicMass,   
                                           G4double GammaEnergy)
//  differential cross section
{
  static const G4double sqrte=sqrt(exp(1.)) ;
  static const G4double bh=202.4,bh1=446.,btf=183.,btf1=1429. ;
  static const G4double rmass=ParticleMass/electron_mass_c2 ;
  static const G4double cc=classic_electr_radius/rmass ;
  static const G4double coeff= 16.*fine_structure_const*cc*cc/3. ;

  G4double dxsection = 0.;

  if( GammaEnergy > KineticEnergy) return dxsection ;

  G4double A = AtomicMass/(g/mole) ;     // !!!!!!!!!!!!!!!!!!!
  G4double E=KineticEnergy+ParticleMass ;
  G4double v=GammaEnergy/E ;
  G4double delta=0.5*ParticleMass*ParticleMass*v/(E-GammaEnergy) ;
  G4double rab0=delta*sqrte ;

  G4double z13=exp(-log(AtomicNumber)/3.) ;
  G4double dn=1.54*exp(0.27*log(A)) ;

  G4double b,b1,dnstar ;

  if(AtomicNumber<1.5)
  {
    b=bh;
    b1=bh1;
    dnstar=dn ;
  }
  else
  {
    b=btf;
    b1=btf1;
    dnstar = exp((1.-1./AtomicNumber)*log(dn)) ;
  }  

  // nucleus contribution logarithm
  G4double rab1=b*z13;
  G4double fn=log(rab1/(dnstar*(electron_mass_c2+rab0*rab1))*
              (ParticleMass+delta*(dnstar*sqrte-2.))) ;
  if(fn <0.) fn = 0. ;
  // electron contribution logarithm
  G4double epmax1=E/(1.+0.5*ParticleMass*rmass/E) ; 
  G4double fe=0.;
  if(GammaEnergy<epmax1)
  {
    G4double rab2=b1*z13*z13 ;
    fe=log(rab2*ParticleMass/((1.+delta*rmass/(electron_mass_c2*sqrte))*
                              (electron_mass_c2+rab0*rab2))) ;
    if(fe<0.) fe=0. ;
  }
  
  dxsection = coeff*(1.-v*(1.-0.75*v))*AtomicNumber*(fn*AtomicNumber+fe)/
                     GammaEnergy ;
  return dxsection ;
}


void G4MuBremsstrahlung::MakeSamplingTables(
                   const G4ParticleDefinition* ParticleType)
{
  G4int nbin;
  G4double AtomicNumber,AtomicWeight,KineticEnergy,
           TotalEnergy,Maxep ;

   ParticleMass = ParticleType->GetPDGMass() ;

   for (G4int iz=0; iz<nzdat; iz++)
   {
     AtomicNumber = zdat[iz];
     AtomicWeight = adat[iz]*g/mole ;

     for (G4int it=0; it<ntdat; it++)
     {
       KineticEnergy = tdat[it];
       TotalEnergy = KineticEnergy + ParticleMass;
       Maxep = KineticEnergy ;

       G4double CrossSection = 0.0 ;

       G4double c,y,ymin,ymax,dy,yy,dx,x,ep,ybin ;

       G4int NbofIntervals ;
       // calculate the differential cross section
       // numerical integration in
       //  log ...............
       c = log(Maxep/CutFixed) ;
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

         ep = CutFixed*exp(c*x) ;

         CrossSection += ep*dx*ComputeDMicroscopicCrossSection(ParticleType,
                                                 KineticEnergy,AtomicNumber,
                                                 AtomicWeight,ep) ;
         if(nbin<NBIN)
         {
           nbin += 1 ;
           ya[nbin]=y ;
           proba[iz][it][nbin] = CrossSection ;
         }
       }

       ya[NBIN] = 0. ;   //   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if(CrossSection > 0.)
       {
         for(G4int ib=0; ib<=nbin; ib++)
         {
           proba[iz][it][ib] /= CrossSection ;
         }
       }
     }
   }
}



G4VParticleChange* G4MuBremsstrahlung::PostStepDoIt(const G4Track& trackData,
                                                  const G4Step& stepData)
{

  static G4double ysmall = -100. ;
  static G4double ytablelow = -5. ;
 
  aParticleChange.Initialize(trackData);
  G4Material* aMaterial=trackData.GetMaterial() ;

  const G4DynamicParticle* aDynamicParticle=trackData.GetDynamicParticle();  

  G4double           KineticEnergy     = aDynamicParticle->GetKineticEnergy();
  G4ParticleMomentum ParticleDirection = 
                                     aDynamicParticle->GetMomentumDirection();

  // Gamma cut in this material
  G4double GammaEnergyCut = 
             (G4Gamma::Gamma()->GetCutsInEnergy())[aMaterial->GetIndex()];

  // check against insufficient energy
  if(KineticEnergy < GammaEnergyCut)
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  // select randomly one element constituing the material  
  G4Element* anElement = SelectRandomAtom(aMaterial);

  G4double TotalEnergy=KineticEnergy+aDynamicParticle->
                             GetDefinition()->GetPDGMass() ;

  G4double dy = 5./G4float(NBIN) ;

  G4double ymin=log(log(GammaEnergyCut/CutFixed)/log(KineticEnergy/CutFixed)) ;

  if(ymin < ysmall)        
    return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  //  sampling using tables 
  G4double v,xc,x,yc,y ;
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
  G4int iymin = G4int((ymin+5.)/dy+0.5) ;  

  if(ymin < ytablelow)
  {
    y = ymin + G4UniformRand()*(ytablelow-ymin) ;
  }
  else
  {
    G4double r = G4UniformRand() ;

    iy = iymin-1 ;
    delmin = proba[izz][itt][NBINminus1]-proba[izz][itt][iymin] ;
    do {
         iy += 1 ;
       } while ((r > (proba[izz][itt][iy]-proba[izz][itt][iymin])/delmin)
                 &&(iy < NBINminus1)) ;

    //sampling is Done uniformly in y in the bin
     y = ya[iy] + G4UniformRand() * ( ya[iy+1] - ya[iy] ) ;
  }

  x = exp(y) ;

  v = CutFixed*exp(x*log(KineticEnergy/CutFixed)) ;            
  if( v <= 0.)
     return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);

  // create G4DynamicParticle object for the Gamma 
  G4double GammaEnergy = v; 

  //  angles of the emitted gamma. ( Z - axis along the parent particle)
  //  Teta = electron_mass_c2/TotalEnergy for the moment .....

  G4double Teta = electron_mass_c2/TotalEnergy ;
  G4double Phi  = twopi * G4UniformRand() ;
  G4double dirx = sin(Teta)*cos(Phi) , diry = sin(Teta)*sin(Phi) ,
           dirz = cos(Teta) ;

  G4ThreeVector GammaDirection ( dirx, diry, dirz);
  GammaDirection.rotateUz(ParticleDirection);   
 
  G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
                                                 GammaDirection, GammaEnergy);

  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(aGamma); 

  // Update the incident particle 
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

  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

G4Element* G4MuBremsstrahlung::SelectRandomAtom(G4Material* aMaterial) const
{
  // select randomly 1 element within the material

  const G4int Index = aMaterial->GetIndex();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  G4double rval = G4UniformRand()*((*PartialSumSigma[Index])[NumberOfElements-1]);
  for ( G4int i=0; i < NumberOfElements; i++ )
    if (rval <= (*PartialSumSigma[Index])[i]) return ((*theElementVector)(i));
  G4cout << " WARNING !!! - The Material " << aMaterial->GetName()
       << " has no elements, NULL pointer returned." << G4endl;
  return NULL;
}

void G4MuBremsstrahlung::PrintInfoDefinition()
{
  G4String comments = "theoretical cross section \n ";
           comments += "         Good description up to 1000 PeV.";

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n    PhysicsTables from " << G4BestUnit(LowerBoundLambda,
                                                     "Energy")
         << " to " << G4BestUnit(UpperBoundLambda,"Energy")
         << " in " << NbinLambda << " bins. \n";
}

