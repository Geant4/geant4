// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuNuclearInteraction.cc,v 1.8 2001/02/05 18:02:12 gcosmo Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// $Id: 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------- G4MuNuclearInteraction physics process ---------
//                by Laszlo Urban, May 1998
//      added simple model for hadronic vertex, J.P. Wellisch, November 1998
// --------------------------------------------------------------
// 26/10/98: new corr.s from R.Kokoulin + cleanup , L.Urban
//
#include "G4MuNuclearInteraction.hh"
#include "G4UnitsTable.hh"

// static members ........
G4int G4MuNuclearInteraction::nzdat =  5 ;
G4double G4MuNuclearInteraction::zdat[]={1.,4.,13.,29.,92.};
G4double G4MuNuclearInteraction::adat[]={1.01,9.01,26.98,63.55,238.03};
G4int G4MuNuclearInteraction::ntdat = 8 ;
G4double G4MuNuclearInteraction::tdat[]={1.e3,1.e4,1.e5,1.e6,1.e7,1.e8,
                                                             1.e9,1.e10};
G4int G4MuNuclearInteraction::NBIN = 1000 ;  
G4double G4MuNuclearInteraction::ya[1001]={0.};
G4double G4MuNuclearInteraction::proba[5][8][1001]={0.};
 
G4MuNuclearInteraction::G4MuNuclearInteraction(const G4String& processName)
  : G4VDiscreteProcess(processName),  
    theCrossSectionTable(NULL),
    theMeanFreePathTable(NULL),
    LowestKineticEnergy (1.*GeV),
    HighestKineticEnergy (1000000.*TeV),
    TotBin(100),
    theMuonMinus ( G4MuonMinus::MuonMinus() ),
    theMuonPlus ( G4MuonPlus::MuonPlus() ),
    thePionZero (G4PionZero::PionZero() ),
    GramPerMole(g/mole),
    CutFixed ( 0.200*GeV)
{  }
 
G4MuNuclearInteraction::~G4MuNuclearInteraction()
{
   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }
   if (theCrossSectionTable) {
      theCrossSectionTable->clearAndDestroy();
      delete theCrossSectionTable;
   }

   if (&PartialSumSigma) {
      PartialSumSigma.clearAndDestroy();
   }
}
 
void G4MuNuclearInteraction::SetPhysicsTableBining(G4double lowE,
                                                  G4double highE, G4int nBins)
{
  LowestKineticEnergy = lowE; HighestKineticEnergy = highE ; TotBin = nBins ;
}

void G4MuNuclearInteraction::BuildPhysicsTable(
                                    const G4ParticleDefinition& aParticleType)
{
  G4double LowEdgeEnergy , Value;
  G4PhysicsLogVector* ptrVector;
   
  if (theCrossSectionTable) {
      theCrossSectionTable->clearAndDestroy() ;
      delete theCrossSectionTable ;
  }

  // make tables for the sampling at initialization
  if (theMeanFreePathTable == NULL) MakeSamplingTables(&aParticleType);

  theCrossSectionTable = new G4PhysicsTable (G4Element::GetNumberOfElements()); 
  const G4ElementTable* theElementTable = G4Element::GetElementTable() ;
  G4double AtomicNumber,AtomicWeight ;

  for (G4int J=0; J < G4Element::GetNumberOfElements(); J++ )
  {
    ptrVector = new G4PhysicsLogVector(LowestKineticEnergy,
                                         HighestKineticEnergy,TotBin) ;
    AtomicNumber = (*theElementTable )(J)->GetZ() ;
    AtomicWeight = (*theElementTable )(J)->GetA() ;

    for ( G4int i = 0 ; i < TotBin ; i++)
    {
      LowEdgeEnergy = ptrVector->GetLowEdgeEnergy(i) ;
      Value = ComputeMicroscopicCrossSection(&aParticleType,
                                                 LowEdgeEnergy,
                                                 AtomicNumber,AtomicWeight) ;
      ptrVector->PutValue(i,Value) ;
    }

    theCrossSectionTable->insertAt( J , ptrVector ) ;
  }

  G4double FixedEnergy = (LowestKineticEnergy + HighestKineticEnergy)/2. ;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;
  if (theMeanFreePathTable) {
     theMeanFreePathTable->clearAndDestroy();
     delete theMeanFreePathTable;
  }
  theMeanFreePathTable = new G4PhysicsTable(G4Material::GetNumberOfMaterials());

  PartialSumSigma.clearAndDestroy();
  PartialSumSigma.resize(G4Material::GetNumberOfMaterials());


  for (G4int K=0 ; K < G4Material::GetNumberOfMaterials(); K++ )  
  { 
    ptrVector = new G4PhysicsLogVector(LowestKineticEnergy,
                                         HighestKineticEnergy,
                                         TotBin ) ;

    const G4Material* material= (*theMaterialTable)[K];

    for ( G4int i = 0 ; i < TotBin ; i++ )      
    {
      LowEdgeEnergy = ptrVector->GetLowEdgeEnergy( i ) ;
      Value = ComputeMeanFreePath( &aParticleType, LowEdgeEnergy,
                                         material );  
        
      ptrVector->PutValue( i , Value ) ;
    }

    theMeanFreePathTable->insertAt( K , ptrVector );

    // Compute the PartialSumSigma table at a given fixed energy
    ComputePartialSumSigma( &aParticleType, FixedEnergy, material);       
  }

  if (&aParticleType == theMuonPlus) PrintInfoDefinition();
}

void G4MuNuclearInteraction::ComputePartialSumSigma(
                                   const G4ParticleDefinition* ParticleType,
                                               G4double KineticEnergy,
                                               const G4Material* aMaterial)

// Build the table of cross section per element. The table is built for MATERIALS.
// This table is used by DoIt to select randomly an element in the material. 
{
   G4int Imate = aMaterial->GetIndex();
   G4int NbOfElements = aMaterial->GetNumberOfElements();
   const G4ElementVector* theElementVector = aMaterial->GetElementVector(); 
   const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();

   PartialSumSigma[Imate] = new G4DataVector();

   G4double SIGMA = 0. ;

   for ( G4int Ielem=0 ; Ielem < NbOfElements ; Ielem++ )
      {             
        SIGMA += theAtomNumDensityVector[Ielem] * 
               ComputeMicroscopicCrossSection( ParticleType, KineticEnergy,
                                        (*theElementVector)(Ielem)->GetZ(),
                                        (*theElementVector)(Ielem)->GetA()) ;
        PartialSumSigma[Imate]->push_back(SIGMA);
      }
}

G4double G4MuNuclearInteraction::ComputeMicroscopicCrossSection(
                                const G4ParticleDefinition* ParticleType,
                                G4double KineticEnergy,
                                G4double AtomicNumber,G4double AtomicWeight)
{
  static const G4double
  xgi[] ={ 0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801 };
  static const G4double
  wgi[] ={ 0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506 };
  static const G4double ak1=6.9 ;
  static const G4double ak2=1.0 ;

  G4double Mass,epmin,epmax,epln,ep,aaa,bbb,hhh,x ;
  G4int kkk ;

  Mass = ParticleType->GetPDGMass() ;

  G4double CrossSection = 0.0 ;

  if ( AtomicNumber < 1. ) return CrossSection;
  if ( KineticEnergy <= CutFixed  ) return CrossSection; 

  epmin = CutFixed ;
  epmax = KineticEnergy + Mass - 0.5*proton_mass_c2 ;

  aaa = log(epmin) ;
  bbb = log(epmax) ;
  kkk = int((bbb-aaa)/ak1)+ak2 ;
  hhh = (bbb-aaa)/kkk ;

  for (G4int l=0 ; l<kkk; l++)
  {
    x = aaa + hhh*l ;
    for (G4int ll=0; ll<8; ll++)
    {
      epln=x+xgi[ll]*hhh ;
      ep = exp(epln) ;
      CrossSection += ep*wgi[ll]* ComputeDMicroscopicCrossSection(ParticleType,
                                                KineticEnergy,  
                                                AtomicNumber,AtomicWeight,     
                                                ep) ;              
    }
  }
  CrossSection *= hhh ;

  if (CrossSection < 0.) CrossSection = 0.;

  return CrossSection;
}

G4double G4MuNuclearInteraction::ComputeDMicroscopicCrossSection(
                                 const G4ParticleDefinition* ParticleType,
                                 G4double KineticEnergy,
                                 G4double AtomicNumber,G4double AtomicWeight,
                                 G4double epsilon)
 // Calculates the differential (D) microscopic cross section 
 //   using the cross section formula of R.P. Kokoulin (18/01/98)
{
  const G4double alam2 = 0.400*GeV*GeV ;
  const G4double alam  = 0.632456*GeV ;
  const G4double coeffn = fine_structure_const/pi ;   

  G4double ep,a,aeff,sigph,v,v1,v2,mass2,up,down ;
  G4double ParticleMass = ParticleType->GetPDGMass() ;
  G4double TotalEnergy = KineticEnergy + ParticleMass ;

  G4double DCrossSection = 0. ;

  if((epsilon >= TotalEnergy - 0.5*proton_mass_c2) 
  ||
     (epsilon <= CutFixed))
       return DCrossSection ;

  ep = epsilon/GeV ;
  a = AtomicWeight/(g/mole) ;

  aeff = 0.22*a+0.78*exp(0.89*log(a)) ;                   //shadowing 
  sigph = (49.2+11.1*log(ep)+151.8/sqrt(ep))*microbarn ; //!!!!!!!!!!! 
  
  v=epsilon/TotalEnergy ;
  v1 = 1.-v ;
  v2 = v*v ;
  mass2 = ParticleMass*ParticleMass ;

  up = TotalEnergy*TotalEnergy*v1/mass2*(1.+mass2*v2/(alam2*v1)) ;
  down = 1.+epsilon/alam*(1.+alam/(2.*proton_mass_c2)+epsilon/alam) ;

  DCrossSection = coeffn*aeff*sigph/epsilon*
                  (-v1+(v1+0.5*v2*(1.+2.*mass2/alam2))*log(up/down)) ;

  if( DCrossSection < 0.) 
      DCrossSection = 0. ; 

  return DCrossSection ;
}

void G4MuNuclearInteraction::MakeSamplingTables(
                                     const G4ParticleDefinition* ParticleType)
{
 
  G4int nbin;
  G4double AtomicNumber,AtomicWeight,KineticEnergy,
           TotalEnergy,Maxep ;  

  G4double ParticleMass = ParticleType->GetPDGMass() ;

  for (G4int iz=0; iz<nzdat; iz++)
  {
    AtomicNumber = zdat[iz];
    AtomicWeight = adat[iz]*GramPerMole ;  

    for (G4int it=0; it<ntdat; it++)
    {
      KineticEnergy = tdat[it];
      TotalEnergy = KineticEnergy + ParticleMass;
      Maxep = TotalEnergy - 0.5*proton_mass_c2 ;

      G4double CrossSection = 0.0 ;

      G4double c,y,ymin,ymax,dy,yy,dx,x,ep ;

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
      ya[NBIN]=0. ;

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

G4VParticleChange* G4MuNuclearInteraction::PostStepDoIt(
                                                  const G4Track& trackData,
                                                  const G4Step& stepData)
{
   static const G4double Mass=theMuonPlus->GetPDGMass() ;
   static const G4double m0=0.2*GeV ;

   aParticleChange.Initialize(trackData);
   G4Material* aMaterial=trackData.GetMaterial() ;

   const G4DynamicParticle* aDynamicParticle=trackData.GetDynamicParticle();  

   G4double           KineticEnergy     = aDynamicParticle->GetKineticEnergy();
   G4ParticleMomentum ParticleDirection = 
                                      aDynamicParticle->GetMomentumDirection();

   // limits of the energy sampling
   G4double TotalEnergy = KineticEnergy + Mass ;
   G4double epmin = CutFixed ;
   G4double epmax = TotalEnergy - 0.5*proton_mass_c2 ;
   // check against insufficient energy
   if (epmax <= epmin )   
   {
     aParticleChange.SetMomentumChange( ParticleDirection );
     aParticleChange.SetEnergyChange( KineticEnergy );
     aParticleChange.SetLocalEnergyDeposit (0.); 
     aParticleChange.SetNumberOfSecondaries(0);
     return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
   }

   // select randomly one element constituing the material  
   G4Element* anElement = SelectRandomAtom(aMaterial);

   // sample  energy of the secondary ( pi0)
   //  sampling using tables 
   G4double ep,xc,x,yc,y ;
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

  //sample energy transfer according to the sampling table

   G4double r = G4UniformRand() ;
 
   iy = -1 ;
   do {
        iy += 1 ;
      } while (((proba[izz][itt][iy]) < r)&&(iy < NBINminus1)) ;

   //sampling is Done uniformly in y in the bin
   if( iy < NBIN )
     y = ya[iy] + G4UniformRand() * ( ya[iy+1] - ya[iy] ) ;
   else
     y = ya[iy] ;

   x = exp(y) ;
   ep = epmin*exp(x*log(epmax/epmin)) ;                              

   // sample scattering angle of mu, but first t should be sampled.
   G4double yy = ep/TotalEnergy ;
   G4double tmin=Mass*Mass*yy*yy/(1.-yy) ;
   G4double tmax=2.*proton_mass_c2*ep ;
   G4double t1,t2,t,w1,w2,w3,y1,y2,y3,rej ;
   if(m0<ep)
   {
     t1=m0*m0;
     t2=ep*ep;
   }
   else
   {
     t1=ep*ep;
     t2=m0*m0;
   }

   w1=tmax*t1;
   w2=tmax+t1 ;
   w3=tmax*(tmin+t1)/(tmin*w2);
   y1=1.-yy;
   y2=0.5*yy*yy;
   y3=y1+y2;

   // now the sampling of t
   G4int ntry=0;
   do
   {
     ntry += 1 ;
     t=w1/(w2*exp(G4UniformRand()*log(w3))-tmax) ;
     rej = (1.-t/tmax)*(y1*(1.-tmin/t)+y2)/(y3*(1.-t/t2)); 
   } while (G4UniformRand() > rej) ;

   // compute angle from t
   G4double sinth2,theta ; //  sinth2 = sin(theta/2)*sin(theta/2) !
   sinth2 = 0.5*(t-tmin)/(2.*(TotalEnergy*(TotalEnergy-ep)-Mass*Mass)-tmin) ;
   theta = acos(1.-2.*sinth2) ;
   
   G4double phi=twopi*G4UniformRand() ;
   G4double sinth=sin(theta) ;
   G4double dirx=sinth*cos(phi) , diry=sinth*sin(phi) , dirz=cos(theta);

   G4ThreeVector finalDirection(dirx,diry,dirz) ;
   finalDirection.rotateUz(ParticleDirection) ;

   G4double NewKinEnergy = KineticEnergy - ep ;
   G4double finalMomentum=sqrt(NewKinEnergy*
                       (NewKinEnergy+2.*Mass)) ;

   G4double Ef=NewKinEnergy+Mass ;
   G4double initMomentum=sqrt(KineticEnergy*(TotalEnergy+Mass)) ;

   G4double Q2=2.*(TotalEnergy*Ef-initMomentum*finalMomentum*cos(theta)-Mass*Mass) ;

   aParticleChange.SetMomentumChange( finalDirection );
   aParticleChange.SetEnergyChange( NewKinEnergy );

   G4LorentzVector primaryMomentum(initMomentum*ParticleDirection, TotalEnergy);
   G4LorentzVector fsMomentum(finalMomentum*finalDirection, Ef);
   G4LorentzVector momentumTransfere = primaryMomentum-fsMomentum;

   G4DynamicParticle* aGamma = 
           new G4DynamicParticle(G4Gamma::Gamma(), momentumTransfere);
   G4Track gammaTrack(aGamma, trackData.GetGlobalTime(), trackData.GetPosition() );
   gammaTrack.SetStep(trackData.GetStep());
   G4Nucleus theTarget(aMaterial);

   G4VParticleChange* aHadronicFS;
   aHadronicFS = theHadronicVertex.ApplyYourself(theTarget, gammaTrack);
   delete aGamma;

   G4int numSecondaries = aHadronicFS->GetNumberOfSecondaries();
   aParticleChange.SetNumberOfSecondaries(numSecondaries);

   G4ParticleMomentum secondaryMomentum = G4ThreeVector(0.,0.,0.);
   for(G4int iSec=0; iSec<numSecondaries; iSec++) 
   {
     secondaryMomentum 
       = secondaryMomentum + aHadronicFS->GetSecondary(iSec)->GetMomentum();
     aParticleChange.AddSecondary(aHadronicFS->GetSecondary(iSec));
   }
   aHadronicFS->Clear();

   return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
}

G4Element* G4MuNuclearInteraction::SelectRandomAtom(G4Material* aMaterial) const
{
  // select randomly 1 element within the material

  const G4int Index = aMaterial->GetIndex();
  const G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  G4double rval = G4UniformRand()*((*PartialSumSigma[Index])[NumberOfElements-1]);
  for ( G4int i=0; i < NumberOfElements; i++ )
    if (rval <= (*PartialSumSigma[Index])[i]) return ((*theElementVector)(i));
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements, NULL pointer returned." << G4endl;
  return NULL;
}
void G4MuNuclearInteraction::PrintInfoDefinition()
{
  G4String comments = "cross sections from R. Kokoulin \n ";
           comments += "         Good description up to 1000 PeV.";

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n    PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
                                                     "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins. \n";

  G4cout << G4endl;
}

