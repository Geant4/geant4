//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4MuNuclearInteraction.cc,v 1.14 2010-09-08 08:59:29 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------- G4MuNuclearInteraction physics process ---------
//                by Laszlo Urban, May 1998
//      added simple model for hadronic vertex, J.P. Wellisch, November 1998
// --------------------------------------------------------------
// 26/10/1998: new corr.s from R.Kokoulin + cleanup , L.Urban
// 23/01/2009  V.Ivanchenko Add deregistration
//

#include "G4MuNuclearInteraction.hh"
#include "G4UnitsTable.hh"
#include "G4HadronicProcessStore.hh"

// static members ........
G4int G4MuNuclearInteraction::nzdat =  5 ;
G4double G4MuNuclearInteraction::zdat[]={1.,4.,13.,29.,92.};
G4double G4MuNuclearInteraction::adat[]={1.01,9.01,26.98,63.55,238.03};
G4int G4MuNuclearInteraction::ntdat = 8 ;
G4double G4MuNuclearInteraction::tdat[]={1.e3,1.e4,1.e5,1.e6,1.e7,1.e8,
                                                             1.e9,1.e10};
G4int G4MuNuclearInteraction::NBIN = 1000 ;  
G4double G4MuNuclearInteraction::ya[1001]={0.};
G4double G4MuNuclearInteraction::proba[5][8][1001]={{{0.}}};
 
G4MuNuclearInteraction::G4MuNuclearInteraction(const G4String& processName)
  : G4VDiscreteProcess(processName),  
    theMeanFreePathTable(0),
    theCrossSectionTable(0),
    LowestKineticEnergy (1.*GeV),
    HighestKineticEnergy (1000000.*TeV),
    TotBin(100),
    CutFixed ( 0.200*GeV),
    GramPerMole(g/mole),
    theMuonMinus ( G4MuonMinus::MuonMinus() ),
    theMuonPlus ( G4MuonPlus::MuonPlus() ),
    thePionZero (G4PionZero::PionZero() )
{  
  SetProcessSubType(fHadronInelastic);
  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
}
 
G4MuNuclearInteraction::~G4MuNuclearInteraction()
{
  G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);

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


G4bool G4MuNuclearInteraction::IsApplicable(const G4ParticleDefinition& particle)
{
   return(   (&particle == theMuonMinus)||(&particle == theMuonPlus)) ;
}

void G4MuNuclearInteraction::PreparePhysicsTable(
                                    const G4ParticleDefinition& aParticleType)
{
  G4HadronicProcessStore::Instance()
    ->RegisterParticleForExtraProcess(this, &aParticleType);
}

void G4MuNuclearInteraction::BuildPhysicsTable(
                                    const G4ParticleDefinition& aParticleType)
{
  G4HadronicProcessStore::Instance()->PrintInfo(&aParticleType);

  G4double LowEdgeEnergy , Value;
  G4PhysicsLogVector* ptrVector;
   
  if (theCrossSectionTable) {
    theCrossSectionTable->clearAndDestroy() ;
    delete theCrossSectionTable ;
  }

  // make tables for the sampling at initialization
  if (theMeanFreePathTable == 0) MakeSamplingTables(&aParticleType);

  theCrossSectionTable = new G4PhysicsTable (G4Element::GetNumberOfElements()); 
  const G4ElementTable* theElementTable = G4Element::GetElementTable() ;
  G4double AtomicNumber,AtomicWeight ;

  for (size_t J=0; J < G4Element::GetNumberOfElements(); J++ )
  {
    ptrVector = new G4PhysicsLogVector(LowestKineticEnergy,
				       HighestKineticEnergy,TotBin) ;
    AtomicNumber = (*theElementTable )[J]->GetZ() ;
    AtomicWeight = (*theElementTable )[J]->GetA() ;

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


  for (size_t K=0 ; K < G4Material::GetNumberOfMaterials(); K++ )  
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
					 (*theElementVector)[Ielem]->GetZ(),
					 (*theElementVector)[Ielem]->GetA()) ;
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
  if (epmax <= epmin ) return CrossSection; // NaN bug correction

  aaa = std::log(epmin) ;
  bbb = std::log(epmax) ;
  kkk = G4int((bbb-aaa)/ak1 +ak2) ;
  hhh = (bbb-aaa)/kkk ;

  for (G4int l=0 ; l<kkk; l++)
  {
    x = aaa + hhh*l ;
    for (G4int ll=0; ll<8; ll++)
    {
      epln=x+xgi[ll]*hhh ;
      ep = std::exp(epln) ;
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
                                 G4double, G4double AtomicWeight,
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

  aeff = 0.22*a+0.78*std::exp(0.89*std::log(a)) ;                   //shadowing 
  sigph = (49.2+11.1*std::log(ep)+151.8/std::sqrt(ep))*microbarn ; //!!!!!!!!!!! 
  
  v=epsilon/TotalEnergy ;
  v1 = 1.-v ;
  v2 = v*v ;
  mass2 = ParticleMass*ParticleMass ;

  up = TotalEnergy*TotalEnergy*v1/mass2*(1.+mass2*v2/(alam2*v1)) ;
  down = 1.+epsilon/alam*(1.+alam/(2.*proton_mass_c2)+epsilon/alam) ;

  DCrossSection = coeffn*aeff*sigph/epsilon*
                  (-v1+(v1+0.5*v2*(1.+2.*mass2/alam2))*std::log(up/down)) ;

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

      // calculate the differential cross section
      // numerical integration in    
      //  log ...............
      c = std::log(Maxep/CutFixed) ;
      ymin = -5. ;
      ymax = 0. ;
      dy = (ymax-ymin)/NBIN ; 

      nbin=-1;              

      y = ymin - 0.5*dy ;
      yy = ymin - dy ;
      for (G4int i=0 ; i<NBIN; i++)
      {
        y += dy ;
        x = std::exp(y) ;
        yy += dy ;
        dx = std::exp(yy+dy)-std::exp(yy) ;
      
        ep = CutFixed*std::exp(c*x) ;

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
     aParticleChange.ProposeMomentumDirection( ParticleDirection );
     aParticleChange.ProposeEnergy( KineticEnergy );
     aParticleChange.ProposeLocalEnergyDeposit (0.); 
     aParticleChange.SetNumberOfSecondaries(0);
     return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
   }

   // select randomly one element constituing the material  
   G4Element* anElement = SelectRandomAtom(aMaterial);

   // sample  energy of the secondary ( pi0)
   //  sampling using tables 
   G4double ep,x,y ;
   G4int iy ;

   // select sampling table ;
   G4double lnZ = std::log(anElement->GetZ()) ;
   G4double delmin = 1.e10 ;
   G4double del ;
   G4int izz=0,itt=0,NBINminus1 ;
   NBINminus1 = NBIN-1 ;
   for (G4int iz=0; iz<nzdat; iz++)
   {
     del = std::abs(lnZ-std::log(zdat[iz])) ;
     if(del<delmin)
     {
        delmin=del ;
        izz=iz ;
     }
   }
   delmin = 1.e10 ;
   for (G4int it=0; it<ntdat; it++)
   {
     del = std::abs(std::log(KineticEnergy)-std::log(tdat[it])) ;
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

   x = std::exp(y) ;
   ep = epmin*std::exp(x*std::log(epmax/epmin)) ;                              

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
     t=w1/(w2*std::exp(G4UniformRand()*std::log(w3))-tmax) ;
     rej = (1.-t/tmax)*(y1*(1.-tmin/t)+y2)/(y3*(1.-t/t2)); 
   } while (G4UniformRand() > rej) ;

   // compute angle from t
   G4double sinth2,theta ; //  sinth2 = std::sin(theta/2)*std::sin(theta/2) !
   sinth2 = 0.5*(t-tmin)/(2.*(TotalEnergy*(TotalEnergy-ep)-Mass*Mass)-tmin) ;
   theta = std::acos(1.-2.*sinth2) ;
   
   G4double phi=twopi*G4UniformRand() ;
   G4double sinth=std::sin(theta) ;
   G4double dirx=sinth*std::cos(phi) , diry=sinth*std::sin(phi) , dirz=std::cos(theta);

   G4ThreeVector finalDirection(dirx,diry,dirz) ;
   finalDirection.rotateUz(ParticleDirection) ;

   G4double NewKinEnergy = KineticEnergy - ep ;
   G4double finalMomentum=std::sqrt(NewKinEnergy*
                       (NewKinEnergy+2.*Mass)) ;

   G4double Ef=NewKinEnergy+Mass ;
   G4double initMomentum=std::sqrt(KineticEnergy*(TotalEnergy+Mass)) ;

   // G4double Q2=2.*(TotalEnergy*Ef-initMomentum*finalMomentum*std::cos(theta)-Mass*Mass) ;

   aParticleChange.ProposeMomentumDirection( finalDirection );
   aParticleChange.ProposeEnergy( NewKinEnergy );

   G4LorentzVector primaryMomentum(initMomentum*ParticleDirection, TotalEnergy);
   G4LorentzVector fsMomentum(finalMomentum*finalDirection, Ef);
   G4LorentzVector momentumTransfere = primaryMomentum-fsMomentum;

   G4DynamicParticle* aGamma = 
           new G4DynamicParticle(G4Gamma::Gamma(), momentumTransfere);
   G4Track gammaTrack(aGamma, trackData.GetGlobalTime(), trackData.GetPosition() );
   gammaTrack.SetStep(trackData.GetStep());
   G4Nucleus theTarget(aMaterial);

   G4VParticleChange* aHadronicFS = theHadronicVertex.ApplyYourself(theTarget, gammaTrack);
   //delete aGamma;

   G4int numSecondaries = aHadronicFS->GetNumberOfSecondaries();
   aParticleChange.SetNumberOfSecondaries(numSecondaries);

   //   G4ParticleMomentum secondaryMomentum = G4ThreeVector(0.,0.,0.);
   for(G4int iSec=0; iSec<numSecondaries; iSec++) 
   {
     //secondaryMomentum 
     //  = secondaryMomentum + aHadronicFS->GetSecondary(iSec)->GetMomentum();
     aParticleChange.AddSecondary(aHadronicFS->GetSecondary(iSec));
   }
   aHadronicFS->Clear();

   return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
}

G4Element* G4MuNuclearInteraction::SelectRandomAtom(G4Material* aMaterial) const
{
  // select randomly 1 element within the material

  G4int Index = aMaterial->GetIndex();
  G4int NumberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  if(1 == NumberOfElements) return ((*theElementVector)[0]);

  G4double rval = G4UniformRand()*((*PartialSumSigma[Index])[NumberOfElements-1]);
  for ( G4int i=0; i < NumberOfElements; i++ ) {
    if (rval <= (*PartialSumSigma[Index])[i]) return ((*theElementVector)[i]);
  }
  G4cout << "G4MuNuclearInteraction WARNING !!! no element selected for '"
	 << aMaterial->GetName()
	 << " 1st element returned." << G4endl;
  return ((*theElementVector)[0]);
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

G4double G4MuNuclearInteraction::GetMeanFreePath(const G4Track& trackData,
						 G4double,
						 G4ForceCondition* condition)
{
   const G4DynamicParticle* aDynamicParticle;
   G4Material* aMaterial;
   G4double MeanFreePath;
   G4bool isOutRange ;
 
   *condition = NotForced ;

   aDynamicParticle = trackData.GetDynamicParticle();
   aMaterial = trackData.GetMaterial();

   G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();

   if (KineticEnergy <  LowestKineticEnergy)
     MeanFreePath = DBL_MAX ;
   else {
     if (KineticEnergy > HighestKineticEnergy) 
                                 KineticEnergy = 0.99*HighestKineticEnergy ;
     MeanFreePath = (*theMeanFreePathTable)(aMaterial->GetIndex())->
                                 GetValue( KineticEnergy, isOutRange );
   }
   return MeanFreePath; 
} 

G4double G4MuNuclearInteraction::ComputeMeanFreePath(
                                     const G4ParticleDefinition* ParticleType,
				     G4double KineticEnergy,
				     const G4Material* aMaterial)
{
  const G4ElementVector* theElementVector = aMaterial->GetElementVector() ;
  const G4double* theAtomNumDensityVector = 
                                      aMaterial->GetAtomicNumDensityVector();

  G4double SIGMA = 0 ;

  for ( size_t i=0 ; i < aMaterial->GetNumberOfElements() ; i++ )
  {             
    SIGMA += theAtomNumDensityVector[i] * 
             ComputeMicroscopicCrossSection( ParticleType, KineticEnergy,
                                          (*theElementVector)[i]->GetZ(),
                                          (*theElementVector)[i]->GetA()) ;
  }       

  return SIGMA<=0.0 ? DBL_MAX : 1./SIGMA ;
}


