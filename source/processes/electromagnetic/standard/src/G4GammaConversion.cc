// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GammaConversion.cc,v 1.2 1999-03-04 16:31:24 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4GammaConversion physics process --------
//                   by Michel Maire, 24 May 1996
// **************************************************************
// 11-06-96, Added SelectRandomAtom() method, M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 24-06-96, simplification in ComputeMicroscopicCrossSection, M.Maire
// 24-06-96, in DoIt : change the particleType stuff, M.Maire
// 25-06-96, modification in the generation of the teta angle, M.Maire
// 16-09-96, minors optimisations in DoIt. Thanks to P.Urban
//           dynamical array PartialSumSigma
// 13-12-96, fast sampling of epsil below 2 MeV, L.Urban
// 14-01-97, crossection table + meanfreepath table.
//           PartialSumSigma removed, M.Maire
// 14-01-97, in DoIt the positron is always created, even with Ekine=0,
//           for further annihilation, M.Maire
// 14-03-97, new Physics scheme for geant4alpha, M.Maire
// 28-03-97, protection in BuildPhysicsTable, M.Maire
// 19-06-97, correction in ComputeMicroscopicCrossSection, L.Urban
// 04-06-98, in DoIt, secondary production condition: range>min(threshold,safety)
// 13-08-98, new methods SetBining() PrintInfo()
// --------------------------------------------------------------

#include "G4GammaConversion.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// constructor
 
G4GammaConversion::G4GammaConversion(const G4String& processName)
  : G4VDiscreteProcess (processName),            // initialization
    theCrossSectionTable(NULL),
    theMeanFreePathTable(NULL),  
    LowestEnergyLimit (2*electron_mass_c2),
    HighestEnergyLimit(100*GeV),
    NumbBinTable(100)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// destructor
 
G4GammaConversion::~G4GammaConversion()
{
   if (theCrossSectionTable) {
      theCrossSectionTable->clearAndDestroy();
      delete theCrossSectionTable;
   }

   if (theMeanFreePathTable) {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
   }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void G4GammaConversion::SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins)
{
  LowestEnergyLimit = lowE; HighestEnergyLimit = highE; NumbBinTable = nBins;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
void G4GammaConversion::BuildPhysicsTable(const G4ParticleDefinition& GammaType)
// Build microscopic cross section table and mean free path table
{
   G4double LowEdgeEnergy, Value;
   G4PhysicsLogVector* ptrVector;

// Build microscopic cross section tables for the e+e- pair creation

   if (theCrossSectionTable) {
         theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable; }

   theCrossSectionTable = new G4PhysicsTable(G4Element::GetNumberOfElements()) ;
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

// Build mean free path table for the e+e- pair creation

   if (theMeanFreePathTable) {
          theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable; }

   theMeanFreePathTable = new G4PhysicsTable( G4Material::GetNumberOfMaterials() ) ;
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
           }

        theMeanFreePathTable->insertAt( J , ptrVector ) ;

      }
   
   PrintInfoDefinition();        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4GammaConversion::ComputeMicroscopicCrossSection
                              (G4double GammaEnergy, G4double AtomicNumber)
 
// Calculates the microscopic cross section in GEANT4 internal units.
// A parametrized formula from L. Urban is used to estimate the total cross section.
// It gives a good description of the data from 1.5 MeV to 100 GeV.
// below 1.5 MeV: sigma=sigma(1.5MeV)*(GammaEnergy-2electronmass)*(GammaEnergy-2electronmass) 

{
 G4double GammaEnergyLimit = 1.5*MeV;
 G4double CrossSection = 0.0 ;
 if ( AtomicNumber < 1. ) return CrossSection;
 if ( GammaEnergy < 2*electron_mass_c2 ) return CrossSection ;
 
 static const G4double
        a0= 8.7842e+2*microbarn, a1=-1.9625e+3*microbarn, a2= 1.2949e+3*microbarn,
        a3=-2.0028e+2*microbarn, a4= 1.2575e+1*microbarn, a5=-2.8333e-1*microbarn;

 static const G4double
        b0=-1.0342e+1*microbarn, b1= 1.7692e+1*microbarn, b2=-8.2381   *microbarn,
        b3= 1.3063   *microbarn, b4=-9.0815e-2*microbarn, b5= 2.3586e-3*microbarn;

 static const G4double
        c0=-4.5263e+2*microbarn, c1= 1.1161e+3*microbarn, c2=-8.6749e+2*microbarn,
        c3= 2.1773e+2*microbarn, c4=-2.0467e+1*microbarn, c5= 6.5372e-1*microbarn;

 G4double GammaEnergySave = GammaEnergy ;
 if (GammaEnergy < GammaEnergyLimit) GammaEnergy = GammaEnergyLimit ;

 G4double X = log(GammaEnergy/electron_mass_c2), X2=X*X, X3=X2*X, X4=X3*X, X5=X4*X;

 G4double F1 = a0 + a1*X + a2*X2 + a3*X3 + a4*X4 + a5*X5,
          F2 = b0 + b1*X + b2*X2 + b3*X3 + b4*X4 + b5*X5,
          F3 = c0 + c1*X + c2*X2 + c3*X3 + c4*X4 + c5*X5;     

 CrossSection = (AtomicNumber+1.)*
                (F1*AtomicNumber + F2*AtomicNumber*AtomicNumber + F3);

 if (GammaEnergySave < GammaEnergyLimit)
   {
     X=GammaEnergySave-2.*electron_mass_c2;
     CrossSection *= X*X ;
   }

 if (CrossSection < 0.) CrossSection = 0.;

 return CrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
 
G4VParticleChange* G4GammaConversion::PostStepDoIt(const G4Track& aTrack,
                                                  const G4Step&  aStep)
//
// The secondaries e+e- energies are sampled using the Bethe - Heitler cross sections
//     with Coulomb correction.
// A modified version of the random number techniques of Butcher & Messel is used 
//    (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : Effects due to the breakdown of the Born approximation at low energy
//          are ignored.
// Note 2 : The differential cross section implicitly takes account of pair creation
//          in both nuclear and atomic electron fields. However triplet prodution is
//          not generated.
 
{
   aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();

   const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
   G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
   G4ParticleMomentum GammaDirection = aDynamicGamma->GetMomentumDirection();


   G4double epsil ;
   G4double epsil0 = electron_mass_c2/GammaEnergy ;

  // do it fast if GammaEnergy < 2. MeV
   const G4double Egsmall=2.*MeV; 
   if (GammaEnergy<Egsmall) { epsil = epsil0 + (0.5-epsil0)*G4UniformRand(); }

   else
   {  // now comes the case with GammaEnergy >= 2. MeV

   // select randomly one element constituing the material  
   G4Element* anElement = SelectRandomAtom(aDynamicGamma, aMaterial);

   // Extract Coulomb factor for this Element
   G4double FZ = 8.*(anElement->GetIonisation()->GetlogZ3());
   if (GammaEnergy > 50.*MeV) FZ += 8.*(anElement->GetfCoulomb());

   // limits of the screening variable
   G4double screenfac = 136.*epsil0/(anElement->GetIonisation()->GetZ3()) ;
   G4double screenmax = exp ((42.24 - FZ)/8.368) - 0.952 ;
   G4double screenmin = min(4.*screenfac,screenmax) ;

   // limits of the energy sampling
   G4double epsil1 = 0.5 - 0.5*sqrt(1. - screenmin/screenmax) ;
   G4double epsilmin = max(epsil0,epsil1) , epsilrange = 0.5 - epsilmin ;

   //
   // sample the energy rate of the created electron (or positron) 
   //
   //G4double epsil, screenvar, greject ;
   G4double  screenvar, greject ;

   G4double F10 = ScreenFunction1(screenmin) - FZ , F20 = ScreenFunction2(screenmin) - FZ;
   G4double NormF1 = max(F10*epsilrange*epsilrange,0.) ,  NormF2 = max(1.5*F20,0.);

   do {
        if ( NormF1/(NormF1+NormF2) > G4UniformRand() )
             { epsil = 0.5 - epsilrange*pow(G4UniformRand(), 1/3) ;
               screenvar = screenfac/(epsil*(1-epsil));
               greject = (ScreenFunction1(screenvar) - FZ)/F10 ;
             } 
        else { epsil = epsilmin + epsilrange*G4UniformRand();
               screenvar = screenfac/(epsil*(1-epsil));
               greject = (ScreenFunction2(screenvar) - FZ)/F20 ;
             }

   } while( greject < G4UniformRand() );

   }   //  end of epsil sampling.........................

   //
   // fixe charges randomly
   //

   G4double ElectTotEnergy, PositTotEnergy;
   if (RandFlat::shootBit())
     {
       ElectTotEnergy = (1.-epsil)*GammaEnergy;
       PositTotEnergy = epsil*GammaEnergy;
     }
   else
     {
       PositTotEnergy = (1.-epsil)*GammaEnergy;
       ElectTotEnergy = epsil*GammaEnergy;
     }

   //
   // scattered electron (positron) angles. ( Z - axis along the parent photon)
   //
   //  universal distribution suggested by L. Urban (Geant3 manual (1993) Phys211),
   //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

   G4double u;
   const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

   if (9./(9.+d) > G4UniformRand()) u = - log(G4UniformRand()*G4UniformRand())/a1 ;
      else                          u = - log(G4UniformRand()*G4UniformRand())/a2 ;

   G4double Teta = u*electron_mass_c2/GammaEnergy ;
   G4double Phi  = twopi * G4UniformRand() ;
   G4double dirx = sin(Teta)*cos(Phi) , diry = sin(Teta)*sin(Phi) , dirz = cos(Teta) ;
 
   //
   // kinematic of the created pair
   //
   // the electron and positron are assumed to have a symetric angular distribution
   // with respect to the Z axis along the parent photon.

   G4double LocalEnerDeposit = 0. ;
   aParticleChange.SetNumberOfSecondaries(2) ; 

   G4double ElectKineEnergy = max(0.,ElectTotEnergy - electron_mass_c2) ;
   //  condition changed !
     if((G4EnergyLossTables::GetRange(G4Electron::Electron(),
          ElectKineEnergy,aMaterial)>aStep.GetPostStepPoint()->GetSafety())
         ||
        (ElectKineEnergy >
        (G4Electron::Electron()->GetCutsInEnergy())[aMaterial->GetIndex()]))

      {
        G4ThreeVector ElectDirection ( dirx, diry, dirz );
        ElectDirection.rotateUz(GammaDirection);   
 
        // create G4DynamicParticle object for the particle1  
        G4DynamicParticle* aParticle1= new G4DynamicParticle (G4Electron::Electron(),
                                                 ElectDirection, ElectKineEnergy);
        aParticleChange.AddSecondary( aParticle1 ) ; 
       }
    else
       { LocalEnerDeposit += ElectKineEnergy ; }

   // the e+ is always created (even with Ekine=0) for further annihilation.

   G4double PositKineEnergy = max(0.,PositTotEnergy - electron_mass_c2) ;

  //  if (G4EnergyLossTables::GetRange(G4Positron::Positron(),PositKineEnergy,aMaterial)
  //      < min(G4Positron::GetCuts(), aStep.GetPostStepPoint()->GetSafety()) )  
     if((G4EnergyLossTables::GetRange(G4Positron::Positron(),
          PositKineEnergy,aMaterial)<aStep.GetPostStepPoint()->GetSafety())
         &&
        (PositKineEnergy <
        (G4Positron::Positron()->GetCutsInEnergy())[aMaterial->GetIndex()]))

      {
        LocalEnerDeposit += PositKineEnergy ;
        PositKineEnergy = 0. ;
      }
   G4ThreeVector PositDirection ( -dirx, -diry, dirz );
   PositDirection.rotateUz(GammaDirection);   
 
   // create G4DynamicParticle object for the particle2 
   G4DynamicParticle* aParticle2= new G4DynamicParticle (G4Positron::Positron(),
                                                 PositDirection, PositKineEnergy);
   aParticleChange.AddSecondary( aParticle2 ) ; 

   aParticleChange.SetLocalEnergyDeposit( LocalEnerDeposit ) ;

   //
   // Kill the incident photon 
   //

   aParticleChange.SetMomentumChange( 0., 0., 0. ) ;
   aParticleChange.SetEnergyChange( 0. ) ; 
   aParticleChange.SetStatusChange( fStopAndKill ) ;

   //  Reset NbOfInteractionLengthLeft and return aParticleChange
   return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Element* G4GammaConversion::SelectRandomAtom(const G4DynamicParticle* aDynamicGamma,
                                               G4Material* aMaterial)
{
  // select randomly 1 element within the material 

  const G4int NumberOfElements            = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  if (NumberOfElements == 1) return (*theElementVector)(0);

  const G4double* theAtomNumDensityVector = aMaterial->GetAtomicNumDensityVector();

  G4double PartialSumSigma = 0. ;
  G4double rval = G4UniformRand()/MeanFreePath;
 
  for ( G4int i=0 ; i < NumberOfElements ; i++ )
      { PartialSumSigma += theAtomNumDensityVector[i] *
                   GetMicroscopicCrossSection(aDynamicGamma, (*theElementVector)(i));
        if (rval <= PartialSumSigma) return ((*theElementVector)(i));
      }
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements, NULL pointer returned." << endl;
  return NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaConversion::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from a parametrisation. ";
           comments += "Good description from 1.5 MeV to 100 GeV for all Z. \n";
           comments += "        e+e- energies according Bethe-Heitler";
                     
  G4cout << endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestEnergyLimit,"Energy")
         << " to " << G4BestUnit(HighestEnergyLimit,"Energy") 
         << " in " << NumbBinTable << " bins. \n";
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
