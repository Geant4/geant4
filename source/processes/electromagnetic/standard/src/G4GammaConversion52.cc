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
//
// $Id: G4GammaConversion52.cc,v 1.4 2006/10/16 15:26:49 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//------------------ G4GammaConversion52 physics process -------------------------
//                   by Michel Maire, 24 May 1996
//
// 11-06-96 Added SelectRandomAtom() method, M.Maire
// 21-06-96 SetCuts implementation, M.Maire
// 24-06-96 simplification in ComputeCrossSectionPerAtom, M.Maire
// 24-06-96 in DoIt : change the particleType stuff, M.Maire
// 25-06-96 modification in the generation of the teta angle, M.Maire
// 16-09-96 minors optimisations in DoIt. Thanks to P.Urban
//          dynamical array PartialSumSigma
// 13-12-96 fast sampling of epsil below 2 MeV, L.Urban
// 14-01-97 crossection table + meanfreepath table.
//          PartialSumSigma removed, M.Maire
// 14-01-97 in DoIt the positron is always created, even with Ekine=0,
//          for further annihilation, M.Maire
// 14-03-97 new Physics scheme for geant4alpha, M.Maire
// 28-03-97 protection in BuildPhysicsTable, M.Maire
// 19-06-97 correction in ComputeCrossSectionPerAtom, L.Urban
// 04-06-98 in DoIt, secondary production condition:
//            range>std::min(threshold,safety)
// 13-08-98 new methods SetBining() PrintInfo()
// 28-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 11-07-01 PostStepDoIt - sampling epsil: power(rndm,0.333333)
// 13-07-01 DoIt: suppression of production cut for the (e-,e+) (mma)
// 06-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01 BuildThePhysicsTable() called from constructor (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 20-09-01 DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01 come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 11-01-02 ComputeCrossSection: correction of extrapolation below EnergyLimit
// 21-03-02 DoIt: correction of the e+e- angular distribution (bug 363) mma
// 08-11-04 Remove of Store/Retrieve tables (V.Ivantchenko)
// 04-05-05 Add 52 to class name (V.Ivanchenko)
// 16-11-05 replace shootBit() by G4UniformRand()  mma
// -----------------------------------------------------------------------------

#include "G4GammaConversion52.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
using namespace std;
 
G4GammaConversion52::G4GammaConversion52(const G4String& processName,
    G4ProcessType type):G4VDiscreteProcess (processName, type),
    theCrossSectionTable(NULL),
    theMeanFreePathTable(NULL),  
    LowestEnergyLimit (2*electron_mass_c2),
    HighestEnergyLimit(100*GeV),
    NumbBinTable(100),
    fminimalEnergy(1*eV)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
// destructor
 
G4GammaConversion52::~G4GammaConversion52()
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4GammaConversion52::IsApplicable( const G4ParticleDefinition& particle)
{
   return ( &particle == G4Gamma::Gamma() ); 
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void G4GammaConversion52::SetPhysicsTableBining(
                                  G4double lowE, G4double highE, G4int nBins)
{
  LowestEnergyLimit = lowE; HighestEnergyLimit = highE; NumbBinTable = nBins;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4GammaConversion52::BuildPhysicsTable(const G4ParticleDefinition&)
// Build cross section and mean free path tables
{
   G4double LowEdgeEnergy, Value;
   G4PhysicsLogVector* ptrVector;

// Build cross section per atom tables for the e+e- pair creation

   if (theCrossSectionTable) {
         theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable;}

   theCrossSectionTable = new G4PhysicsTable(G4Element::GetNumberOfElements());
   const G4ElementTable* theElementTable = G4Element::GetElementTable();
   G4double AtomicNumber;
   size_t J;

   for ( J=0 ; J < G4Element::GetNumberOfElements(); J++ )  
       { 
        //create physics vector then fill it ....
        ptrVector = new G4PhysicsLogVector(LowestEnergyLimit,HighestEnergyLimit,
                                           NumbBinTable );
        AtomicNumber = (*theElementTable)[J]->GetZ();
 
        for ( G4int i = 0 ; i < NumbBinTable ; i++ )      
           {
             LowEdgeEnergy = ptrVector->GetLowEdgeEnergy( i ) ;
             Value = ComputeCrossSectionPerAtom( LowEdgeEnergy, AtomicNumber);  
             ptrVector->PutValue( i , Value ) ;
           }

        theCrossSectionTable->insertAt( J , ptrVector ) ;

      }

// Build mean free path table for the e+e- pair creation

   if (theMeanFreePathTable) 
     { theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable;}

   theMeanFreePathTable= new G4PhysicsTable(G4Material::GetNumberOfMaterials());
   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
   G4Material* material;

   for ( J=0 ; J < G4Material::GetNumberOfMaterials(); J++ )
       {
        //create physics vector then fill it ....
        ptrVector = new G4PhysicsLogVector(LowestEnergyLimit,HighestEnergyLimit,
                                           NumbBinTable);
        material = (*theMaterialTable)[J];

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GammaConversion52::ComputeCrossSectionPerAtom
                              (G4double GammaEnergy, G4double AtomicNumber)

// Calculates the microscopic cross section in GEANT4 internal units.
// A parametrized formula from L. Urban is used to estimate
// the total cross section.
// It gives a good description of the data from 1.5 MeV to 100 GeV.
// below 1.5 MeV: sigma=sigma(1.5MeV)*(GammaEnergy-2electronmass)
//                                   *(GammaEnergy-2electronmass) 

{
 G4double GammaEnergyLimit = 1.5*MeV;
 G4double CrossSection = 0.0 ;
 if ( AtomicNumber < 1. ) return CrossSection;
 if ( GammaEnergy < 2*electron_mass_c2 ) return CrossSection;

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

 G4double X=log(GammaEnergy/electron_mass_c2),X2=X*X, X3=X2*X, X4=X3*X, X5=X4*X;

 G4double F1 = a0 + a1*X + a2*X2 + a3*X3 + a4*X4 + a5*X5,
          F2 = b0 + b1*X + b2*X2 + b3*X3 + b4*X4 + b5*X5,
          F3 = c0 + c1*X + c2*X2 + c3*X3 + c4*X4 + c5*X5;     

 CrossSection = (AtomicNumber+1.)*
                (F1*AtomicNumber + F2*AtomicNumber*AtomicNumber + F3);

 if (GammaEnergySave < GammaEnergyLimit)
   {
     X = (GammaEnergySave - 2.*electron_mass_c2)
        /(GammaEnergyLimit- 2.*electron_mass_c2);
     CrossSection *= X*X;
   }

 if (CrossSection < 0.) CrossSection = 0.;

 return CrossSection;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GammaConversion52::ComputeMeanFreePath(G4double GammaEnergy,
                                                G4Material* aMaterial)

// computes and returns the photon mean free path in GEANT4 internal units

{
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();   

  G4double SIGMA = 0 ;

  for ( size_t i=0 ; i < aMaterial->GetNumberOfElements() ; i++ )
      {             
            SIGMA += NbOfAtomsPerVolume[i] * 
                     ComputeCrossSectionPerAtom(GammaEnergy,
                                               (*theElementVector)[i]->GetZ());
      }       

  return SIGMA > DBL_MIN ? 1./SIGMA : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GammaConversion52::GetCrossSectionPerAtom(
                                   const G4DynamicParticle* aDynamicGamma,
                                         G4Element*         anElement)
 
// gives the total cross section per atom in GEANT4 internal units

{
   G4double crossSection;
   G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
   G4bool isOutRange ;

   if (GammaEnergy <  LowestEnergyLimit)
     crossSection = 0. ;
   else {
     if (GammaEnergy > HighestEnergyLimit) GammaEnergy=0.99*HighestEnergyLimit;
     crossSection = (*theCrossSectionTable)(anElement->GetIndex())->
                    GetValue( GammaEnergy, isOutRange );
   }

   return crossSection; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GammaConversion52::GetMeanFreePath(const G4Track& aTrack,
                                                     G4double,
                                                     G4ForceCondition*)

// returns the photon mean free path in GEANT4 internal units
// (MeanFreePath is a private member of the class)

{
   const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
   G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
   G4Material* aMaterial = aTrack.GetMaterial();

   G4bool isOutRange;

   if (GammaEnergy <  LowestEnergyLimit)
     MeanFreePath = DBL_MAX;
   else {
     if (GammaEnergy > HighestEnergyLimit) GammaEnergy=0.99*HighestEnergyLimit;
     MeanFreePath = (*theMeanFreePathTable)(aMaterial->GetIndex())->
                    GetValue( GammaEnergy, isOutRange );
   }

   return MeanFreePath; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VParticleChange* G4GammaConversion52::PostStepDoIt(const G4Track& aTrack,
                                                  const G4Step&  aStep)
//
// The secondaries e+e- energies are sampled using the Bethe - Heitler
// cross sections with Coulomb correction.
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : Effects due to the breakdown of the Born approximation at
//          low energy are ignored.
// Note 2 : The differential cross section implicitly takes account of 
//          pair creation in both nuclear and atomic electron fields.
//          However triplet prodution is not generated.

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
   G4double screenfac = 136.*epsil0/(anElement->GetIonisation()->GetZ3());
   G4double screenmax = exp ((42.24 - FZ)/8.368) - 0.952 ;
   G4double screenmin = min(4.*screenfac,screenmax);

   // limits of the energy sampling
   G4double epsil1 = 0.5 - 0.5*sqrt(1. - screenmin/screenmax) ;
   G4double epsilmin = max(epsil0,epsil1) , epsilrange = 0.5 - epsilmin;

   //
   // sample the energy rate of the created electron (or positron) 
   //
   //G4double epsil, screenvar, greject ;
   G4double  screenvar, greject ;

   G4double F10 = ScreenFunction1(screenmin) - FZ;
   G4double F20 = ScreenFunction2(screenmin) - FZ;
   G4double NormF1 = max(F10*epsilrange*epsilrange,0.); 
   G4double NormF2 = max(1.5*F20,0.);

   do {
        if ( NormF1/(NormF1+NormF2) > G4UniformRand() )
             { epsil = 0.5 - epsilrange*pow(G4UniformRand(), 0.333333);
               screenvar = screenfac/(epsil*(1-epsil));
               greject = (ScreenFunction1(screenvar) - FZ)/F10;
             } 
        else { epsil = epsilmin + epsilrange*G4UniformRand();
               screenvar = screenfac/(epsil*(1-epsil));
               greject = (ScreenFunction2(screenvar) - FZ)/F20;
             }

   } while( greject < G4UniformRand() );

   }   //  end of epsil sampling
   
   //
   // fixe charges randomly
   //

   G4double ElectTotEnergy, PositTotEnergy;
   if (G4UniformRand() > 0.5)
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
   //  universal distribution suggested by L. Urban 
   // (Geant3 manual (1993) Phys211),
   //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

   G4double u;
   const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

   if (9./(9.+d) >G4UniformRand()) u= - log(G4UniformRand()*G4UniformRand())/a1;
   else                            u= - log(G4UniformRand()*G4UniformRand())/a2;

   G4double TetEl = u*electron_mass_c2/ElectTotEnergy;
   G4double TetPo = u*electron_mass_c2/PositTotEnergy;
   G4double Phi  = twopi * G4UniformRand();
   G4double dxEl= sin(TetEl)*cos(Phi),dyEl= sin(TetEl)*sin(Phi),dzEl=cos(TetEl);
   G4double dxPo=-sin(TetPo)*cos(Phi),dyPo=-sin(TetPo)*sin(Phi),dzPo=cos(TetPo);
   
   //
   // kinematic of the created pair
   //
   // the electron and positron are assumed to have a symetric
   // angular distribution with respect to the Z axis along the parent photon.

   aParticleChange.SetNumberOfSecondaries(2); 

   G4double ElectKineEnergy = max(0.,ElectTotEnergy - electron_mass_c2);
   G4double localEnergyDeposit = 0.;

   if (ElectKineEnergy > fminimalEnergy)
     {
       G4ThreeVector ElectDirection (dxEl, dyEl, dzEl);
       ElectDirection.rotateUz(GammaDirection);   

       // create G4DynamicParticle object for the particle1  
       G4DynamicParticle* aParticle1= new G4DynamicParticle(
                        G4Electron::Electron(),ElectDirection,ElectKineEnergy);
       aParticleChange.AddSecondary(aParticle1); 
     }
   else
     { localEnergyDeposit += ElectKineEnergy;}   

   // the e+ is always created (even with Ekine=0) for further annihilation.

   G4double PositKineEnergy = max(0.,PositTotEnergy - electron_mass_c2);
   if (PositKineEnergy < fminimalEnergy)
     { localEnergyDeposit += PositKineEnergy; PositKineEnergy = 0.;}

   G4ThreeVector PositDirection (dxPo, dyPo, dzPo);
   PositDirection.rotateUz(GammaDirection);   

   // create G4DynamicParticle object for the particle2 
   G4DynamicParticle* aParticle2= new G4DynamicParticle(
                      G4Positron::Positron(),PositDirection,PositKineEnergy);
   aParticleChange.AddSecondary(aParticle2); 

   aParticleChange.ProposeLocalEnergyDeposit(localEnergyDeposit);

   //
   // Kill the incident photon 
   //

   aParticleChange.ProposeEnergy( 0. ); 
   aParticleChange.ProposeTrackStatus( fStopAndKill );

   //  Reset NbOfInteractionLengthLeft and return aParticleChange
   return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Element* G4GammaConversion52::SelectRandomAtom(
                                         const G4DynamicParticle* aDynamicGamma,
                                               G4Material* aMaterial)
{
  // select randomly 1 element within the material

  const G4int NumberOfElements            = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  if (NumberOfElements == 1) return (*theElementVector)[0];

  const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

  G4double PartialSumSigma = 0. ;
  G4double rval = G4UniformRand()/MeanFreePath;

  for ( G4int i=0 ; i < NumberOfElements ; i++ )
      { PartialSumSigma += NbOfAtomsPerVolume[i] *
                  GetCrossSectionPerAtom(aDynamicGamma, (*theElementVector)[i]);
        if (rval <= PartialSumSigma) return ((*theElementVector)[i]);
      }
  G4cout << " WARNING !!! - The Material '"<< aMaterial->GetName()
       << "' has no elements, NULL pointer returned." << G4endl;
  return NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4GammaConversion52::StorePhysicsTable(const G4ParticleDefinition* particle,
                                            const G4String& directory,
				                  G4bool          ascii)
{
  G4String filename;

  // store cross section table
  filename = GetPhysicsTableFileName(particle,directory,"CrossSection",ascii);
  if ( !theCrossSectionTable->StorePhysicsTable(filename, ascii) ){
    G4cout << " FAIL theCrossSectionTable->StorePhysicsTable in " << filename
           << G4endl;
    return false;
  }

  // store mean free path table
  filename = GetPhysicsTableFileName(particle,directory,"MeanFreePath",ascii);
  if ( !theMeanFreePathTable->StorePhysicsTable(filename, ascii) ){
    G4cout << " FAIL theMeanFreePathTable->StorePhysicsTable in " << filename
           << G4endl;
    return false;
  }

  G4cout << GetProcessName() << " for " << particle->GetParticleName()
         << ": Success to store the PhysicsTables in "
         << directory << G4endl;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
G4bool G4GammaConversion52::RetrievePhysicsTable(const G4ParticleDefinition* particle,
					       const G4String& directory,
				                     G4bool          ascii)
{
  // delete theCrossSectionTable and theMeanFreePathTable
  if (theCrossSectionTable != 0) {
    theCrossSectionTable->clearAndDestroy();
    delete theCrossSectionTable;
  }
  if (theMeanFreePathTable != 0) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }

  G4String filename;

  // retreive cross section table
  filename = GetPhysicsTableFileName(particle,directory,"CrossSection",ascii);
  theCrossSectionTable = new G4PhysicsTable(G4Element::GetNumberOfElements());
  if ( !G4PhysicsTableHelper::RetrievePhysicsTable(filename, ascii) ){
    G4cout << " FAIL theCrossSectionTable->RetrievePhysicsTable in " << filename
           << G4endl;
    return false;
  }

  // retreive mean free path table
  filename = GetPhysicsTableFileName(particle,directory,"MeanFreePath",ascii);
  theMeanFreePathTable = new G4PhysicsTable(G4Material::GetNumberOfMaterials());
  if ( !G4PhysicsTableHelper::RetrievePhysicsTable(filename, ascii) ){
    G4cout << " FAIL theMeanFreePathTable->RetrievePhysicsTable in " << filename
           << G4endl;
    return false;
  }

  G4cout << GetProcessName() << " for " << particle->GetParticleName()
         << ": Success to retrieve the PhysicsTables from "
         << directory << G4endl;
  return true;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GammaConversion52::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from a parametrisation. ";
           comments += "Good description from 1.5 MeV to 100 GeV for all Z. \n";
           comments += "        e+e- energies according Bethe-Heitler";

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from "
	           << G4BestUnit(LowestEnergyLimit, "Energy")
         << " to " << G4BestUnit(HighestEnergyLimit,"Energy")
         << " in " << NumbBinTable << " bins. \n";
  G4cout << "        WARNING: This process is obsolete and will be soon removed" 
	 << G4endl;
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
