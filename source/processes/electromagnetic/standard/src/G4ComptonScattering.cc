// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ComptonScattering.cc,v 1.1 1999-01-07 16:11:19 gunter Exp $
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
//      ------------ G4ComptonScattering physics process --------
//                   by Michel Maire, April 1996
// **************************************************************
// 28-05-96, DoIt() small change in ElecDirection, by M.Maire
// 10-06-96, simplification in ComputeMicroscopicCrossSection(), by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 13-09-96, small changes in DoIt for better efficiency. Thanks to P.Urban
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 05-03-97, new Physics scheme, M.Maire
// 28-03-97, protection in BuildPhysicsTable, M.Maire
// 07-04-98, remove 'tracking cut' of the scattered gamma, MMa
// 04-06-98, in DoIt, secondary production condition: range>min(threshold,safety)
// 13-08-98, new methods SetBining()  PrintInfo()
// 15-12-98, cross section=0 below 10 keV
// --------------------------------------------------------------

#include "G4ComptonScattering.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// constructor
 
G4ComptonScattering::G4ComptonScattering(const G4String& processName)
  : G4VDiscreteProcess (processName),
    theCrossSectionTable(NULL),
    theMeanFreePathTable(NULL),  
    LowestEnergyLimit ( 10*keV),              // initialization
    HighestEnergyLimit(100*GeV),
    NumbBinTable(100)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// destructor
 
G4ComptonScattering::~G4ComptonScattering()
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

void G4ComptonScattering::SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins)
{
  LowestEnergyLimit = lowE; HighestEnergyLimit = highE; NumbBinTable = nBins;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
void G4ComptonScattering::BuildPhysicsTable(const G4ParticleDefinition& GammaType)

// Build microscopic cross section table and mean free path table
{
   G4double LowEdgeEnergy, Value;
   G4PhysicsLogVector* ptrVector;

// Build microscopic cross section tables for the Compton Scattering process

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
             Value = ComputeCrossSectionPerAtom( LowEdgeEnergy, AtomicNumber);  
             ptrVector->PutValue( i , Value ) ;
           }

        theCrossSectionTable->insertAt( J , ptrVector ) ;

      }

// Build mean free path table for the Compton Scattering process

   if (theMeanFreePathTable) {
             theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable; }
 
   theMeanFreePathTable = new G4PhysicsTable(G4Material::GetNumberOfMaterials());
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

G4double G4ComptonScattering::ComputeCrossSectionPerAtom
                              (G4double GammaEnergy, G4double Z)
 
// Calculates the microscopic cross section in GEANT4 internal units.
// A parametrized formula from L. Urban is used to estimate the total cross section.
// It gives a good description of the data from 10 keV to 100/Z GeV.
 
{
 G4double CrossSection = 0.0 ;
 if ( Z < 1. )                     return CrossSection;
 if ( GammaEnergy < 10.*keV      ) return CrossSection;
 if ( GammaEnergy > (100.*GeV/Z) ) return CrossSection;

 static const G4double a = 20.0 , b = 230.0 , c = 440.0 ;
  
 static const G4double
  d1= 2.7965e-1*barn, d2=-1.8300e-1*barn, d3= 6.7527   *barn, d4=-1.9798e+1*barn,
  e1= 1.9756e-5*barn, e2=-1.0205e-2*barn, e3=-7.3913e-2*barn, e4= 2.7079e-2*barn,
  f1=-3.9178e-7*barn, f2= 6.8241e-5*barn, f3= 6.0480e-5*barn, f4= 3.0274e-4*barn;
     
 G4double p1Z = Z*(d1 + e1*Z + f1*Z*Z), p2Z = Z*(d2 + e2*Z + f2*Z*Z),
          p3Z = Z*(d3 + e3*Z + f3*Z*Z), p4Z = Z*(d4 + e4*Z + f4*Z*Z);

 G4double X = GammaEnergy / electron_mass_c2 ;

 return CrossSection = p1Z*log(1.+2*X)/X
                       + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4VParticleChange* G4ComptonScattering::PostStepDoIt(const G4Track& aTrack,
                                                    const G4Step&  aStep)
//
// The scattered gamma energy is sampled according to Klein - Nishina formula.
// The random number techniques of Butcher & Messel are used (Nuc Phys 20(1960),15).
// GEANT4 internal units
//
// Note : Effects due to binding of atomic electrons are negliged.
 
{
   aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();

   const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
   G4double GammaEnergy0 = aDynamicGamma->GetKineticEnergy();
   G4double E0_m = GammaEnergy0 / electron_mass_c2 ;

   G4ParticleMomentum GammaDirection0 = aDynamicGamma->GetMomentumDirection();


   //
   // sample the energy rate of the scattered gamma 
   //
 
   G4double epsilon, epsilonsq, onecost, sint2, greject ;

   G4double epsilon0 = 1./(1. + 2*E0_m) , epsilon0sq = epsilon0*epsilon0;
   G4double alpha1   = - log(epsilon0)  , alpha2 = 0.5*(1.- epsilon0sq);

   do {
       if ( alpha1/(alpha1+alpha2) > G4UniformRand() )
            { epsilon   = exp(-alpha1*G4UniformRand());  // pow(epsilon0,G4UniformRand())
              epsilonsq = epsilon*epsilon; }
       else {
             epsilonsq = epsilon0sq + (1.- epsilon0sq)*G4UniformRand();
             epsilon   = sqrt(epsilonsq);
       };
       onecost = (1.- epsilon)/(epsilon*E0_m);
       sint2   = onecost*(2.-onecost);
       greject = 1. - epsilon*sint2/(1.+ epsilonsq);
   } while (greject < G4UniformRand());
 
   //
   // scattered gamma angles. ( Z - axis along the parent gamma)
   //

   G4double cosTeta = 1. - onecost , sinTeta = sqrt (sint2);
   G4double Phi     = twopi * G4UniformRand() ;
   G4double dirx = sinTeta*cos(Phi) , diry = sinTeta*sin(Phi) , dirz = cosTeta ;

   //
   // update G4VParticleChange for the scattered gamma 
   //
   
   G4ThreeVector GammaDirection1 ( dirx,diry,dirz );
   GammaDirection1.rotateUz(GammaDirection0);
   aParticleChange.SetMomentumChange( GammaDirection1 ) ;
   G4double GammaEnergy1 = epsilon*GammaEnergy0;
   if (GammaEnergy1 > 0.)
     {
       aParticleChange.SetEnergyChange( GammaEnergy1 ) ;
     }
   else
     {    
       aParticleChange.SetEnergyChange(0.) ;
       aParticleChange.SetStatusChange(fStopAndKill);
     }
       
   //
   // kinematic of the scattered electron
   //

   G4double ElecKineEnergy = GammaEnergy0 - GammaEnergy1 ;
      
     // condition changed!
     if((G4EnergyLossTables::GetRange(G4Electron::Electron(),
          ElecKineEnergy,aMaterial)>aStep.GetPostStepPoint()->GetSafety())
         ||
        (ElecKineEnergy > 
        (G4Electron::Electron()->GetCutsInEnergy())[aMaterial->GetIndex()]))
      {
        G4double ElecMomentum = sqrt(ElecKineEnergy*(ElecKineEnergy+2.*electron_mass_c2));
        G4ThreeVector ElecDirection (
        (GammaEnergy0*GammaDirection0 - GammaEnergy1*GammaDirection1)*(1./ElecMomentum) );
 
        // create G4DynamicParticle object for the electron.  
        G4DynamicParticle* aElectron= new G4DynamicParticle (G4Electron::Electron(),
                                                        ElecDirection, ElecKineEnergy) ;

        aParticleChange.SetNumberOfSecondaries(1) ;
        aParticleChange.AddSecondary( aElectron ) ;
        aParticleChange.SetLocalEnergyDeposit (0.) ; 

       }
    else
       {
         aParticleChange.SetNumberOfSecondaries(0) ;
         aParticleChange.SetLocalEnergyDeposit (ElecKineEnergy) ;
       }

   //  Reset NbOfInteractionLengthLeft and return aParticleChange
   return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ComptonScattering::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from a parametrisation(L.Urban). ";
           comments += "Good description from 10 KeV to (100/Z) GeV. \n";
           comments += "        Scattered gamma energy according Klein-Nishina.";
                     
  G4cout << endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestEnergyLimit,"Energy")
         << " to " << G4BestUnit(HighestEnergyLimit,"Energy") 
         << " in " << NumbBinTable << " bins. \n";
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
