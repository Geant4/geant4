// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4eplusAnnihilation.cc,v 1.3 2001-02-22 18:26:09 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// 10-01-97, crossection table + mean free path table, M.Maire
// 17-03-97, merge 'in fly' and 'at rest', M.Maire
// 23-03-97, protection in BuildPhysicsTable, M.Maire
// 31-08-98, new methods SetBining() and PrintInfo()
// 22-02-01, postStepDoIt: fStopButAlive instead of kineEnergy == 0.  
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eplusAnnihilation.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// constructor
 
G4eplusAnnihilation::G4eplusAnnihilation(const G4String& processName)
  : G4VRestDiscreteProcess (processName),
    theCrossSectionTable(NULL),
    theMeanFreePathTable(NULL),   
    LowestEnergyLimit ( 10*keV),      // initialization
    HighestEnergyLimit( 10*TeV),
    NumbBinTable(100)
 { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
// destructor
 
G4eplusAnnihilation::~G4eplusAnnihilation()
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


void G4eplusAnnihilation::SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins)
{
  LowestEnergyLimit = lowE; HighestEnergyLimit = highE; NumbBinTable = nBins;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
void G4eplusAnnihilation::BuildPhysicsTable(const G4ParticleDefinition& PositronType)

// Build microscopic total cross section tables and mean free path table
{
   G4double LowEdgeEnergy, Value;
   G4PhysicsLogVector* ptrVector;

// Build microscopic cross section tables for the e+e- annihilation

   if (theCrossSectionTable) {
          theCrossSectionTable->clearAndDestroy(); delete theCrossSectionTable; }

   theCrossSectionTable = new G4PhysicsTable( G4Element::GetNumberOfElements()) ;
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

// Build mean free path table for the e+e- annihilation

   if (theMeanFreePathTable) {
          theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable; }

   theMeanFreePathTable = new G4PhysicsTable( G4Material::GetNumberOfMaterials() );
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

G4double G4eplusAnnihilation::ComputeCrossSectionPerAtom
                              (G4double PositKinEnergy, G4double AtomicNumber)
 
// Calculates the microscopic cross section of annihilation into two photons
// from the Heilter formula.
// GEANT4 internal units.

{
 static const G4double pi_rcl2 = pi*classic_electr_radius*classic_electr_radius;

 G4double gama = 1. + PositKinEnergy/electron_mass_c2;
 G4double gama2 = gama*gama, sqgama2 = sqrt(gama2-1.);

 return pi_rcl2*AtomicNumber*((gama2+4*gama+1.)*log(gama+sqgama2) - (gama+3.)*sqgama2) 
                            /((gama2-1.)*(gama+1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
 
G4VParticleChange* G4eplusAnnihilation::PostStepDoIt(const G4Track& aTrack,
                                                    const G4Step&  aStep)
//
// The secondaries Gamma energies are sampled using the Heitler cross section.
//  
// A modified version of the random number techniques of Butcher & Messel is used 
//    (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : The initial electron is assumed free and at rest.
//          
// Note 2 : The annihilation processes producing one or more than two photons are
//          ignored, as negligible compared to the two photons process.         
 
{
   aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();
     

   const G4DynamicParticle* aDynamicPositron = aTrack.GetDynamicParticle();
   G4double PositKinEnergy = aDynamicPositron->GetKineticEnergy();
   G4ParticleMomentum PositDirection = aDynamicPositron->GetMomentumDirection();

   aParticleChange.Initialize(aTrack);

   // Do not make anything if particle is stopped, the annihilation then
   // should be performed by the AtRestDoIt!
   if (aTrack.GetTrackStatus() == fStopButAlive) return &aParticleChange;

   G4double gamam1 = PositKinEnergy/electron_mass_c2;
   G4double gama   = gamam1+1. , gamap1 = gamam1+2. , sqgrate = sqrt(gamam1/gamap1)/2. ,
                                                      sqg2m1  = sqrt(gamam1*gamap1);

   // limits of the energy sampling
   G4double epsilmin = 0.5 - sqgrate , epsilmax = 0.5 + sqgrate;
   G4double epsilqot = epsilmax/epsilmin;

   //
   // sample the energy rate of the created gammas 
   //
   G4double epsil, greject ;

   do {
        epsil = epsilmin*pow(epsilqot,G4UniformRand());
        greject = 1. - epsil + (2*gama*epsil-1.)/(epsil*gamap1*gamap1);
   } while( greject < G4UniformRand() );

   //
   // scattered Gamma angles. ( Z - axis along the parent positron)
   //
   
   G4double cost = (epsil*gamap1-1.)/(epsil*sqg2m1) , sint = sqrt((1.+cost)*(1.-cost));


   G4double phi  = twopi * G4UniformRand() ;
   G4double dirx = sint*cos(phi) , diry = sint*sin(phi) , dirz = cost;
 
   //
   // kinematic of the created pair
   //

   G4double LocalEnerDeposit = 0. ;
   aParticleChange.SetNumberOfSecondaries(2) ; 

   G4double TotalAvailableEnergy = PositKinEnergy + 2*electron_mass_c2;
   G4double Phot1Energy = epsil*TotalAvailableEnergy;

   G4double GammaCut= (G4Gamma::GetCutsInEnergy())[aMaterial->GetIndex()];

   if (Phot1Energy > GammaCut)
      {
        G4ThreeVector Phot1Direction ( dirx, diry, dirz );
        Phot1Direction.rotateUz(PositDirection);   
 
        // create G4DynamicParticle object for the particle1  
        G4DynamicParticle* aParticle1= new G4DynamicParticle (G4Gamma::Gamma(),
                                                 Phot1Direction, Phot1Energy);
        aParticleChange.AddSecondary( aParticle1 ) ; 
       }
    else
       { LocalEnerDeposit += Phot1Energy; }

    G4double Phot2Energy =(1.-epsil)*TotalAvailableEnergy; 

   if (Phot2Energy > GammaCut)
      {
        G4double Eratio = Phot1Energy/Phot2Energy;
        G4double PositP = sqrt(PositKinEnergy*(PositKinEnergy+2.*electron_mass_c2));
        G4ThreeVector Phot2Direction (-dirx*Eratio, -diry*Eratio,
                                         (PositP-dirz*Phot1Energy)/Phot2Energy); 
        Phot2Direction.rotateUz(PositDirection);
 
        // create G4DynamicParticle object for the particle2 
        G4DynamicParticle* aParticle2= new G4DynamicParticle (G4Gamma::Gamma(),
                                                 Phot2Direction, Phot2Energy);
        aParticleChange.AddSecondary( aParticle2 ) ; 
       }
    else
       { LocalEnerDeposit += Phot2Energy; }

   aParticleChange.SetLocalEnergyDeposit( LocalEnerDeposit ) ;

   //
   // Kill the incident positron 
   //

   aParticleChange.SetMomentumChange( 0., 0., 0. ) ;
   aParticleChange.SetEnergyChange( 0. ) ; 
   aParticleChange.SetStatusChange( fStopAndKill ) ;

   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4VParticleChange* G4eplusAnnihilation::AtRestDoIt(const G4Track& aTrack,
                                                  const G4Step&  aStep)
//
// Performs the e+ e- annihilation when both particles are assumed at rest.
// It generates two back to back photons with energy = electron_mass.
// The angular distribution is isotropic. 
// GEANT4 internal units
//
// Note : Effects due to binding of atomic electrons are negliged.
 
{
   aParticleChange.Initialize(aTrack);
   G4Material* aMaterial = aTrack.GetMaterial();

   aParticleChange.SetNumberOfSecondaries(2) ; 

   if (electron_mass_c2 > (G4Gamma::GetCutsInEnergy())[aMaterial->GetIndex()])
      {
        G4double cosTeta = 2*G4UniformRand()-1. , sinTeta = sqrt(1.-cosTeta*cosTeta);
        G4double Phi     = twopi * G4UniformRand() ;
        G4ThreeVector Direction (sinTeta*cos(Phi), sinTeta*sin(Phi), cosTeta);   
 
        aParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                                 Direction, electron_mass_c2) );
        aParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                                -Direction, electron_mass_c2) ); 

        aParticleChange.SetLocalEnergyDeposit(0.);
       }
   else
      { aParticleChange.SetLocalEnergyDeposit( 2*electron_mass_c2 ); }

   // Kill the incident positron 
   //
   aParticleChange.SetStatusChange( fStopAndKill );
      
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::PrintInfoDefinition()
{
  G4String comments = "Total cross section from Heilter formula (annihilation into 2 photons).\n";
           comments += "        gamma energies sampled according Heitler";
                     
  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestEnergyLimit,"Energy")
         << " to " << G4BestUnit(HighestEnergyLimit,"Energy") 
         << " in " << NumbBinTable << " bins. \n";
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
