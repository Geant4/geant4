// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4ionLowEnergyIonisation physics process -----
//                by Vladimir Ivanchenko, 6 September 1999 
//                was made on the base of G4hLowEnergyIonisation class
// ************************************************************
// It is the extention of the ionisation process for the slow 
// charged ions.
// ************************************************************
//  6 September 1999 V.Ivanchenko create
// 30 September 1999 V.Ivanchenko minor upgrade
// 20 January   2000 V.Ivanchenko minor bag fixed
// ------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ionLowEnergyIonisation.hh"
#include "G4UnitsTable.hh"
#include "G4EnergyLossTables.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionLowEnergyIonisation::G4ionLowEnergyIonisation(const G4String& processName)
  : G4hLowEnergyIonisation(processName)
{ 
  LowestKineticEnergy = 10.*eV ;
  HighestKineticEnergy = 100.*TeV ;
  TotBin = 200 ;
  MassRatio = 1.0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionLowEnergyIonisation::~G4ionLowEnergyIonisation() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionLowEnergyIonisation::SetIonDefinition(G4ParticleDefinition* theIonType)
{
  theParticle  = theIonType ;
  MassRatio = proton_mass_c2/(theParticle->GetPDGMass()) ;
  Charge    = (theParticle->GetPDGCharge())/eplus ;
  G4cout << "New ion with Q = " << Charge << "; MassR = " << MassRatio << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionLowEnergyIonisation::GetLowEnergyForParametrisation(const G4Material* material) 

{
  // The low limit of paramerisation of ionisation energy from: 
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985
  // Below this limit the free electron gas model is used

  // hadrons or ions with charge +-1
  if(Charge < 1.5) return ParamLowEnergy/MassRatio ;

  // helium or ions with charge = +2
  if(Charge < 2.5) return ParamLowEnergy ;

  // get elements in the actual material,
  const G4ElementVector* theElementVector = material->GetElementVector() ;
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector() ;
  const G4int NumberOfElements = material->GetNumberOfElements() ;
  G4double Z = 0.0, Norm = 0.0 ; 
  
  // only 1 element in the material
  if( 1 == NumberOfElements ) {
    Z = material->GetZ() ;

  //  loop for the elements in the material
  //  to find out average value of Z
  } else {
    for (G4int iel=0; iel<NumberOfElements; iel++)
      {
        const G4Element* element = (*theElementVector)(iel) ;
        G4double Z2 = element->GetZ() ;
        const G4double W2 = theAtomicNumDensityVector[iel] ;
        Norm += W2 ;
        Z    += Z2 * W2 ;
      }
    Z  /= Norm ;
  }
  G4double E1 = 3.25 * keV ;
  G4double E2 = 25.0 * keV / pow(Z, 0.667) ;
  E1 = G4std::max (E1, E2) ;
  E1 = G4std::max(ParamLowEnergy, E1) / MassRatio ; 
  return E1 ; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionLowEnergyIonisation::GetIonLossWithFluct(const G4DynamicParticle* aParticle,
                                                             G4Material* aMaterial,
                                                             G4double    MeanLoss)

//  calculate actual loss from the mean loss
//  The model used to get the fluctuation is the same as in Glandz in Geant3.
{
  G4double ChargeSquare = Charge*Charge ;

  G4double loss = GetLossWithFluct(aParticle, aMaterial, MeanLoss/ChargeSquare) * ChargeSquare ;

  return loss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionLowEnergyIonisation::PrintInfoDefinition()
{
  G4String comments = "  Knock-on electron cross sections . ";
  comments += "\n         Good description above the mean excitation energy.\n";
  comments += "         delta ray energy sampled from  differential Xsection.";
  
  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
							  "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins."
         << "\n        Low energy losses approximation is taken from  " << DEDXtable
         << "\n        from " << G4BestUnit(ParamLowEnergy,"Energy")
         << " to " << G4BestUnit(ParamHighEnergy,"Energy") << "." << G4endl ;
  if(nStopping) {
    G4cout << "        Simulation of nuclear stopping is switched on.  \n" << G4endl ; 
  }
}











































