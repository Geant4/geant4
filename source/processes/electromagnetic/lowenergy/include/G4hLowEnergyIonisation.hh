// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hLowEnergyIonisation.hh,v 1.2 1999-07-23 16:37:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hLowEnergyIonisation physics process -----
//                by Vladimir Ivanchenko, 14 July 1999 
// ************************************************************
// It is the first implementation of the NEW IONISATION     
// PROCESS. ( delta rays + continuous energy loss)
// It calculates the ionisation for slow charged hadrons.      
// ************************************************************
// ------------------------------------------------------------
 
#ifndef G4hLowEnergyIonisation_h
#define G4hLowEnergyIonisation_h 1
 
#include "G4hIonisation.hh"
  
class G4hLowEnergyIonisation : public G4hIonisation 
 
{
  public:
 
     G4hLowEnergyIonisation(const G4String& processName = "hLowEIoni"); 

    ~G4hLowEnergyIonisation();

    void BuildLossTable(const G4ParticleDefinition& aParticleType);

    void SetDEDXTableName(const G4String& dedxTable);
    
    G4double GetZieglerLoss(const G4Material* material, G4double KinEnergy);

    G4double GetBetheBlochLoss(const G4Material* material, G4double KinEnergy);

    G4double GetFreeElectronGasLoss(G4double paramA, G4double KinEnergy);

    G4double GetStoppingPower1977H(G4int iz, G4double E);

    G4double GetStoppingPower1977He(G4int iz, G4double E);

    G4double GetStoppingPower1977n(G4double Z1, G4double Z2, 
                                   G4double M1, G4double M2, G4double E);

    G4double GetUrbanModel(const G4Element* element, G4double KinEnergy);

    G4double GetDeltaRaysEnergy(const G4Material* material, G4double KinEnergy);

    G4double GetChemicalFactor(const G4Material* material, G4double KinEnergy);

    void PrintInfoDefinition();

  private:

  // hide assignment operator 
    G4hLowEnergyIonisation & operator=(const G4hLowEnergyIonisation &right);
    G4hLowEnergyIonisation(const G4hLowEnergyIonisation&);

  private:
  //  private data members ...............................

    G4PhysicsTable* theMeanFreePathTable;

    // LowestKineticEnergy = lower limit of particle kinetic energy
    // HighestKineticEnergy = upper limit of particle kinetic energy 
    // TotBin = number of bins 
    //  ---------in the energy ionisation loss table-------------------
    G4double LowestKineticEnergy;
    G4double HighestKineticEnergy;
    G4double ZieglerLowEnergy;
    G4double ZieglerHighEnergy;
    G4String DEDXtable;
    G4int TotBin;

    // constants needed for the energy loss calculation

    const G4double twoln10;
    const G4double Factor;
    const G4double bg2lim;
    const G4double taulim;          // energy to start to switch off shell corrections
    const G4double RateMass;

    // particles , cuts in kinetic energy ........
    const G4Electron* theElectron;
    const G4Proton* theProton;
    const G4AntiProton* theAntiProton;

    const G4double* DeltaCutInKineticEnergy ; 
 
    G4double DeltaCutInKineticEnergyNow ;

    G4double ProtonMassAMU;
    G4double ZieglerFactor; // Factor to convert the Stopping Power 
                            // unit [ev/(10^15 atoms/cm^2]
                            // into the Geant4 dE/dx unit
 
};
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4hLowEnergyIonisation::GetFreeElectronGasLoss(G4double paramA, G4double tau)
{
     G4double ionl = paramA * sqrt(tau) ;
     return ionl ;
}

#endif
 







