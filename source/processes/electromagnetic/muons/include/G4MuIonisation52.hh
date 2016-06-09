//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4MuIonisation52.hh,v 1.2 2004/11/10 08:49:09 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// --------------- G4MuIonisation52 physics process ------------------------------
//                 by Laszlo Urban, September 1997
// -----------------------------------------------------------------------------
//
// 10-02-00 modifications , new e.m. structure, L.Urban
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 29-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 19-09-01 come back to the old process name 'MuIoni'
// 29-10-01 all static functions no more inlined (mma)
// 16-01-03 Migrade to cut per region (V.Ivanchenko)
// 08-08-03 This class is frozen at the release 5.2 (V.Ivanchenko)
// 08-11-04 Remove interface of Store/Retrieve tables (V.Ivantchenko)
//
// -----------------------------------------------------------------------------

// Class description
//
// This class manages the ionisation process for muons.
// it inherites from G4VContinuousDiscreteProcess via G4VMuEnergyLoss.
//
// Class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MuIonisation52_h
#define G4MuIonisation52_h 1

#include "G4VMuEnergyLoss.hh"
#include "G4MaterialCutsCouple.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MuIonisation52 : public G4VMuEnergyLoss

{
  public:   // with description

    G4MuIonisation52(const G4String& processName = "MuIoni");

   ~G4MuIonisation52();

    G4bool IsApplicable(const G4ParticleDefinition&);
      // return true for charged particles, false otherwise

    void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
      // this function overloads a virtual function of the base class.
      // It is invoked by the G4ParticleWithCuts::SetCut() method.
      // It invokes BuildLambdaTable(), BuildLossTable(), BuildDEDXTable()

    void BuildLossTable(const G4ParticleDefinition& aParticleType);
      // build the dE/dx tables due to the ionisation, for every materials.
      // (restricted stopping power, Bethe-Bloch formula)

    void BuildLambdaTable(const G4ParticleDefinition& aParticleType);
      // build mean free path tables for the delta rays production.
      // the tables are built for every materials.

    virtual void PrintInfoDefinition();
      // Print few lines of informations about the process: validity range,
      // origine ..etc..
      // Invoked by BuildPhysicsTable().

    G4double GetMeanFreePath(const G4Track& track,
                             G4double previousStepSize,
                             G4ForceCondition* condition );
      // It returns the MeanFreePath of the process for the current track :
      // (energy, material)
      // The previousStepSize and G4ForceCondition* are not used.
      // This function overloads a virtual function of the base class.
      // It is invoked by the ProcessManager of the Particle.

    G4VParticleChange* PostStepDoIt(const G4Track& track,
                                          const G4Step& Step  );
      // It computes the final state of the process (at end of step),
      // returned as a ParticleChange object.
      // This function overloads a virtual function of the base class.
      // It is invoked by the ProcessManager of the Particle.


  protected:   // with description

    virtual G4double ComputeRestrictedMeandEdx(
                            const G4ParticleDefinition& aParticleType,
                            G4double KineticEnergy,
                            const G4Material* material,
                            G4double DeltaThreshold);
      // computes restricted mean dE/dx in Geant4 internal units.

    virtual G4double ComputeCrossSectionPerAtom(
                            const G4ParticleDefinition& aParticleType,
                            G4double KineticEnergy,
                            G4double AtomicNumber,
			    G4double DeltaThreshold);
      // computes total cross section per atom in Geant4 internal units.

    virtual G4double ComputeDifCrossSectionPerAtom(
                            const G4ParticleDefinition& aParticleType,
                            G4double KineticEnergy,
                            G4double AtomicNumber,
			    G4double KnockonEnergy);
      // computes differential cross section per atom.

  protected:

    virtual G4double SecondaryEnergyThreshold(size_t index);

    G4PhysicsTable* theMeanFreePathTable;

  private:

    // hide assignment operator
    G4MuIonisation52 & operator=(const G4MuIonisation52 &right);
    G4MuIonisation52(const G4MuIonisation52&);

  private:

    static G4double LowerBoundLambda;      // bining for lambda table
    static G4double UpperBoundLambda;
    static G4int    NbinLambda;

    G4double LowestKineticEnergy;         // binning for dE/dx table
    G4double HighestKineticEnergy;
    G4int    TotBin;

    const std::vector<G4double>* secondaryEnergyCuts;

  public:  // with description

    static void SetLowerBoundLambda(G4double val);
    static void SetUpperBoundLambda(G4double val);
    static void SetNbinLambda(G4int n);
        // set the parameters of the mean free path table.

    static G4double GetLowerBoundLambda();
    static G4double GetUpperBoundLambda();
    static G4int GetNbinLambda();
      // get the parameters of the mean free path table.
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuIonisation52.icc"

#endif








