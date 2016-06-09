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
// $Id: G4eIonisation52.hh,v 1.2 2004/11/10 08:53:19 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//--------------- G4eIonisation52 physics process --------------------------------
//                by Laszlo Urban, 20 March 1997
//
// 10-02-00 modifications , new e.m. structure, L.Urban
// 03-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 13-08-01 new function ComputeRestrictedMeandEdx() (mma)
// 19-09-01 come back to previous ProcessName "eIoni"
// 29-10-01 all static functions no more inlined (mma)
// 15-01-03 Migrade to cut per region (V.Ivanchenko)
// 08-08-03 This class is frozen at the release 5.2 (V.Ivanchenko)
// 09-11-04 Remove Store/Retrieve tables (V.Ivantchenko)
//
//------------------------------------------------------------------------------

// Class description
//
// This class manages the ionisation process for e-/e+
// it inherites from G4VContinuousDiscreteProcess via G4VeEnergyLoss.
//
// Class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4eIonisation52_h
#define G4eIonisation52_h 1

#include "G4VeEnergyLoss.hh"
#include "G4MaterialCutsCouple.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4eIonisation52 : public G4VeEnergyLoss

{
  public:   // with description

    G4eIonisation52(const G4String& processName = "eIoni");

   ~G4eIonisation52();

    G4bool IsApplicable(const G4ParticleDefinition&);
      // return true for e+/e-, false otherwise

    void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
      // this function overloads a virtual function of the base class.
      // It is invoked by the G4ParticleWithCuts::SetCut() method.
      // It invokes BuildLambdaTable(), BuildLossTable(), BuildDEDXTable()

    void BuildLossTable(const G4ParticleDefinition& aParticleType);
      // build the dE/dx tables due to the ionisation, for every materials.
      // (restricted stopping power, Berger-Seltzer formula)

    void BuildLambdaTable(const G4ParticleDefinition& aParticleType);
      // build mean free path tables for the delta rays production.
      // the tables are built for every materials.

    void PrintInfoDefinition();
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

    G4VParticleChange *PostStepDoIt(const G4Track& track,
                                    const G4Step& Step );
      // It computes the final state of the process (at end of step),
      // returned as a ParticleChange object.
      // This function overloads a virtual function of the base class.
      // It is invoked by the ProcessManager of the Particle.

    G4double GetLambda(G4double KineticEnergy, const G4MaterialCutsCouple* couple);
      // It returns the MeanFreePath of the process for a (energy, material)

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

  protected:

    virtual G4double SecondaryEnergyThreshold(size_t index);

    G4PhysicsTable* theMeanFreePathTable;

  private:

    // hide assignment operator
    G4eIonisation52 & operator=(const G4eIonisation52 &right);
    G4eIonisation52(const G4eIonisation52&);

  private:

    static G4double LowerBoundLambda;    // binning for mean free path table
    static G4double UpperBoundLambda;
    static G4int    NbinLambda;

    G4double LowestKineticEnergy;        // binning for dE/dx table
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

#include "G4eIonisation52.icc"

#endif

