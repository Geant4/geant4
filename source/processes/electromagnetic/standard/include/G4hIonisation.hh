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
// $Id: G4hIonisation.hh,v 1.21 2003-01-17 18:55:43 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------- G4hIonisation physics process -------------------------------
//                 by Laszlo Urban, 30 May 1997 
// -----------------------------------------------------------------------------
//
// corrected by L.Urban on 24/09/97
// corrected by L.Urban on 13/01/98
// bugs fixed by L.Urban on 02/02/99
// 10/02/00 modifications , new e.m. structure, L.Urban
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 14-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 19-09-01 come back to previous process name "hIoni"
// 29-10-01 all static functions no more inlined
// 15-01-03 Migrade to cut per region (V.Ivanchenko)
//
// -----------------------------------------------------------------------------

// Class description
//
// This class manages the ionisation process for hadrons.
// it inherites from G4VContinuousDiscreteProcess via G4VhEnergyLoss.
//
// Class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4hIonisation_h
#define G4hIonisation_h 1

#include "G4VhEnergyLoss.hh"
#include "G4MaterialCutsCouple.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4hIonisation : public G4VhEnergyLoss

{
  public:   // with description

    G4hIonisation(const G4String& processName = "hIoni");

   ~G4hIonisation();

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

    G4bool StorePhysicsTable(G4ParticleDefinition* ,
		             const G4String& directory, G4bool);
      // store eLoss and MeanFreePath tables into an external file
      // specified by 'directory' (must exist before invokation)

    G4bool RetrievePhysicsTable(G4ParticleDefinition* ,
			        const G4String& directory, G4bool);
      // retrieve eLoss and MeanFreePath tables from an external file
      // specified by 'directory'

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

    G4VParticleChange *PostStepDoIt(const G4Track& track,
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

  protected:

    virtual G4double SecondaryEnergyThreshold(size_t index);

    G4PhysicsTable* theMeanFreePathTable;

  private:

    // hide assignment operator
    G4hIonisation & operator=(const G4hIonisation &right);
    G4hIonisation(const G4hIonisation&);

  private:

    static G4double LowerBoundLambda;      // bining for lambda table
    static G4double UpperBoundLambda;
    static G4int    NbinLambda;

    G4double LowestKineticEnergy;          // binning for dE/dx table
    G4double HighestKineticEnergy;
    G4int    TotBin;

    G4double Tmincut;                      // min energy of d-rays
    G4double initialMass;                  // mass for Lambda table

    const G4std::vector<G4double>* secondaryEnergyCuts;

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

#include "G4hIonisation.icc"

#endif








