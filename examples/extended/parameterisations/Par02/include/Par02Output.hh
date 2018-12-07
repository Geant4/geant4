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
/// \file Par02Output.hh
/// \brief Definition of the Par02Output class

#ifndef PAR02_OUTPUT_H
#define PAR02_OUTPUT_H

#include "G4ThreeVector.hh"
#include "globals.hh"

/// Handling the saving to the file.
///
/// A singleton class that manages creation, writing to and closing of the
/// Root output file.
/// @author Anna Zaborowska

class Par02Output {
  public:
    
    /// Indicates to which ntuple to save the information.
    enum SaveType { eNoSave, eSaveMC, eSaveTracker, eSaveEMCal, eSaveHCal };

    /// Allows the access to the unique Par02Output object.
    /// @return A pointer to the Par02Output class.
    static Par02Output* Instance();
    
    /// Sets the file name of the output root file.
    /// @param name The name of the file.
    void SetFileName( G4String name );
    
    /// Gets the file name of the output root file.
    /// @return The name of the file.
    G4String GetFileName();
    
    /// Sets fFileNameWithRunNo that indicates whether to add the run number
    /// to the file name.
    /// @param app If add the run number.
    void AppendName( G4bool app );
    
    /// Calls the G4AnalysisManager::Instance(). It sets the file name of the
    /// output file based on fFileName and fFileNameWithRunNo and opens the file.
    /// @param runID A run number (to be added to file name if fFileNameWithRunNo
    ///              is true).
    void StartAnalysis( G4int runID );
    
    /// Calls the G4AnalysisManager::Instance(). 
    /// It writes to the output file and close it.
    void EndAnalysis();
    
    /// Creates Ntuples used to store information about particle (its ID, PDG code,
    /// energy deposits, etc.). To be called for each event in Par02EventAction.
    void CreateNtuples();
    
    /// Creates histograms to combine information from all the events in the run.
    /// To be called for each run in Par02RunAction.
    void CreateHistograms();
    
    /// Saves the information about the particle (track).
    /// @param aWhatToSave enum indicating what kind of information to store 
    ///                    (in which ntuple).
    /// @param aPartID A unique ID within event (taken Geant TrackID).
    /// @param aPDG A PDG code of a particle.
    /// @param aVector A vector to be stored (particle momentum in tracker or
    ///                position of energy deposit in calorimeter).
    /// @param aResolution A resolution of the detector that was used.
    /// @param aEfficiency An efficiency of the detector that was used.
    /// @param aEnergy An energy deposit (for calorimeters only: 
    ///                Par02Output::SaveType::eEMCal or Par02Output::SaveType::eHCal).
    void SaveTrack( SaveType aWhatToSave, G4int aPartID,  G4int aPDG,
                    G4ThreeVector aVector, G4double aResolution = 0,
                    G4double aEfficiency = 1, G4double aEnergy = 0 ) ;
    
    /// Fills the histogram.
    /// @param HNo Number of a histogram (decided by the order of creation
    ///            in CreateHistograms(), the first one is 0).
    /// @param value A value to be filled into the histogram.
    void FillHistogram( G4int HNo, G4double value ) const;

    ~Par02Output();

  protected:
    
    /// A default, protected constructor (due to singleton pattern).
    Par02Output();

  private:
    
    /// The pointer to the only Par02Output class object.
    static Par02Output* fPar02Output;

    /// Current ntuple Id 
    static G4ThreadLocal G4int fCurrentNtupleId;
    
    /// A name of the output root file.
    G4String fFileName;
    
    /// If true, a run number should be added to the file. Default: false.
    G4bool fFileNameWithRunNo;
    
    /// A control value of particle ID to ensure that data saved to various ntuples
    /// match the same particle. It is set when Monte Carlo information is saved
    /// and checked for all the detectors.
    static G4ThreadLocal G4int fCurrentID;
};

#endif

