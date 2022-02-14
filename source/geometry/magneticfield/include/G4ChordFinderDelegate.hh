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
// G4ChordFinderDelegate
//
// Class description:
//
// Implementation of common algorithm of finding step size 
// with distance to chord less then provided value. 

// Created: D.Sorokin
// --------------------------------------------------------------------
#ifndef G4CHORD_FINDER_DELEGATE_HH
#define G4CHORD_FINDER_DELEGATE_HH

#include "G4VIntegrationDriver.hh"

template <class Driver>
class G4ChordFinderDelegate
{
  public:

    virtual ~G4ChordFinderDelegate();

    G4double AdvanceChordLimitedImpl(G4FieldTrack& track,
                                     G4double hstep,
                                     G4double eps,
                                     G4double chordDistance);
    void ResetStepEstimate();

    void TestChordPrint(G4int noTrials, 
                        G4int lastStepTrial, 
                        G4double dChordStep, 
                        G4double fDeltaChord,
                        G4double nextStepTrial);

    // Get statistics about number of calls & trials in FindNextChord
    G4int GetNoCalls(); 
    G4int GetNoTrials();        // Total number of trials
    G4int GetNoMaxTrials();     // Maximum # of trials for one call

    // Parameters for  performance ... change with great care
    void SetFractions_Last_Next(G4double fractLast = 0.90, 
                                G4double fractNext = 0.95); 
    void SetFirstFraction(G4double fractFirst);

    // Printing for monitoring ...
    G4double GetFirstFraction();         // Originally 0.999
    G4double GetFractionLast();          // Originally 1.000
    G4double GetFractionNextEstimate();  // Originally 0.980

    G4double GetLastStepEstimateUnc(); 
    void SetLastStepEstimateUnc(G4double stepEst);

    void  StreamDelegateInfo( std::ostream& os ) const;
     // Write out the parameters / state of the driver
   
  private:

    Driver& GetDriver();

    G4double FindNextChord(const G4FieldTrack& yStart,
                           G4double stepMax,
                           G4double epsStep,
                           G4double chordDistance,
                           G4FieldTrack& yEnd, // Endpoint
                           G4double& dyErrPos, // Error of endpoint
                           G4double& pStepForAccuracy);

    G4double NewStep(G4double stepTrialOld, 
                     G4double dChordStep, // Curr. dchord achieved
                     G4double fDeltaChord,
                     G4double& stepEstimate_Unconstrained);

    void AccumulateStatistics(G4int noTrials);

    void PrintStatistics();

    G4double fFirstFraction = 0.999;
    G4double fFractionLast = 1.0;
    G4double fFractionNextEstimate = 0.98;
    G4double fLastStepEstimate_Unconstrained = DBL_MAX;

    G4int fTotalNoTrials = 0;
    G4int fNoCalls = 0;
    G4int fmaxTrials = 0;
};

#include "G4ChordFinderDelegate.icc"

#endif
