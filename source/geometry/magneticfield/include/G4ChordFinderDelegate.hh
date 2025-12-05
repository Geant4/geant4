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
// with distance to chord less than provided value. 

// Author: Dmitry Sorokin (CERN, Google Summer of Code 2017), 12.09.2018
// --------------------------------------------------------------------
#ifndef G4CHORD_FINDER_DELEGATE_HH
#define G4CHORD_FINDER_DELEGATE_HH

#include <iomanip>
#include "G4VIntegrationDriver.hh"

/**
 * @brief G4ChordFinderDelegate is a templated class for a common algorithm
 * of finding step size with distance to the chord less than the provided value.
 */

template <class Driver>
class G4ChordFinderDelegate
{
  public:

    /**
     * Virtual Destructor.
     */
    virtual ~G4ChordFinderDelegate();

    /**
     * Computes the step to take, based on chord limits.
     *  @param[in,out] track The current track in field.
     *  @param[in] hstep Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] chordDistance Maximum sagitta distance.
     *  @returns The length of step taken.
     */
    G4double AdvanceChordLimitedImpl(G4FieldTrack& track,
                                     G4double hstep,
                                     G4double eps,
                                     G4double chordDistance);

    /**
     * Resets last step estimate to DBL_MAX.
     */
    void ResetStepEstimate();

    /**
     * Getter and setter for step estimate.
     */
    G4double GetLastStepEstimateUnc(); 
    void SetLastStepEstimateUnc(G4double stepEst);

    /**
     * Gets statistics about number of calls & trials in FindNextChord().
     */
    G4int GetNoCalls(); 
    G4int GetNoTrials();        // Total number of trials
    G4int GetNoMaxTrials();     // Maximum # of trials for one call

    /**
     * Setters of performance parameters... change with great care!
     */
    void SetFractions_Last_Next(G4double fractLast = 0.90, 
                                G4double fractNext = 0.95); 
    void SetFirstFraction(G4double fractFirst);

    /**
     * Printing for monitoring ...
     */
    G4double GetFirstFraction();         // Originally 0.999
    G4double GetFractionLast();          // Originally 1.000
    G4double GetFractionNextEstimate();  // Originally 0.980

    /**
     * Writes out to stream the parameters/state of the driver.
     */
    void StreamDelegateInfo( std::ostream& os ) const;
   
    /**
     * statistics printout for testing.
     */
    void TestChordPrint(G4int noTrials, 
                        G4int lastStepTrial, 
                        G4double dChordStep, 
                        G4double fDeltaChord,
                        G4double nextStepTrial);
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

  private:

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
