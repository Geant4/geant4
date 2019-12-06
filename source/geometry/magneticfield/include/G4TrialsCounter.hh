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
// G4TrialsCounter
//
// Class description:
//
// Keeps statistics of the number of trials, including the Maximum
// and how many times it was reached.

// Author: John Apostolakis, CERN - 08.12.2006   
// -------------------------------------------------------------------
#ifndef G4TRIALS_COUNTER_HH
#define G4TRIALS_COUNTER_HH

#include "G4Types.hh"
#include "G4String.hh"

class  G4TrialsCounter
{
  public:  // with description

    G4TrialsCounter( const G4String& nameStats,
                     const G4String& description, G4bool printOnExit = false ); 
   ~G4TrialsCounter();

    inline void AccumulateCounts( G4int noTrials ); 
       //  Add this number to stats
    void ClearCounts(); 
       //  Reset all counts
    G4int ReturnTotals( G4int& calls, G4int& maxTrials, G4int& numMaxT ) ; 
       //  Return number of count/trials, calls, max & no-max

    void PrintStatistics(); 

  private:

    G4int fTotalNoTrials = 0;    //  Counts sum of trials 
    G4int fNumberCalls = 0;      //  Total # of calls to accumulate
    G4int fmaxTrials = 0;        // Max value of trials
    G4int fNoTimesMaxTrials = 0; // How many times maximum is reached

    G4String fName;        //  Identifies stats, and is printed 
    G4String fDescription; //  Explanation of stats

    G4bool fStatsVerbose = false;
      // If verbose and not printed, print on destruction
    G4bool fPrinted = false;
      // Flag, to avoid reprinting on destruction
}; 

#include "G4TrialsCounter.icc"

#endif  /* End of ifndef G4TRIALS_COUNTER_HH */
