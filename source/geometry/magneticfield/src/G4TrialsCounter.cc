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
// $Id: $
// GEANT4 tag $Name:  $
//
// class G4TrialsCounter
//
// Class inline implementation
//
// Author: Dec 8, 2006  John Apostolakis
// -------------------------------------------------------------------

#include "G4TrialsCounter.hh"
#include "G4ios.hh"

G4TrialsCounter::G4TrialsCounter( const G4String& nameStats,
                                  const G4String& description,
                                        G4bool printOnExit )
  : fName(nameStats), fDescription(description),
    fStatsVerbose(printOnExit), fPrinted(false) 
{ 
  ClearCounts(); 
}
   
G4TrialsCounter::~G4TrialsCounter() 
{ 
  if( (fStatsVerbose) && (!fPrinted) )  { PrintStatistics(); }
}

void
G4TrialsCounter::PrintStatistics()
{
  // Print Statistics
  G4cout << "G4TrialsCounter::PrintStatistics()" << G4endl
         << "Report of counts for " << fDescription  << " : " << G4endl;
  G4cout << "Stats for '" <<  fName << "' > "
         << "  No-trials= " << fTotalNoTrials
         << "  No-calls= "  << fNumberCalls
         << "  Max-trial= " << fmaxTrials
         << "  no-max= "    << fNoTimesMaxTrials 
         << G4endl; 
  fPrinted= true; 
}

void G4TrialsCounter::ClearCounts()
{
  fTotalNoTrials= 0; 
  fNumberCalls  = 0; 
  fmaxTrials    = 0;        // Maximum --> so only unsigned ints expected
  fNoTimesMaxTrials=0; 
}

G4int
G4TrialsCounter::ReturnTotals( G4int& calls, G4int& maxTrials, G4int& numMaxT ) 
{
  calls    = fNumberCalls; 
  maxTrials= fmaxTrials;
  numMaxT  = fNoTimesMaxTrials; 

  return fTotalNoTrials; 
}
