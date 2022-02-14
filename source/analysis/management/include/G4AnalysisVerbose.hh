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

// Utility class for analysis category messages.

// Author: Ivana Hrivnacova, 17/10/2011  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisVerbose_h
#define G4AnalysisVerbose_h 1

#include "globals.hh"

#include <array>

class G4AnalysisVerbose
{
  public:
    G4AnalysisVerbose();
    ~G4AnalysisVerbose() = default;

    void Message(G4int verboseLevel,
                 const G4String& action,
                 const G4String& object,
                 const G4String& objectName,
                 G4bool success = true) const;

  private:
    // Static data members
    static constexpr int fkMaxLevel = 4;

    // Data members
    std::array<G4String, fkMaxLevel> fDoneText 
      { G4String("- done"), G4String("- done"), G4String(), G4String() };
    std::array<G4String, fkMaxLevel> fToBeDoneText 
      { G4String(), G4String(), G4String("done "), G4String("going to ") };
    G4String fFailureText
      { "has failed" };
};

#endif
