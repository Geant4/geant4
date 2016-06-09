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
// $Id$

// Author: Ivana Hrivnacova, 17/10/2011  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisVerbose_h
#define G4AnalysisVerbose_h 1

#include "globals.hh"

class G4AnalysisVerbose
{
  public:
    G4AnalysisVerbose(const G4String& type, G4int verboseLevel);
    virtual ~G4AnalysisVerbose();

    void Message(const G4String& action, 
                 const G4String& object, 
                 const G4String& objectName,
                 G4bool success = true);

    void Message(const G4String& action, 
                 const G4String& object, 
                 G4ExceptionDescription& description,
                 G4bool success = true);
                 
    G4String GetType() const;                 
        
  private:
    // data members
    //
    G4String fType;
    G4String fToBeDoneText;
    G4String fDoneText;
    G4String fFailureText;
};

// inline functions

inline  G4String G4AnalysisVerbose::GetType() const {
  return fType;
}  

#endif

