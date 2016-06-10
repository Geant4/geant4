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
// $Id: G4AnalysisVerbose.cc 72292 2013-07-15 12:01:43Z ihrivnac $

// Author: Ivana Hrivnacova, 17/10/2011  (ivana@ipno.in2p3.fr)

#include "G4AnalysisVerbose.hh"
#include "G4UnitsTable.hh"

#include <iostream>

//_____________________________________________________________________________
G4AnalysisVerbose::G4AnalysisVerbose(const G4String& type, G4int verboseLevel)
 : fType(type),
   fToBeDoneText(),
   fDoneText(),
   fFailureText()
{
   if ( verboseLevel == 1 ) fDoneText = "- done";
   if ( verboseLevel == 2 ) fDoneText = "- done";
   if ( verboseLevel == 3 ) fToBeDoneText = "done ";
   if ( verboseLevel == 4 ) fToBeDoneText = "going to ";
   fFailureText = "has failed";
}

//_____________________________________________________________________________
G4AnalysisVerbose::~G4AnalysisVerbose()
{  
}

// 
// public method
//

//_____________________________________________________________________________
void G4AnalysisVerbose::Message(const G4String& action, 
                                const G4String& object, 
                                const G4String& objectName,
                                G4bool success) const
{
  G4cout << "... "
         << fToBeDoneText
         << action
         << " "
         << fType
         << " "
         << object 
         << " : "
         << objectName 
         << " ";

  if ( success )
     G4cout << fDoneText;
  else   
     G4cout << fFailureText;
        
  G4cout << G4endl;
}  
  
//_____________________________________________________________________________
void G4AnalysisVerbose::Message(const G4String& action, 
                                const G4String& object, 
                                G4ExceptionDescription& description,
                                G4bool success) const
{
  G4cout << "... "
         << fToBeDoneText
         << action
         << " "
         << fType
         << " "
         << object 
         << " : "
         << description.str() 
         << " ";

  if ( success )
     G4cout << fDoneText;
  else   
     G4cout << fFailureText;
        
  G4cout << G4endl;
}  
  
