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
// $Id: G4ConversionFatalError.hh 78955 2014-02-05 09:45:46Z gcosmo $
//
// Fatal G4Exception error
//
// Jane Tinslay, September 2006
//
#ifndef G4CONVERSIONFATALERROR_HH
#define G4CONVERSIONFATALERROR_HH

#include "globals.hh"
#include "G4String.hh"

struct G4ConversionFatalError {
  
  void ReportError(const G4String& input, const G4String& message) const {
    G4ExceptionDescription ed; 
    ed <<input<<": "<<message<<G4endl;
    G4Exception("G4ConversionFatalError::ReportError",
                "greps0101", FatalErrorInArgument, ed);
  }
  
};

#endif
