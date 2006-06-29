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
//
 // G4 Process: KaonZeroL Inelastic Process
 // J.L. Chuma, TRIUMF, 11-Feb-1997
 // Last modified: 03-Apr-1997

#ifndef G4KaonZeroLInelasticProcess_h
#define G4KaonZeroLInelasticProcess_h 1

// Class Description
// Process for KaonZeroLong Inelastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End

//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4KaonZeroLInelasticProcess : public G4HadronicInelasticProcess
 class G4KaonZeroLInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4KaonZeroLInelasticProcess(
     const G4String& processName = "KaonZeroLInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4KaonZeroLong::KaonZeroLong() )
      G4HadronInelasticProcess( processName, G4KaonZeroLong::KaonZeroLong() )
    { }
    
    ~G4KaonZeroLInelasticProcess()
    { }
    
 };
 
#endif
