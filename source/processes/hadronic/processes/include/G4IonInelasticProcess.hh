//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
 // Hadronic Process: Ion Inelastic Process
 // J.P. Wellisch, CERN, Apr. 14 2000
 // Last modified: 03-Apr-1997

#ifndef G4IonInelasticProcess_h
#define G4IonInelasticProcess_h 1
 
// Class Description
// Process for Ion Inelastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End

#include "G4HadronInelasticProcess.hh"
#include "G4GenericIon.hh"
 

 class G4IonInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4IonInelasticProcess(
     const G4String& processName = "IonInelastic" ) :

      G4HadronInelasticProcess( processName, G4GenericIon::GenericIon() )
    { }
        
    ~G4IonInelasticProcess()
    { }
 };

#endif

