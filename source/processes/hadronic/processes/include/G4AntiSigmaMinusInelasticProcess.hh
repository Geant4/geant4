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
// $Id: G4AntiSigmaMinusInelasticProcess.hh,v 1.5 2001-08-01 17:12:09 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: AntiSigmaMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 18-Feb-1997
 // Last modified: 03-Apr-1997
 
 // Note:  there is no .cc file
 
#ifndef G4AntiSigmaMinusInelasticProcess_h
#define G4AntiSigmaMinusInelasticProcess_h 1
 
// Class Description
// Process for AntiSigmaMinus Inelastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End

//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4AntiSigmaMinusInelasticProcess : public G4HadronicInelasticProcess
 class G4AntiSigmaMinusInelasticProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4AntiSigmaMinusInelasticProcess(
     const G4String& processName = "AntiSigmaMinusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4AntiSigmaMinus::AntiSigmaMinus() )
      G4HadronInelasticProcess( processName, G4AntiSigmaMinus::AntiSigmaMinus() )
    { }
    
    ~G4AntiSigmaMinusInelasticProcess()
    { }
 };
 
#endif
 

