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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PionMinusInelasticProcess.hh,v 1.4 2001-07-11 10:08:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // PionMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 05-Nov-1996
 // Last modified: 03-Apr-1997

#ifndef G4PionMinusInelasticProcess_h
#define G4PionMinusInelasticProcess_h 1
 
// Class Description
// Process for PionMinus Inelastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End

//#include "G4HadronicInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
 
// class G4PionMinusInelasticProcess : public G4HadronicInelasticProcess
 class G4PionMinusInelasticProcess : public G4HadronInelasticProcess
 {
    
 public:
    
    G4PionMinusInelasticProcess(
     const G4String& processName = "PionMinusInelastic" ) :
      //      G4HadronicInelasticProcess( processName, G4PionMinus::PionMinus() )
      G4HadronInelasticProcess( processName, G4PionMinus::PionMinus() )
    { }
    
    ~G4PionMinusInelasticProcess()
    { }
    
 };
 
#endif
