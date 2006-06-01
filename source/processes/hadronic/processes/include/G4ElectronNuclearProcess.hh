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
// $Id: G4ElectronNuclearProcess.hh,v 1.4 2006-06-01 15:32:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class Description:
//
// Process for electron nuclear inelastic scattering; 
// to be used in your physics list in case you need this physics.
//

// Hadronic Process: Ion Inelastic Process
// J.P. Wellisch, CERN, Apr. 14 2000
// Last modified: 03-Apr-1997

#ifndef G4ElectronNuclearProcess_h
#define G4ElectronNuclearProcess_h 1
 
#include "G4HadronInelasticProcess.hh"
#include "G4ElectroNuclearCrossSection.hh"
 

 class G4ElectronNuclearProcess : public G4HadronInelasticProcess
 {
 public:
    
    G4ElectronNuclearProcess( const G4String& processName = "ElectroNuclear" );
    ~G4ElectronNuclearProcess();

 private:
 
   G4ElectroNuclearCrossSection theData;
 };

#endif

