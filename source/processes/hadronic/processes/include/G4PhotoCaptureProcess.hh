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
// $Id: G4PhotoCaptureProcess.hh 69176 2013-04-22 06:51:19Z gcosmo $
//
// Class Description
// Process for photon nuclear inelastic scattering; 
// to be used in your physics list in case you need this physics.
//

// Hadronic Process: Ion Inelastic Process
// J.P. Wellisch, CERN, Apr. 14 2000
// Last modified: 03-Apr-1997

#ifndef G4PhotoCaptureProcess_h
#define G4PhotoCaptureProcess_h 1
 
#include "G4HadronicProcess.hh"
//#include "G4PhotoCaptureCrossSection.hh"


class G4PhotoCaptureProcess : public G4HadronicProcess
{
  public:
    
    G4PhotoCaptureProcess(const G4String& processName = "photonCapture")
      : G4HadronicProcess( processName,fCapture )
    {
    }

  ~G4PhotoCaptureProcess() {}

  virtual void ProcessDescription(std::ostream& outFile) const;
};

#endif

