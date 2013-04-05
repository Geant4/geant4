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
//  27 Feb 2013: Andrea Dotti, first implementation
//
// MT initializaiton code. Each thread will execute this code that contains the
// per-thread initialization code (e.g. instantiaion of user actions, etc)

#ifndef ExN02WorkerInitialization_hh
#define ExN02WorkerInitialization_hh

#include "G4VUserWorkerInitialization.hh"

class ExN02DetectorConstruction;

class ExN02WorkerInitialization : public G4VUserWorkerInitialization {
protected:
    void WorkerStart() const;
//Application specific stuff
private:
    ExN02DetectorConstruction* detector;
    G4String macroFileName;
public:
    void SetDetectorConstruction( ExN02DetectorConstruction* det ) { detector = det; }
    void SetMacroFileName( const G4String& macFileName ) { macroFileName = macFileName; }
};

#endif //ExN02WorkerInitialization_hh