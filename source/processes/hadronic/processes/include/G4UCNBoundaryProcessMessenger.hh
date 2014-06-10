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
// $Id: G4UCNBoundaryProcessMessenger 71487 2013-06-17 08:19:40Z gcosmo $
//
//
////////////////////////////////////////////////////////////////////////
// Ultra Cold Neutron (UCN) Boundary Process Messenger Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:         G4UCNBoundaryProcessMessenger
// Description:  Messenger class for UCN material boundaries
// Version:      1.0
// Created:      2014-05-14
// Author:       Peter Gumplinger
// Adopted from: UCNMaterialBoundaryMessenger by Peter Fierlinger 9.9.2004
// Updated:
// mail:         gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4UCNBOUNDARYPROCESSMESSENGER_HH
#define G4UCNBOUNDARYPROCESSMESSENGER_HH 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4UCNBoundaryProcess;

class G4UIdirectory;

class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4UCNBoundaryProcessMessenger : public G4UImessenger
{
  public:

    G4UCNBoundaryProcessMessenger(G4UCNBoundaryProcess* );
    virtual ~G4UCNBoundaryProcessMessenger();

    void SetNewValue(G4UIcommand* ,G4String );

  private:

    G4UCNBoundaryProcess* theUCNBoundaryProcess;

    G4UIdirectory*     boundaryDir;

    G4UIcmdWithAnInteger* VerboseCmd;
    G4UIcmdWithABool* MicroRoughnessCmd;
};

#endif /* G4UCNBOUNDARYPROCESSMESSENGER_HH */
