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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRPrimGenActionMessenger.hh
//   Header file of a messenger that handles primary
//   generator action class.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRPrimGenActionMessenger_h
#define GRPrimGenActionMessenger_h 1

class G4UIdirectory;
class G4UIcommand;
class G4UIparameter;

#include "G4UImessenger.hh"

class GRPrimGenActionMessenger: public G4UImessenger
{

  public:
    GRPrimGenActionMessenger(G4int verboseLvl = 0);
    ~GRPrimGenActionMessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  public:
    void UpdateParticleList();

  private:
    G4int verboseLevel;
    G4UIdirectory* srcDir;
    G4UIcommand* defCmd;
    G4UIparameter* particlePara;
};

#endif
