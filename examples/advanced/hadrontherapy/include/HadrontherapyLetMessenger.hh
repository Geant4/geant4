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
// $Id: HadrontherapyLetMessenger.hh; May 2007
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), M. Sallemi, A. Salvia
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifndef HadrontherapyLetMessenger_h
#define HadrontherapyLetMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class HadrontherapyLet;
class G4UIdirectory;
class G4UIcmdWithAString;
//class G4UIcmdWithADouble;

class HadrontherapyLetMessenger: public G4UImessenger
{
  public:
  HadrontherapyLetMessenger(HadrontherapyLet* );
  ~HadrontherapyLetMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:

  // Pointer to the let component
  HadrontherapyLet* hadrontherapyletp;

  G4UIdirectory* letDir;
  G4UIcmdWithAString* letDepthCmd;

};
#endif
