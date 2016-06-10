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
// $Id: G4HnMessenger.hh 66310 2012-12-17 11:56:35Z ihrivnac $

// The messenger class for histogram information management.
// It implements commands in /analysis/h1 directory.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4HnMessenger_h
#define G4HnMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>

class G4HnManager;
class G4AnalysisMessengerHelper;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;

class G4HnMessenger : public G4UImessenger
{
  public:
    explicit G4HnMessenger(G4HnManager& manager);
    virtual ~G4HnMessenger();
   
    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String value) final;
    
  private:
    void SetHnAsciiCmd();
    void SetHnActivationCmd();
    void SetHnActivationToAllCmd();
    void SetHnPlottingCmd();
    void SetHnPlottingToAllCmd();
 
    G4HnManager&  fManager; ///< Associated class
    std::unique_ptr<G4AnalysisMessengerHelper>   fHelper;

    std::unique_ptr<G4UIcmdWithAnInteger>  fSetHnAsciiCmd;   
    std::unique_ptr<G4UIcommand>       fSetHnActivationCmd;   
    std::unique_ptr<G4UIcmdWithABool>  fSetHnActivationAllCmd;   
    std::unique_ptr<G4UIcommand>       fSetHnPlottingCmd;   
    std::unique_ptr<G4UIcmdWithABool>  fSetHnPlottingAllCmd;   
};
  
#endif

