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
/// \file SAXSSteppingAction.hh
/// \brief Definition of the SAXSSteppingAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SAXSSteppingAction_h
#define SAXSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <map>

class SAXSEventAction;
class G4LogicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Stepping Action.
/// It is used to retrive information about every interaction of the photons
/// occuring iside the phantom, and in particular, scattering events. 
/// The total number of scattering events (subdivided by tyoe) is passed to 
/// the EventAction class and stored at the end of the event.

class SAXSSteppingAction : public G4UserSteppingAction
{
public:
    SAXSSteppingAction(SAXSEventAction*);
    virtual ~SAXSSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

private:
    SAXSEventAction* fEventAction;   

        G4LogicalVolume* fPhantom;
    
    G4int fEventNumber;
    G4int fNSe;   
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

