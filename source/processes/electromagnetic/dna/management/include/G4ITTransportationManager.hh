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
// $Id: G4ITTransportationManager.hh 65022 2012-11-12 16:43:12Z gcosmo $
//
/// \brief {Duplicated version of G4TransportationManager.
///         This class just contains the pointer to the navigator object of the
///         simulation.}
//
// -------------------------------------------------------------------

#ifndef G4ITTRANSPORTATIONMANAGER_HH
#define G4ITTRANSPORTATIONMANAGER_HH


#include "globals.hh"

class G4ITNavigator;

class G4ITTransportationManager
{
public:
    static void DeleteInstance();
    static G4ITTransportationManager* GetTransportationManager();
    G4ITNavigator* GetNavigatorForTracking();

private:
    G4ITTransportationManager();
    ~G4ITTransportationManager();
    static G4ITTransportationManager* fpInstance;
    void Initialize();
    G4ITNavigator* fpNavigator;
};

#endif // G4ITTRANSPORTATIONMANAGER_HH
