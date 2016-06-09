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
// $Id: G4ITModelManager.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4ITMODELMANAGER_H
#define G4ITMODELMANAGER_H

#include "globals.hh"
#include <map>
#include "G4VITModel.hh"

/**
  * G4ITModelManager chooses which model to use according
  * to the global simulation time.
  */
class G4ITModelManager
{

public:
    G4ITModelManager();
    ~G4ITModelManager();
    void Initialize();
    G4ITModelManager(const G4ITModelManager& other);
    G4ITModelManager& operator=(const G4ITModelManager& rhs);
    void SetModel(G4VITModel* aModel, G4double startingTime);
    G4VITModel* GetModel(const G4double globalTime);

protected :
    typedef std::map<G4double /*startingTime*/, G4VITModel* /*aModel*/>  mapModels ;
    mapModels fModels ;
    G4bool fIsInitialized ;
};

#endif // G4ITMODELMANAGER_H
