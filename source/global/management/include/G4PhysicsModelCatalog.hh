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
// $Id: G4PhysicsModelCatalog.hh 67281 2013-02-13 14:41:55Z gcosmo $
//
// 
// -----------------------------------------------------------------
//
//      ------------------- class G4PhysicsModelCatalog -----------------
//
// Class description:
//

#ifndef G4PhysicsModelCatalog_HH
#define G4PhysicsModelCatalog_HH

#include "globals.hh"
#include <vector>
#include "G4String.hh"

typedef std::vector<G4String> modelCatalog;

class G4PhysicsModelCatalog
{
private:  // with description

    G4PhysicsModelCatalog();
    G4PhysicsModelCatalog(const G4PhysicsModelCatalog&);
    G4PhysicsModelCatalog& operator=(const G4PhysicsModelCatalog&);

public:  // with description
    
    ~G4PhysicsModelCatalog();
    static G4int Register(const G4String&);
    static const G4String& GetModelName(G4int);

public:  // without description

    static G4int GetIndex(const G4String&);
    static G4int Entries();
    static void Destroy();
    
private:

    static modelCatalog* catalog;

};

#endif
