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
// G4tgbBorderSurface
//
// Class description:
//
// Transient class of a optical surface; 
// builds a G4tgbBorderSurface.
//
#ifndef G4tgbBorderSurface_hh
#define G4tgbBorderSurface_hh 1

#include <vector>
#include <string>

#include "globals.hh"
#include "G4tgrBorderSurface.hh"
#include "G4tgbVolumeMgr.hh"

class G4tgbBorderSurface
{
  public:

    G4tgbBorderSurface(G4tgrBorderSurface*);
    ~G4tgbBorderSurface();

    G4String GetName() const { return theName; }

    // Builds a G4BorderSurface or a G4SkinBorderSurface
    void BuildG4BorderSurface();
    void BuildG4SkinBorderSurface();
    void BuildG4LogicalBorderSurface();

  private:

    G4String theName;
    G4tgrBorderSurface* theTgrBorderSurface = nullptr;
    G4tgbVolumeMgr* fTgbVolmgr = G4tgbVolumeMgr::GetInstance();
    G4OpticalSurface *fOptic;
};

#endif

