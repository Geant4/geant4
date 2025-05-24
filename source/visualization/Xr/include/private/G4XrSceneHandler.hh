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
//
//
// John Allison  5th April 2001
// A template for a simplest possible graphics driver.
//?? Lines or sections marked like this require specialisation for your driver.

#ifndef G4XRSCENEHANDLER_HH
#define G4XRSCENEHANDLER_HH

#include "G4VSceneHandler.hh"
#include "G4XrViewer.hh"

#include <map>
#include <vector>

class G4XrSceneHandler : public G4VSceneHandler
{
  public:
    G4XrSceneHandler(G4VGraphicsSystem& system, const G4String& name);
    ~G4XrSceneHandler() override = default;

    ////////////////////////////////////////////////////////////////
    // Required implementation of pure virtual functions...

    using G4VSceneHandler::AddPrimitive;
    void AddPrimitive(const G4Polyline&) override;
    void AddPrimitive(const G4Text&) override;
    void AddPrimitive(const G4Circle&) override;
    void AddPrimitive(const G4Square&) override;
    void AddPrimitive(const G4Polyhedron&) override;

  protected:
    static G4int fSceneIdCount;  // Counter for Vtk scene handlers.

  private:
    friend class G4XrViewer;
};

#endif
