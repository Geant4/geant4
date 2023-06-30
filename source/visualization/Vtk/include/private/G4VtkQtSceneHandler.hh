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

#ifndef G4VTKQTSCENEHANDLER_HH
#define G4VTKQTSCENEHANDLER_HH

#include "G4VSceneHandler.hh"
#include "G4VtkQtViewer.hh"
#include "G4VtkSceneHandler.hh"

#include <map>
#include <vector>

class G4VtkQtSceneHandler : public G4VtkSceneHandler
{
  public:
    G4VtkQtSceneHandler(G4VGraphicsSystem& system, const G4String& name);
    ~G4VtkQtSceneHandler() override = default;

    ////////////////////////////////////////////////////////////////
    // Implementation of pure virtual functions...
    void AddPrimitive(const G4Polyline&) override;
    void AddPrimitive(const G4Text&) override;
    void AddPrimitive(const G4Circle&) override;
    void AddPrimitive(const G4Square&) override;
    void AddPrimitive(const G4Polyhedron&) override;

  protected:
    static G4int fSceneIdCount;  // Counter for Vtk scene handlers.

  private:
    friend class G4VtkQtViewer;
};

#endif
