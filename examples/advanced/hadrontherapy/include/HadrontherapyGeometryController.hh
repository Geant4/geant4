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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy
#ifndef HadrontherapyGeometryController_hh
#define HadrontherapyGeometryController_hh 1

#include "globals.hh"
#include "G4String.hh"
#include "G4VUserDetectorConstruction.hh"

/**
 * Controller for geometry selection
 *
 * This controller is called by the geometry messenger and used to
 * select the geometry. Each available geometry must have unique name
 * and it must be known by the geometry controller.
 */
class HadrontherapyGeometryController
{
public:
  HadrontherapyGeometryController();
  ~HadrontherapyGeometryController();

  /**
   * Select a geometry by name.
   */
  void SetGeometry(G4String);

private:
  void registerGeometry(G4VUserDetectorConstruction *detector);
};

#endif
