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
// File name:     RadmonMaterialsMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMaterialsMessenger.hh,v 1.5 2006/06/29 16:16:49 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONMATERIALSMESSENGER_HH
 #define  RADMONMATERIALSMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"

 // Forward declarations
 class RadmonMaterialsManager;
 
 class RadmonMaterialsMessenger : public RadmonMessenger
 {
  public:
                                                RadmonMaterialsMessenger(RadmonMaterialsManager * manager);
                                               ~RadmonMaterialsMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonMaterialsMessenger();
                                                RadmonMaterialsMessenger(const RadmonMaterialsMessenger & copy);
   RadmonMaterialsMessenger &                   operator=(const RadmonMaterialsMessenger & copy);

  // Private Attributes
   RadmonMaterialsManager *                     materialsManager;

  // Commands
   RADMON_DECLARE_COMMAND(CreateElement);
   RADMON_DECLARE_COMMAND(CreateMaterial);
   RADMON_DECLARE_COMMAND(AddComponentByAtoms);
   RADMON_DECLARE_COMMAND(AddComponentByFraction);
   RADMON_DECLARE_COMMAND(SetMaterialColor);
   RADMON_DECLARE_COMMAND(SetMaterialTrasparency);
   RADMON_DECLARE_COMMAND(SetMaterialVisibility);
   RADMON_DECLARE_COMMAND(SetMaterialStyle);
   RADMON_DECLARE_COMMAND(Dump);
   RADMON_DECLARE_COMMAND(Insert);
   RADMON_DECLARE_COMMAND(Save);
 };
#endif /* RADMONMATERIALSMESSENGER_HH */
