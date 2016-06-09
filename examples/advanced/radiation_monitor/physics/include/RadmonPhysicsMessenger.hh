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
// File name:     RadmonPhysicsMessenger.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsMessenger.hh,v 1.3 2006/06/29 16:17:31 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONPHYSICSMESSENGER_HH
 #define  RADMONPHYSICSMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"

 // Forward declarations
 class RadmonVPhysicsLayout;
 
 class RadmonPhysicsMessenger : public RadmonMessenger
 {
  public:
                                                RadmonPhysicsMessenger(RadmonVPhysicsLayout * layout);
                                               ~RadmonPhysicsMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonPhysicsMessenger();
                                                RadmonPhysicsMessenger(const RadmonPhysicsMessenger & copy);
   RadmonPhysicsMessenger &                     operator=(const RadmonPhysicsMessenger & copy);

  // Private Attributes
   RadmonVPhysicsLayout *                       physicsLayout;

  // Commands
   RADMON_DECLARE_COMMAND(AddPhysicsList);
   RADMON_DECLARE_COMMAND(RemovePhysicsList);

   RADMON_DECLARE_COMMAND(SetPhysicsListAttribute);
   RADMON_DECLARE_COMMAND(ClearPhysicsListAttribute);

   RADMON_DECLARE_COMMAND(DumpLayout);
   RADMON_DECLARE_COMMAND(Load);
   RADMON_DECLARE_COMMAND(Save);
 };
#endif /* RADMONPHYSICSMESSENGER_HH */
