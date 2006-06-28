//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonPhysicsMessenger.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsMessenger.hh,v 1.2 2006-06-28 13:54:57 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
