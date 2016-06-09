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
// File name:     RadmonSubPhysicsListWithLabelMessenger.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelMessenger.hh,v 1.4 2006/06/29 16:18:21 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   UI commands for dumping the available constructors
//

#ifndef   RADMONSUBPHYSICSLISTWITHLABELMESSENGER_HH
 #define  RADMONSUBPHYSICSLISTWITHLABELMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"
 #include <set>

 class RadmonSubPhysicsListWithLabelMessenger : public RadmonMessenger
 {
  public:
   static RadmonSubPhysicsListWithLabelMessenger * Instance(void);
  
   void                                         AddAvailablePhysicsList(const G4String & name);
   void                                         RemoveAvailablePhysicsList(const G4String & name);
  
   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonSubPhysicsListWithLabelMessenger();
                                                RadmonSubPhysicsListWithLabelMessenger(const RadmonSubPhysicsListWithLabelMessenger & copy);
                                               ~RadmonSubPhysicsListWithLabelMessenger();
   RadmonSubPhysicsListWithLabelMessenger &     operator=(const RadmonSubPhysicsListWithLabelMessenger & copy);

  // Private Data Types
   typedef std::set<G4String>                   AvailablePhysicsLists;
   
  // Private variables
   AvailablePhysicsLists                        availablePhysicsLists;
   
   static RadmonSubPhysicsListWithLabelMessenger * instance;
   
  // Commands
   RADMON_DECLARE_COMMAND(Dump);
 };
#endif /* RADMONSUBPHYSICSLISTWITHLABELMESSENGER_HH */
