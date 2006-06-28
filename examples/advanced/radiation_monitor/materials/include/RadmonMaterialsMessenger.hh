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
// File name:     RadmonMaterialsMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMaterialsMessenger.hh,v 1.4 2006-06-28 13:54:15 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
