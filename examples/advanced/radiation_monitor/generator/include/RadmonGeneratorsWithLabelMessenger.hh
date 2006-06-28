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
// File name:     RadmonGeneratorsWithLabelMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelMessenger.hh,v 1.2 2006-06-28 13:53:23 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for dumping the available constructors
//

#ifndef   RADMONGENERATORSWITHLABELMESSENGER_HH
 #define  RADMONGENERATORSWITHLABELMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"
 #include <set>

 class RadmonGeneratorsWithLabelMessenger : public RadmonMessenger
 {
  public:
   static RadmonGeneratorsWithLabelMessenger *  Instance(void);
  
   void                                         AddAvailableGenerator(const G4String & name);
   void                                         RemoveAvailableGenerator(const G4String & name);
  
   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonGeneratorsWithLabelMessenger();
                                                RadmonGeneratorsWithLabelMessenger(const RadmonGeneratorsWithLabelMessenger & copy);
                                               ~RadmonGeneratorsWithLabelMessenger();
   RadmonGeneratorsWithLabelMessenger &         operator=(const RadmonGeneratorsWithLabelMessenger & copy);

  // Private Data Types
   typedef std::set<G4String>                   AvailableGenerators;
   
  // Private variables
   AvailableGenerators                          availableGenerators;
   
   static RadmonGeneratorsWithLabelMessenger * instance;
   
  // Commands
   RADMON_DECLARE_COMMAND(Dump);
 };
#endif /* RADMONGENERATORSWITHLABELMESSENGER_HH */
