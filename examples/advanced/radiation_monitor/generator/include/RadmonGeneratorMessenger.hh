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
// File name:     RadmonGeneratorLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorMessenger.hh,v 1.2 2006-06-28 13:53:03 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Messenger to change sources properties
//

#ifndef   RADMONGENERATORMESSENGER_HH
 #define  RADMONGENERATORMESSENGER_HH
 
 // Include files
 #include "RadmonMessenger.hh"
 
 // Forward declarations
 class RadmonVGeneratorLayout;
 
 class RadmonGeneratorMessenger : public RadmonMessenger
 {
  public:
                                                RadmonGeneratorMessenger(RadmonVGeneratorLayout * layout);
                                               ~RadmonGeneratorMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorMessenger(void);
                                                RadmonGeneratorMessenger(const RadmonGeneratorMessenger & copy);
   RadmonGeneratorMessenger &                   operator=(const RadmonGeneratorMessenger & copy);

  // Private Attributes
   RadmonVGeneratorLayout *                     generatorLayout;

  // Commands
   RADMON_DECLARE_COMMAND(InsertSource);
   RADMON_DECLARE_COMMAND(SetRelativeSourceIntensity);
   RADMON_DECLARE_COMMAND(RemoveSource);

   RADMON_DECLARE_COMMAND(AppendSourceAlgorithm);
   RADMON_DECLARE_COMMAND(SetSourceAlgorithmType);
   RADMON_DECLARE_COMMAND(RemoveSourceAlgorithm);

   RADMON_DECLARE_COMMAND(SetSourceAlgorithmAttribute);
   RADMON_DECLARE_COMMAND(ClearSourceAlgorithmAttribute);

   RADMON_DECLARE_COMMAND(Load);
   RADMON_DECLARE_COMMAND(Save);
   RADMON_DECLARE_COMMAND(DumpLayout);
 };
#endif /* RADMONGENERATORMESSENGER_HH */
