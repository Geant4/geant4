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
// File name:     RadmonGeneratorLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorMessenger.hh,v 1.3 2006/06/29 16:15:13 gunter Exp $
// Tag:           $Name: geant4-08-01 $
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
