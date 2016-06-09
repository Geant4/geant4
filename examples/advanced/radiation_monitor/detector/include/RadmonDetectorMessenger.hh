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
// File name:     RadmonDetectorMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMessenger.hh,v 1.5 2006/06/29 16:10:39 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONDETECTORMESSENGER_HH
 #define  RADMONDETECTORMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"

 // Forward declarations
 class RadmonVDetectorLayout;
 
 class RadmonDetectorMessenger : public RadmonMessenger
 {
  public:
                                                RadmonDetectorMessenger(RadmonVDetectorLayout * layout);
                                               ~RadmonDetectorMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorMessenger();
                                                RadmonDetectorMessenger(const RadmonDetectorMessenger & copy);
   RadmonDetectorMessenger &                    operator=(const RadmonDetectorMessenger & copy);

  // Private Attributes
   RadmonVDetectorLayout *                      detectorLayout;

  // Commands
   RADMON_DECLARE_COMMAND(EnableEnvironment);
   RADMON_DECLARE_COMMAND(DisableEnvironment);
   RADMON_DECLARE_COMMAND(SetEnvironmentType);
   RADMON_DECLARE_COMMAND(SetEnvironmentAttribute);
   RADMON_DECLARE_COMMAND(ClearEnvironmentAttribute);

   RADMON_DECLARE_COMMAND(CreateMultilayer);
   RADMON_DECLARE_COMMAND(RemoveMultilayer);
   RADMON_DECLARE_COMMAND(SetMultilayerWidth);
   RADMON_DECLARE_COMMAND(SetMultilayerHeight);
   RADMON_DECLARE_COMMAND(AppendLayerToMultilayer);
   RADMON_DECLARE_COMMAND(RemoveLayerFromMultilayer);
   RADMON_DECLARE_COMMAND(RemoveAllLayersFromMultilayer);
   
   RADMON_DECLARE_COMMAND(SetLayerThickness);
   RADMON_DECLARE_COMMAND(SetLayerType);
   RADMON_DECLARE_COMMAND(SetLayerAttribute);
   RADMON_DECLARE_COMMAND(ClearLayerAttribute);

   RADMON_DECLARE_COMMAND(CreatePlacement);
   RADMON_DECLARE_COMMAND(RemovePlacement);
   RADMON_DECLARE_COMMAND(SetPlacementPosition);
   RADMON_DECLARE_COMMAND(SetPlacementRotation);
   RADMON_DECLARE_COMMAND(SetRelativePlacementPosition);
   RADMON_DECLARE_COMMAND(SetRelativePlacementRotation);
   
   RADMON_DECLARE_COMMAND(DumpLayout);
   RADMON_DECLARE_COMMAND(Load);
   RADMON_DECLARE_COMMAND(Save);
 };
#endif /* RADMONDETECTORMESSENGER_HH */
