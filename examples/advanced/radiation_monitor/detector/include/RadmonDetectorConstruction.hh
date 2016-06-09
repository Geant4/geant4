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
// File name:     RadmonDetectorConstruction.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorConstruction.hh,v 1.3.2.2 2006/06/29 16:09:21 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   Implementation of the G4VUserDetectorConstruction
//

#ifndef   RADMONDETECTORCONSTRUCTION_HH
 #define  RADMONDETECTORCONSTRUCTION_HH

 // Include files
 #include "globals.hh"
 #include "G4VUserDetectorConstruction.hh"
 #include "RadmonVLayoutObserver.hh"
 #include <stack>
 #include <utility>
 
 // Forward declaration
 class G4VPhysicalVolume;
 class G4VSolid;
 class G4LogicalVolume;
 class RadmonVDetectorLayout;
 class RadmonVDetectorEntityConstructor;
 class RadmonVDetectorEntitiesConstructorsFactory;
 class G4String;

 class RadmonDetectorConstruction : public G4VUserDetectorConstruction, public RadmonVLayoutObserver
 {
  public:
                                                RadmonDetectorConstruction(RadmonVDetectorLayout * layout, RadmonVDetectorEntitiesConstructorsFactory * factory);
   virtual                                     ~RadmonDetectorConstruction();

   virtual G4VPhysicalVolume *                  Construct(void);

   virtual void                                 OnLayoutChange(void);

  private:
  // Private methods
   void                                         Destruct(void);
   
   void                                         BuildEnvironmentFromType(const G4String & type);
   void                                         BuildEnvironmentSphere(void);
   void                                         BuildMultilayer(G4int index);

  // Hidden constructors and operators
                                                RadmonDetectorConstruction();
                                                RadmonDetectorConstruction(const RadmonDetectorConstruction & copy);
   RadmonDetectorConstruction &                 operator=(const RadmonDetectorConstruction & copy);

  // Private data types
   typedef std::pair<RadmonVDetectorEntityConstructor *, G4VPhysicalVolume *> LayerItem;
   typedef std::stack<LayerItem>                LayersStack;

  // Private attributes
   RadmonVDetectorLayout *                      detectorLayout;
   RadmonVDetectorEntitiesConstructorsFactory * constructorsFactory;
    
   LayersStack                                  layersStack;

   RadmonVDetectorEntityConstructor *           environmentConstructor;
   G4VPhysicalVolume *                          environmentPhysicalVolume;
   G4LogicalVolume *                            environmentLogicalVolume;
   G4VSolid *                                   environmentSolid;
 };
#endif /* RADMONDETECTORCONSTRUCTION_HH */
