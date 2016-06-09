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
// File name:     RadmonDetectorLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayout.hh,v 1.4.2.2 2006/06/29 16:10:37 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//
// Description:   Class to keep track of the configured detector layout
//

#ifndef   RADMONDETECTORLAYOUT_HH
 #define  RADMONDETECTORLAYOUT_HH
 
 // Include files
 #include "RadmonVDetectorLayout.hh"

 #include "RadmonDetectorEnvironmentLayout.hh"
 #include "RadmonDetectorMultilayersLayoutCollection.hh"
 #include "RadmonDetectorMultilayerPlacementsLayoutCollection.hh"
 
 class RadmonDetectorLayout : public RadmonVDetectorLayout
 {
  public:
                                                RadmonDetectorLayout();
                                               ~RadmonDetectorLayout();

   virtual void                                 EnableEnvironment(void);
   virtual void                                 DisableEnvironment(void);
   virtual G4bool                               IsEnabledEnvironment(void) const;

   virtual void                                 SetEnvironmentType(const G4String & type);
   virtual const G4String &                     GetEnvironmentType() const;
   virtual G4int                                GetEnvironmentNAttributes(void) const;
   virtual const G4String &                     GetEnvironmentAttributeName(G4int index) const;
   virtual void                                 SetEnvironmentAttribute(const G4String & attributeName, const G4String & attributeValue);
   virtual const G4String                       GetEnvironmentAttribute(const G4String & attributeName, const G4String & defaultAttributeValue) const;
   virtual void                                 ClearEnvironmentAttribute(const G4String & attributeName);

   virtual void                                 CreateMultilayer(const G4String & multilayerLabel);
   virtual void                                 RemoveMultilayer(const G4String & multilayerLabel);
   virtual void                                 SetMultilayerWidth(const G4String & multilayerLabel, G4double width);
   virtual G4double                             GetMultilayerWidth(const G4String & multilayerLabel) const;
   virtual void                                 SetMultilayerHeight(const G4String & multilayerLabel, G4double height);
   virtual G4double                             GetMultilayerHeight(const G4String & multilayerLabel) const;
   virtual G4double                             GetMultilayerTotalThickness(const G4String & multilayerLabel) const;

   virtual void                                 AppendLayerToMultilayer(const G4String & multilayerLabel, const G4String & layerLabel);
   virtual void                                 RemoveLayerFromMultilayer(const G4String & multilayerLabel, const G4String & layerLabel);
   virtual void                                 RemoveAllLayersFromMultilayer(const G4String & multilayerLabel);
   virtual G4int                                GetMultilayerNLayers(const G4String & multilayerLabel) const;
   virtual const G4String &                     GetMultilayerLayerLabel(const G4String & multilayerLabel, G4int index) const;

   virtual void                                 SetLayerThickness(const G4String & multilayerLabel, const G4String & layerLabel, G4double thickness);
   virtual G4double                             GetLayerThickness(const G4String & multilayerLabel, const G4String & layerLabel) const;
   virtual void                                 SetLayerType(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & type);
   virtual const G4String &                     GetLayerType(const G4String & multilayerLabel, const G4String & layerLabel) const;

   virtual G4int                                GetLayerNAttributes(const G4String & multilayerLabel, const G4String & layerLabel) const;
   virtual const G4String &                     GetLayerAttributeName(const G4String & multilayerLabel, const G4String & layerLabel, G4int index) const;
   virtual void                                 SetLayerAttribute(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & attributeName, const G4String & attributeValue);
   virtual const G4String                       GetLayerAttribute(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & attributeName, const G4String & defaultAttributeValue) const;
   virtual void                                 ClearLayerAttribute(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & attributeName);

   virtual void                                 CreatePlacement(const G4String & placementLabel, const G4String & multilayerName);
   virtual G4int                                GetNPlacements() const;
   virtual const G4String &                     GetPlacementLabel(G4int index) const;
   virtual void                                 RemovePlacement(const G4String & placementLabel);

   virtual const G4String &                     GetPlacementMultilayerType(const G4String & placementLabel) const;
   virtual void                                 SetPlacementPosition(const G4String & placementLabel, const G4ThreeVector & position);
   virtual const G4ThreeVector &                GetPlacementPosition(const G4String & placementLabel) const;
   virtual void                                 SetPlacementPosition(const G4String & placementLabel, const G4String & originLabel, const G4ThreeVector & offset);
   virtual void                                 SetPlacementRotation(const G4String & placementLabel, const G4RotationMatrix & rotation);
   virtual const G4RotationMatrix &             GetPlacementRotation(const G4String & placementLabel) const;
   virtual void                                 SetPlacementRotation(const G4String & placementLabel, const G4String & originLabel, const G4RotationMatrix & relativeRotation);

   virtual void                                 DumpLayout(std::ostream & out) const;

   virtual G4bool                               Load(std::istream & in);
   virtual G4bool                               Save(std::ostream & out) const;


  private:
   inline bool                                  IsPlaced(const G4String & multilayerLabel);
   inline RadmonDetectorMultilayerLayout *      FindMultilayer(const G4String & multilayerLabel);
   inline const RadmonDetectorMultilayerLayout *FindMultilayer(const G4String & multilayerLabel) const;
   inline RadmonDetectorMultilayerPlacementLayout * FindPlacement(const G4String & placementLabel);
   inline const RadmonDetectorMultilayerPlacementLayout *FindPlacement(const G4String & placementLabel) const;
   inline RadmonDetectorLayerLayout *           FindLayer(const G4String & multilayerLabel, const G4String & layerLabel);
   inline const RadmonDetectorLayerLayout *     FindLayer(const G4String & multilayerLabel, const G4String & layerLabel) const;

   inline G4String &                            GetNullStr() const;
   inline G4ThreeVector &                       GetNullPosition() const;
   inline G4RotationMatrix &                    GetNullRotationMatrix() const;

  // Hidden constructors and operators
                                                RadmonDetectorLayout(const RadmonDetectorLayout & copy);
   RadmonDetectorLayout &                       operator=(const RadmonDetectorLayout & copy);

  // Private attributes
   RadmonDetectorMultilayerPlacementsLayoutCollection multilayerPlacementsCollection;
   RadmonDetectorMultilayersLayoutCollection    multilayersCollection;
   RadmonDetectorEnvironmentLayout              environment;
 };
#endif /* RADMONDETECTORLAYOUT_HH */
