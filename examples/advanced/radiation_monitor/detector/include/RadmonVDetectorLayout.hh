//
// File name:     RadmonVDetectorLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLayout.hh,v 1.2 2005-09-12 17:14:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class to keep track of the configured detector
//                layout
//

#ifndef   RADMONVDETECTORLAYOUT_HH
 #define  RADMONVDETECTORLAYOUT_HH
 
 // Include files
 #include "RadmonVDetectorLayoutSubject.hh"

 #include "globals.hh"
 #include "G4ThreeVector.hh"
 #include "G4RotationMatrix.hh"
 
 class RadmonVDetectorLayout : public RadmonVDetectorLayoutSubject
 {
  public:
   virtual void                                 EnableEnvironment(void) = 0;
   virtual void                                 DisableEnvironment(void) = 0;
   virtual G4bool                               IsEnabledEnvironment(void) const = 0;

   virtual void                                 SetEnvironmentType(const G4String & type) = 0;
   virtual const G4String &                     GetEnvironmentType() const = 0;
   virtual void                                 SetEnvironmentAttribute(const G4String & attributeName, const G4String & attributeValue) = 0;
   virtual const G4String                       GetEnvironmentAttribute(const G4String & attributeName, const G4String & defaultAttributeValue) const = 0;
   virtual void                                 ClearEnvironmentAttribute(const G4String & attributeName) = 0;

   virtual void                                 CreateMultilayer(const G4String & multilayerLabel) = 0;
   virtual void                                 RemoveMultilayer(const G4String & multilayerLabel) = 0;
   virtual void                                 SetMultilayerWidth(const G4String & multilayerLabel, G4double width) = 0;
   virtual G4double                             GetMultilayerWidth(const G4String & multilayerLabel) const = 0;
   virtual void                                 SetMultilayerHeight(const G4String & multilayerLabel, G4double height) = 0;
   virtual G4double                             GetMultilayerHeight(const G4String & multilayerLabel) const = 0;

   virtual void                                 AppendLayerToMultilayer(const G4String & multilayerLabel, const G4String & layerLabel) = 0;
   virtual void                                 RemoveLayerFromMultilayer(const G4String & multilayerLabel, const G4String & layerLabel) = 0;
   virtual void                                 RemoveAllLayers(const G4String & multilayerLabel) = 0;

   virtual void                                 SetLayerThickness(const G4String & multilayerLabel, const G4String & layerLabel, G4double thickness) = 0;
   virtual G4double                             GetLayerThickness(const G4String & multilayerLabel, const G4String & layerLabel) const = 0;
   virtual void                                 SetLayerType(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & type) = 0;
   virtual const G4String &                     GetLayerType(const G4String & multilayerLabel, const G4String & layerLabel) const = 0;

   virtual void                                 SetLayerAttribute(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & attributeName, const G4String & attributeValue) = 0;
   virtual const G4String                       GetLayerAttribute(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & attributeName, const G4String & defaultAttributeValue) const = 0;
   virtual void                                 ClearLayerAttribute(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & attributeName) = 0;

   virtual void                                 CreatePlacement(const G4String & placementLabel, const G4String & multilayerName) = 0;
   virtual G4int                                GetNPlacements() const = 0;
   virtual const G4String &                     GetPlacementLabel(G4int index) const = 0;
   virtual void                                 RemoveMultilayerPlacement(const G4String & placementLabel) = 0;

   virtual const G4String &                     GetPlacementMultilayerType(const G4String & placementLabel) const = 0;
   virtual void                                 SetPlacementPosition(const G4String & placementLabel, const G4ThreeVector & position) = 0;
   virtual const G4ThreeVector &                GetPlacementPosition(const G4String & placementLabel) const = 0;
   virtual void                                 SetPlacementPosition(const G4String & placementLabel, const G4String & originLabel, const G4ThreeVector & offset) = 0;
   virtual void                                 SetPlacementRotation(const G4String & placementLabel, const G4RotationMatrix & rotation) = 0;
   virtual const G4RotationMatrix &             GetPlacementRotation(const G4String & placementLabel) const = 0;
   virtual void                                 SetPlacementRotation(const G4String & placementLabel, const G4String & originLabel, const G4RotationMatrix & relativeRotation) = 0;

   virtual void                                 DumpLayout(std::ostream & out) const = 0;

   virtual G4bool                               Load(std::istream & in) = 0;
   virtual G4bool                               Save(std::ostream & out) const = 0;

  protected:
   inline                                       RadmonVDetectorLayout();
   inline                                      ~RadmonVDetectorLayout();

  private:
  // Hidden constructors and operators
                                                RadmonVDetectorLayout(const RadmonVDetectorLayout & copy);
    RadmonVDetectorLayout &                     operator=(const RadmonVDetectorLayout & copy);
 };
 
 // Inline implementations
 #include "RadmonVDetectorLayout.icc"
#endif /* RADMONVDETECTORLAYOUT_HH */
