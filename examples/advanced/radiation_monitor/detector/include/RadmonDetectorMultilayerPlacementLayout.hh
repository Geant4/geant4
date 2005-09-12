//
// File name:     RadmonDetectorMultilayerPlacementLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementLayout.hh,v 1.2 2005-09-12 17:14:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class to manage placed multilayer
//

#ifndef   RADMONDETECTORMULTILAYERPLACEMENTLAYOUT_HH
 #define  RADMONDETECTORMULTILAYERPLACEMENTLAYOUT_HH
 
 // Include files
 #include "G4RotationMatrix.hh"
 #include "G4ThreeVector.hh"
 #include "G4String.hh"
 
 class RadmonDetectorMultilayerPlacementLayout
 {
  public:
   inline                                       RadmonDetectorMultilayerPlacementLayout();
   inline                                       RadmonDetectorMultilayerPlacementLayout(const RadmonDetectorMultilayerPlacementLayout & copy);
   inline                                      ~RadmonDetectorMultilayerPlacementLayout();

   inline RadmonDetectorMultilayerPlacementLayout & operator=(const RadmonDetectorMultilayerPlacementLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline const G4String &                      GetMultilayerLabel(void) const;

   inline void                                  SetLabel(const G4String & label);
   inline void                                  SetMultilayerLabel(const G4String & label);

   inline const G4ThreeVector &                 GetAbsolutePosition(void) const;
   inline const G4RotationMatrix &              GetAbsoluteRotation(void) const;
   G4ThreeVector                                GetRelativePosition(const RadmonDetectorMultilayerPlacementLayout & reference) const;
   G4RotationMatrix                             GetRelativeRotation(const RadmonDetectorMultilayerPlacementLayout & reference) const;

   inline void                                  SetAbsolutePosition(const G4ThreeVector & position);
   inline void                                  SetAbsoluteRotation(const G4RotationMatrix & rotation);
   void                                         SetRelativePosition(const RadmonDetectorMultilayerPlacementLayout & reference, const G4ThreeVector & position);
   void                                         SetRelativeRotation(const RadmonDetectorMultilayerPlacementLayout & reference, const G4RotationMatrix & rotation);

   void                                         DumpLayout(std::ostream & out, const G4String & indent=G4String()) const;

  private:
  // Private attributes
   G4String                                     placementLabel;
   G4String                                     multilayerLabel;
   G4RotationMatrix                             multilayerRotation;
   G4ThreeVector                                multilayerPosition;
 };

 // Inline implementations
 #include "RadmonDetectorMultilayerPlacementLayout.icc"
#endif /* RADMONDETECTORMULTILAYERPLACEMENTLAYOUT_HH */
