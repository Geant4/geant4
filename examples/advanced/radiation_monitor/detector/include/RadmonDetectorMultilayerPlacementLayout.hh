//
// File name:     RadmonDetectorMultilayerPlacementLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementLayout.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
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
                                                RadmonDetectorMultilayerPlacementLayout();
                                                RadmonDetectorMultilayerPlacementLayout(const RadmonDetectorMultilayerPlacementLayout & copy);
                                               ~RadmonDetectorMultilayerPlacementLayout();

   RadmonDetectorMultilayerPlacementLayout &    operator=(const RadmonDetectorMultilayerPlacementLayout & copy);

   const G4String &                             GetLabel(void) const;
   const G4String &                             GetMultilayerLabel(void) const;

   void                                         SetLabel(const G4String & label);
   void                                         SetMultilayerLabel(const G4String & label);

   const G4ThreeVector &                        GetAbsolutePosition(void) const;
   const G4RotationMatrix &                     GetAbsoluteRotation(void) const;
   G4ThreeVector                                GetRelativePosition(const RadmonDetectorMultilayerPlacementLayout & reference) const;
   G4RotationMatrix                             GetRelativeRotation(const RadmonDetectorMultilayerPlacementLayout & reference) const;

   void                                         SetAbsolutePosition(const G4ThreeVector & position);
   void                                         SetAbsoluteRotation(const G4RotationMatrix & rotation);
   void                                         SetRelativePosition(const RadmonDetectorMultilayerPlacementLayout & reference, const G4ThreeVector & position);
   void                                         SetRelativeRotation(const RadmonDetectorMultilayerPlacementLayout & reference, const G4RotationMatrix & rotation);

   void                                         DumpLayout(std::ostream & out) const;

  private:
  // Private attributes
   G4String                                     placementLabel;
   G4String                                     multilayerLabel;
   G4RotationMatrix                             multilayerRotation;
   G4ThreeVector                                multilayerPosition;
 };
#endif /* RADMONDETECTORMULTILAYERPLACEMENTLAYOUT_HH */
