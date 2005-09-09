//
// File name:     RadmonDetectorLayerLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerLayout.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class that describes layer informations
//

#ifndef   RADMONDETECTORLAYERLAYOUT_HH
 #define  RADMONDETECTORLAYERLAYOUT_HH

 // Include files
 #include "RadmonDetectorLayoutEntityWithAttributes.hh"

 class RadmonDetectorLayerLayout : public RadmonDetectorLayoutEntityWithAttributes
 {
  public:
   inline                                       RadmonDetectorLayerLayout();
   inline                                       RadmonDetectorLayerLayout(const RadmonDetectorLayerLayout & copy);
   inline                                      ~RadmonDetectorLayerLayout();

   inline RadmonDetectorLayerLayout &           operator=(const RadmonDetectorLayerLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline const G4String &                      GetType(void) const;
   inline G4double                              GetThickness(void) const;
   inline void                                  SetLabel(const G4String & label);
   inline void                                  SetType(const G4String & type);
   inline void                                  SetThickness(G4double thickness);

   inline void                                  DumpLayout(std::ostream & out) const;

  private:
  // Private attributes
   G4String                                     layerLabel;
   G4String                                     layerType;
   G4double                                     layerThickness;
 };

 // Inline implementations
 #include "RadmonDetectorLayerLayout.icc"
#endif /* RADMONDETECTORLAYERLAYOUT_HH */

