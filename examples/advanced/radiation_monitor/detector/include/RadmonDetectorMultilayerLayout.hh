//
// File name:     RadmonDetectorMultilayerLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerLayout.hh,v 1.2 2005-09-12 17:14:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class to manage multilayers
//

#ifndef   RADMONDETECTORMULTILAYERLAYOUT_HH
 #define  RADMONDETECTORMULTILAYERLAYOUT_HH

 // Include files
 #include "RadmonDetectorLayerLayout.hh"
 #include "RadmonTLabelledCollection.hh"
 #include "G4String.hh"
 #include "globals.hh"
 
 class RadmonDetectorMultilayerLayout
 {
  public:
   inline                                       RadmonDetectorMultilayerLayout();
                                                RadmonDetectorMultilayerLayout(const RadmonDetectorMultilayerLayout & copy);
   inline                                      ~RadmonDetectorMultilayerLayout();
   
   RadmonDetectorMultilayerLayout &             operator=(const RadmonDetectorMultilayerLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline void                                  SetLabel(const G4String & label);

   inline G4double                              GetWidth(void) const;
   inline G4double                              GetHeight(void) const;
   G4double                                     GetTotalThickness(void) const;
   inline void                                  SetWidth(G4double width);
   inline void                                  SetHeight(G4double height);

   G4int                                        GetNLayers(void) const;
   G4bool                                       Empty(void) const;

   const RadmonDetectorLayerLayout &            GetLayer(G4int index) const;
   RadmonDetectorLayerLayout &                  GetLayer(G4int index);

   G4bool                                       ExistsLayerByLabel(const G4String & layerLabel) const;
   G4int                                        MultiplicityLayerByLabel(const G4String & layerLabel) const;

   const RadmonDetectorLayerLayout &            FindLayerByLabel(const G4String &layerLabel, G4int count = 0) const;
   RadmonDetectorLayerLayout &                  FindLayerByLabel(const G4String & layerLabel, G4int count = 0);

   RadmonDetectorLayerLayout &                  AppendLayer(void);
   RadmonDetectorLayerLayout &                  PrependLayer(void);

   void                                         RemoveLayerByLabel(const G4String & layerLabel, G4int count = 0);
   void                                         RemoveLayersByLabel(const G4String & layerLabel);
   void                                         RemoveLayer(G4int index);
   void                                         RemoveLayersByRange(G4int first, G4int last);
   void                                         RemoveAllLayers(void);

   void                                         DumpLayout(std::ostream & out, const G4String & indent=G4String()) const;

  private:
  // Private attributes
   G4String                                     multilayerLabel;
   G4double                                     multilayerWidth;
   G4double                                     multilayerHeight;
   RadmonTLabelledCollection<RadmonDetectorLayerLayout> multilayerLayersCollection;
 };

 // Inline implementations
 #include "RadmonDetectorMultilayerLayout.icc"
#endif /* RADMONDETECTORMULTILAYERLAYOUT_HH */
