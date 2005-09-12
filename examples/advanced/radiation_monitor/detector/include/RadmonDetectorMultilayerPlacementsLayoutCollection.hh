//
// File name:     RadmonDetectorMultilayerPlacementsLayoutCollection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementsLayoutCollection.hh,v 1.2 2005-09-12 17:14:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class to collect placed multilayer
//

#ifndef   RADMONDETECTORMULTILAYERPLACEMENTSLAYOUTCOLLECTION_HH
 #define  RADMONDETECTORMULTILAYERPLACEMENTSLAYOUTCOLLECTION_HH

 // Include files
 #include "RadmonTLabelledCollection.hh"
 
 #include "globals.hh"
 
 // Forward declarations
 class G4String;
 class RadmonDetectorMultilayerPlacementLayout;

 class RadmonDetectorMultilayerPlacementsLayoutCollection
 {
  public:
   inline                                       RadmonDetectorMultilayerPlacementsLayoutCollection();
   inline                                      ~RadmonDetectorMultilayerPlacementsLayoutCollection();

   G4int                                        GetNPlacements(void) const;
   G4bool                                       Empty(void) const;

   const RadmonDetectorMultilayerPlacementLayout & GetPlacement(G4int index) const;
   RadmonDetectorMultilayerPlacementLayout &    GetPlacement(G4int index);

   G4bool                                       ExistsPlacementByLabel(const G4String & label) const;
   G4int                                        MultiplicityPlacementByLabel(const G4String & label) const;

   const RadmonDetectorMultilayerPlacementLayout & FindPlacementByLabel(const G4String & label, G4int count=0) const;
   RadmonDetectorMultilayerPlacementLayout &    FindPlacementByLabel(const G4String & label, G4int count=0);

   RadmonDetectorMultilayerPlacementLayout &    CreatePlacement(void);

   void                                         RemovePlacementByLabel(const G4String & label, G4int count=0);
   void                                         RemovePlacementsByLabel(const G4String & label);
   void                                         RemovePlacement(G4int index);
   void                                         RemoveAllPlacements(void);
 
   void                                         DumpLayout(std::ostream & out, const G4String & indent) const;

  private:
  // Hidden constructors and operators
                                                RadmonDetectorMultilayerPlacementsLayoutCollection(const RadmonDetectorMultilayerPlacementsLayoutCollection & copy);
   RadmonDetectorMultilayerPlacementsLayoutCollection & operator=(const RadmonDetectorMultilayerPlacementsLayoutCollection & copy);

  // Private attributes
   RadmonTLabelledCollection<RadmonDetectorMultilayerPlacementLayout> multilayerPlacementsCollection;
 };

 // Inline implementations
 #include "RadmonDetectorMultilayerPlacementsLayoutCollection.icc"
#endif /* RADMONDETECTORMULTILAYERPLACEMENTSLAYOUTCOLLECTION_HH */
