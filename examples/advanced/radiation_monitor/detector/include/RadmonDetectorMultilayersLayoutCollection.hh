//
// File name:     RadmonDetectorMultilayersLayoutCollection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayersLayoutCollection.hh,v 1.2 2005-09-12 17:14:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class to collect multilayers
//

#ifndef   RADMONDETECTORMULTILAYERSLAYOUTCOLLECTION_HH
 #define  RADMONDETECTORMULTILAYERSLAYOUTCOLLECTION_HH

 // Include files
 #include "RadmonDetectorMultilayerLayout.hh"
 #include "RadmonTLabelledCollection.hh"

 class RadmonDetectorMultilayersLayoutCollection
 {
  public:
   inline                                       RadmonDetectorMultilayersLayoutCollection();
   inline                                      ~RadmonDetectorMultilayersLayoutCollection();

   G4int                                        GetNMultilayers(void) const;
   G4bool                                       Empty(void) const;

   const RadmonDetectorMultilayerLayout &       GetMultilayer(G4int index) const;
   RadmonDetectorMultilayerLayout &             GetMultilayer(G4int index);

   G4bool                                       ExistsMultilayerByLabel(const G4String & label) const;
   G4int                                        MultiplicityMultilayerByLabel(const G4String & label) const;

   const RadmonDetectorMultilayerLayout &       FindMultilayerByLabel(const G4String & label, G4int count = 0) const;
   RadmonDetectorMultilayerLayout &             FindMultilayerByLabel(const G4String & label, G4int count = 0);

   RadmonDetectorMultilayerLayout &             CreateMultilayer(void);

   void                                         RemoveMultilayerByLabel(const G4String & label, G4int count = 0);
   void                                         RemoveMultilayersByLabel(const G4String & label);
   void                                         RemoveMultilayer(G4int index);
   void                                         RemoveAllMultilayers(void);
 
   void                                         DumpLayout(std::ostream & out, const G4String &indent=G4String()) const;

  private:
  // Hidden constructors and operators
                                                RadmonDetectorMultilayersLayoutCollection(const RadmonDetectorMultilayersLayoutCollection & copy);
   RadmonDetectorMultilayersLayoutCollection &  operator=(const RadmonDetectorMultilayersLayoutCollection & copy);

  // Private attributes
   RadmonTLabelledCollection<RadmonDetectorMultilayerLayout> multilayersCollection;
 };

 // Inline implementations
 #include "RadmonDetectorMultilayersLayoutCollection.icc"
#endif /* RADMONDETECTORMULTILAYERSLAYOUTCOLLECTION_HH */
