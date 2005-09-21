//
// File name:     RadmonTLabelledCollection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTLabelledCollection.hh,v 1.2 2005-09-21 14:55:28 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Container for objects with the following methods
//                 - Constructor
//                 - Copy-constructor
//                 - Destructor
//                 - Copy operation
//                 - const G4String & GetLabel(void) const
//

#ifndef   RADMONTLABELLEDCOLLECTION_HH
 #define  RADMONTLABELLEDCOLLECTION_HH
 
 // Include files
 #include "globals.hh"

 #include <vector>
 
 // Forward declarations
 class G4String;

 template <typename LabelledObject>
 class RadmonTLabelledCollection
 {
  public:
                                                RadmonTLabelledCollection();
                                                RadmonTLabelledCollection(const RadmonTLabelledCollection<LabelledObject> & copy);
                                               ~RadmonTLabelledCollection();

   RadmonTLabelledCollection<LabelledObject> &  operator=(const RadmonTLabelledCollection<LabelledObject> & copy);

   G4int                                        GetNItems() const;
   G4bool                                       Empty() const;

   G4bool                                       ExistsItemByLabel(const G4String & label) const;
   G4int                                        MultiplicityItemByLabel(const G4String & label) const;

   const LabelledObject &                       GetItem(G4int index) const;
   LabelledObject &                             GetItem(G4int index);
   const LabelledObject &                       FindItemByLabel(const G4String & label, G4int count = 0) const;
   LabelledObject &                             FindItemByLabel(const G4String & label, G4int count = 0);

   LabelledObject &                             InsertItemAfter(G4int index);
   LabelledObject &                             InsertItemBefore(G4int index);
   LabelledObject &                             AppendItem(void);
   LabelledObject &                             PrependItem(void);

   void                                         RemoveItemByLabel(const G4String & label, G4int count = 0);
   void                                         RemoveItemsByLabel(const G4String & label);
   void                                         RemoveItem(G4int index);
   void                                         RemoveItemsByRange(G4int first, G4int last);
   void                                         RemoveAllItems();

  private:
   typedef std::vector<LabelledObject>          Collection;
   Collection                                   collection;
 };

 // Inline implementations 
 #include "RadmonTLabelledCollection.icc"
#endif /* RADMONTLABELLEDCOLLECTION_HH */
