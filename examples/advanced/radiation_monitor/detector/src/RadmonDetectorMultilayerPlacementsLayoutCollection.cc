//
// File name:     RadmonDetectorMultilayerPlacementsLayoutCollection.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementsLayoutCollection.cc,v 1.1 2005-09-12 17:13:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMultilayerPlacementsLayoutCollection.hh"



G4int                                           RadmonDetectorMultilayerPlacementsLayoutCollection :: GetNPlacements(void) const
{
 return multilayerPlacementsCollection.GetNItems();
}



G4bool                                          RadmonDetectorMultilayerPlacementsLayoutCollection :: Empty(void) const
{
 return multilayerPlacementsCollection.Empty();
}





const RadmonDetectorMultilayerPlacementLayout & RadmonDetectorMultilayerPlacementsLayoutCollection :: GetPlacement(G4int index) const
{
 return multilayerPlacementsCollection.GetItem(index);
}



RadmonDetectorMultilayerPlacementLayout &       RadmonDetectorMultilayerPlacementsLayoutCollection :: GetPlacement(G4int index)
{
 return multilayerPlacementsCollection.GetItem(index);
}





G4bool                                          RadmonDetectorMultilayerPlacementsLayoutCollection :: ExistsPlacementByLabel(const G4String & label) const
{
 return multilayerPlacementsCollection.ExistsItemByLabel(label);
}



G4int                                           RadmonDetectorMultilayerPlacementsLayoutCollection :: MultiplicityPlacementByLabel(const G4String & label) const
{
 return multilayerPlacementsCollection.MultiplicityItemByLabel(label);
}





const RadmonDetectorMultilayerPlacementLayout & RadmonDetectorMultilayerPlacementsLayoutCollection :: FindPlacementByLabel(const G4String & label, G4int count) const
{
 return multilayerPlacementsCollection.FindItemByLabel(label, count);
}



RadmonDetectorMultilayerPlacementLayout &       RadmonDetectorMultilayerPlacementsLayoutCollection :: FindPlacementByLabel(const G4String & label, G4int count)
{
 return multilayerPlacementsCollection.FindItemByLabel(label, count);
}





RadmonDetectorMultilayerPlacementLayout &       RadmonDetectorMultilayerPlacementsLayoutCollection :: CreatePlacement(void)
{
 return multilayerPlacementsCollection.AppendItem();
}





void                                            RadmonDetectorMultilayerPlacementsLayoutCollection :: RemovePlacementByLabel(const G4String & label, G4int count)
{
 multilayerPlacementsCollection.RemoveItemByLabel(label, count);
}



void                                            RadmonDetectorMultilayerPlacementsLayoutCollection :: RemovePlacementsByLabel(const G4String & label)
{
 multilayerPlacementsCollection.RemoveItemsByLabel(label);
}



void                                            RadmonDetectorMultilayerPlacementsLayoutCollection :: RemovePlacement(G4int index)
{
 multilayerPlacementsCollection.RemoveItem(index);
}



void                                            RadmonDetectorMultilayerPlacementsLayoutCollection :: RemoveAllPlacements(void)
{
 multilayerPlacementsCollection.RemoveAllItems();
}





void                                            RadmonDetectorMultilayerPlacementsLayoutCollection :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 G4String indent2(indent);
 indent2.prepend("  ");

 const G4int n(multilayerPlacementsCollection.GetNItems());
 
 for(G4int i(0); i<n; i++)
 {
  if (i!=0)
   out << '\n';
   
  out << indent << "Placement # " << i;
  
  GetPlacement(i).DumpLayout(out, indent2);
 }
}
