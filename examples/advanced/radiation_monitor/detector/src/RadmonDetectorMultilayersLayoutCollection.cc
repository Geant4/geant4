//
// File name:     RadmonDetectorMultilayersLayoutCollection.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayersLayoutCollection.cc,v 1.3 2005-10-24 14:51:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMultilayersLayoutCollection.hh"
#include "RadmonDumpStyle.hh"
#include <iomanip>



G4int                                           RadmonDetectorMultilayersLayoutCollection :: GetNMultilayers(void) const
{
 return multilayersCollection.GetNItems();
}



G4bool                                          RadmonDetectorMultilayersLayoutCollection :: Empty(void) const
{
 return multilayersCollection.Empty();
}





const RadmonDetectorMultilayerLayout &          RadmonDetectorMultilayersLayoutCollection :: GetMultilayer(G4int index) const
{
 return multilayersCollection.GetItem(index);
}



RadmonDetectorMultilayerLayout &                RadmonDetectorMultilayersLayoutCollection :: GetMultilayer(G4int index)
{
 return multilayersCollection.GetItem(index);
}





G4bool                                          RadmonDetectorMultilayersLayoutCollection :: ExistsMultilayerByLabel(const G4String & label) const
{
 return multilayersCollection.ExistsItemByLabel(label);
}



G4int                                           RadmonDetectorMultilayersLayoutCollection :: MultiplicityMultilayerByLabel(const G4String & label) const
{
 return multilayersCollection.MultiplicityItemByLabel(label);
}





const RadmonDetectorMultilayerLayout &          RadmonDetectorMultilayersLayoutCollection :: FindMultilayerByLabel(const G4String & label, G4int count) const
{
 return multilayersCollection.FindItemByLabel(label, count);
}



RadmonDetectorMultilayerLayout &                RadmonDetectorMultilayersLayoutCollection :: FindMultilayerByLabel(const G4String & label, G4int count)
{
 return multilayersCollection.FindItemByLabel(label, count);
}





RadmonDetectorMultilayerLayout &                RadmonDetectorMultilayersLayoutCollection :: CreateMultilayer(void)
{
 return multilayersCollection.AppendItem();
}





void                                            RadmonDetectorMultilayersLayoutCollection :: RemoveMultilayerByLabel(const G4String & label, G4int count)
{
 multilayersCollection.RemoveItemByLabel(label, count);
}



void                                            RadmonDetectorMultilayersLayoutCollection :: RemoveMultilayersByLabel(const G4String & label)
{
 multilayersCollection.RemoveItemsByLabel(label);
}



void                                            RadmonDetectorMultilayersLayoutCollection :: RemoveMultilayer(G4int index)
{
 multilayersCollection.RemoveItem(index);
}



void                                            RadmonDetectorMultilayersLayoutCollection :: RemoveAllMultilayers(void)
{
 multilayersCollection.RemoveAllItems();
}





void                                            RadmonDetectorMultilayersLayoutCollection :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 const G4int n(multilayersCollection.GetNItems());
 
 if (n==0)
 {
  out << indent << "No multilayers defined.\n";
  return;
 }

 G4String indent2(indent);
 indent2.prepend("  ");

 for(G4int i(0); i<n; i++)
 {
  if (i!=0)
   out << '\n';
   
  out << indent << "Multilayer # " << i << '\n';
  
  GetMultilayer(i).DumpLayout(out, indent2);
 }
}
