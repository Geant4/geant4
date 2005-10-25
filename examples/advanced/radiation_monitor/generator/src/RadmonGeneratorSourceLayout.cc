//
// File name:     RadmonGeneratorSourceLayout.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorSourceLayout.cc,v 1.1 2005-10-25 16:36:43 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonGeneratorSourceLayout.hh"
#include "RadmonDumpStyle.hh"

#include <iomanip>



                                                RadmonGeneratorSourceLayout :: RadmonGeneratorSourceLayout(const RadmonGeneratorSourceLayout & copy)
:
 algorithmsCollection(copy.algorithmsCollection),
 sourceLabel(copy.sourceLabel),
 sourceIntensity(copy.sourceIntensity)
{
}





RadmonGeneratorSourceLayout &                   RadmonGeneratorSourceLayout :: operator=(const RadmonGeneratorSourceLayout & copy)
{
 algorithmsCollection=copy.algorithmsCollection;
 sourceLabel=copy.sourceLabel;
 sourceIntensity=copy.sourceIntensity;
 
 return (*this);
}





RadmonGeneratorSourceAlgorithmLayout &          RadmonGeneratorSourceLayout :: AppendAlgorithm(void)
{
 return algorithmsCollection.AppendItem();
}





RadmonGeneratorSourceAlgorithmLayout &          RadmonGeneratorSourceLayout :: GetAlgorithm(G4int index)
{
 return algorithmsCollection.GetItem(index);
}



const RadmonGeneratorSourceAlgorithmLayout &    RadmonGeneratorSourceLayout :: GetAlgorithm(G4int index) const
{
 return algorithmsCollection.GetItem(index);
}



G4bool                                          RadmonGeneratorSourceLayout :: ExistsAlgorithmByLabel(const G4String & label) const
{
 return algorithmsCollection.ExistsItemByLabel(label);
}



G4int                                           RadmonGeneratorSourceLayout :: MultiplicityAlgorithmByLabel(const G4String & label) const
{
 return algorithmsCollection.MultiplicityItemByLabel(label);
}



RadmonGeneratorSourceAlgorithmLayout &          RadmonGeneratorSourceLayout :: FindAlgorithmByLabel(const G4String & label, G4int count)
{
 return algorithmsCollection.FindItemByLabel(label, count);
}



const RadmonGeneratorSourceAlgorithmLayout &    RadmonGeneratorSourceLayout :: FindAlgorithmByLabel(const G4String & label, G4int count) const
{
 return algorithmsCollection.FindItemByLabel(label, count);
}





void                                            RadmonGeneratorSourceLayout :: RemoveAlgorithmByLabel(const G4String & label, G4int count)
{
 return algorithmsCollection.RemoveItemByLabel(label, count);
}



void                                            RadmonGeneratorSourceLayout :: RemoveAlgorithmsByLabel(const G4String & label)
{
 return algorithmsCollection.RemoveItemsByLabel(label);
}



void                                            RadmonGeneratorSourceLayout :: RemoveAlgorithm(G4int index)
{
 return algorithmsCollection.RemoveItem(index);
}



void                                            RadmonGeneratorSourceLayout :: RemoveAlgorithmsByRange(G4int first, G4int last)
{
 return algorithmsCollection.RemoveItemsByRange(first, last);
}



void                                            RadmonGeneratorSourceLayout :: RemoveAllAlgorithms(void)
{
 return algorithmsCollection.RemoveAllItems();
}





void                                            RadmonGeneratorSourceLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 G4int width(RADMONDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;

 out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << sourceLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Insensity";  out.setf(std::ostream::right, std::ostream::adjustfield); out << " = " << sourceIntensity << '\n';

 G4String indent2(indent);
 indent2.prepend("  ");

 const G4int n(algorithmsCollection.GetNItems());

 if (n==0)
 {
  out << indent2 << "No algorithms defined.\n";
  return;
 }

 G4String indent3(indent2);
 indent3.prepend("  ");

 for(G4int i(0); i<n; i++)
 {
  out << indent2 << "Algorithm # " << i << '\n';
   
  GetAlgorithm(i).DumpLayout(out, indent3);
 }
}
