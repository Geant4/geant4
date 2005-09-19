//
// File name:     RadmonDetectorMultilayerLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerLayout.cc,v 1.4 2005-09-19 19:42:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMultilayerLayout.hh"
#include "RadmonDetectorDumpStyle.hh"
#include "G4UnitsTable.hh"

#include <iomanip>



                                                RadmonDetectorMultilayerLayout :: RadmonDetectorMultilayerLayout(const RadmonDetectorMultilayerLayout & copy)
:
 multilayerLabel(copy.multilayerLabel),
 multilayerWidth(copy.multilayerWidth),
 multilayerHeight(copy.multilayerHeight),
 multilayerLayersCollection(copy.multilayerLayersCollection)
{
}





RadmonDetectorMultilayerLayout &                RadmonDetectorMultilayerLayout :: operator=(const RadmonDetectorMultilayerLayout & copy)
{
 multilayerLabel=copy.multilayerLabel;
 multilayerWidth=copy.multilayerWidth;
 multilayerHeight=copy.multilayerHeight;
 multilayerLayersCollection=copy.multilayerLayersCollection;
 
 return (*this);
}





G4double                                        RadmonDetectorMultilayerLayout :: GetTotalThickness(void) const
{
 G4int i(multilayerLayersCollection.GetNItems());
 G4double thickness(0);
 
 while (i>0)
 {
  i--;
  thickness+=GetLayer(i).GetThickness();
 }
 
 return thickness;
}





G4int                                           RadmonDetectorMultilayerLayout :: GetNLayers(void) const
{
 return multilayerLayersCollection.GetNItems();
}



G4bool                                          RadmonDetectorMultilayerLayout :: Empty(void) const
{
 return multilayerLayersCollection.Empty();
}





const RadmonDetectorLayerLayout &               RadmonDetectorMultilayerLayout :: GetLayer(G4int index) const
{
 return multilayerLayersCollection.GetItem(index);
}



RadmonDetectorLayerLayout &                     RadmonDetectorMultilayerLayout :: GetLayer(G4int index)
{
 return multilayerLayersCollection.GetItem(index);
}





G4bool                                          RadmonDetectorMultilayerLayout :: ExistsLayerByLabel(const G4String & layerLabel) const
{
 return multilayerLayersCollection.ExistsItemByLabel(layerLabel);
}



G4int                                           RadmonDetectorMultilayerLayout :: MultiplicityLayerByLabel(const G4String & layerLabel) const
{
 return multilayerLayersCollection.MultiplicityItemByLabel(layerLabel);
}





const RadmonDetectorLayerLayout &               RadmonDetectorMultilayerLayout :: FindLayerByLabel(const G4String & layerLabel, G4int count) const
{
 return multilayerLayersCollection.FindItemByLabel(layerLabel, count);
}



RadmonDetectorLayerLayout &                     RadmonDetectorMultilayerLayout :: FindLayerByLabel(const G4String & layerLabel, G4int count)
{
 return multilayerLayersCollection.FindItemByLabel(layerLabel, count);
}





RadmonDetectorLayerLayout &                     RadmonDetectorMultilayerLayout :: AppendLayer(void)
{
 return multilayerLayersCollection.AppendItem();
}



RadmonDetectorLayerLayout &                     RadmonDetectorMultilayerLayout :: PrependLayer(void)
{
 return multilayerLayersCollection.PrependItem();
}





void                                            RadmonDetectorMultilayerLayout :: RemoveLayerByLabel(const G4String & layerLabel, G4int count)
{
 multilayerLayersCollection.RemoveItemByLabel(layerLabel, count);
}



void                                            RadmonDetectorMultilayerLayout :: RemoveLayersByLabel(const G4String & layerLabel)
{
 multilayerLayersCollection.RemoveItemsByLabel(layerLabel);
}



void                                            RadmonDetectorMultilayerLayout :: RemoveLayer(G4int index)
{
 multilayerLayersCollection.RemoveItem(index);
}



void                                            RadmonDetectorMultilayerLayout :: RemoveLayersByRange(G4int first, G4int last)
{
 multilayerLayersCollection.RemoveItemsByRange(first, last);
}



void                                            RadmonDetectorMultilayerLayout :: RemoveAllLayers(void)
{
 multilayerLayersCollection.RemoveAllItems();
}





void                                            RadmonDetectorMultilayerLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 G4int width(RADMONDETECTORDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;

 out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << multilayerLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Size";  out.setf(std::ostream::right, std::ostream::adjustfield); out << " = (W) " << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerWidth, "Length") 
                                                                                                                                                                           << " x (T) " << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << G4BestUnit(GetTotalThickness(), "Length")
                                                                                                                                                                           << " x (H) " << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerHeight, "Length") << '\n';

 G4String indent2(indent);
 indent2.prepend("  ");

 const G4int n(multilayerLayersCollection.GetNItems());

 if (n==0)
 {
  out << indent2 << "No layers defined.\n";
  return;
 }

 G4String indent3(indent2);
 indent3.prepend("  ");

 for(G4int i(0); i<n; i++)
 {
  out << indent2 << "Layer # " << i << '\n';
   
  GetLayer(i).DumpLayout(out, indent3);
 }
}
