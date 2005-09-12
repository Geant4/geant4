//
// File name:     RadmonDetectorMultilayerLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerLayout.cc,v 1.2 2005-09-12 17:14:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMultilayerLayout.hh"
#include "RadmonDetectorDumpStyle.hh"
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
 size_t width(RADMONDETECTORDUMPWIDTH-indent.length());
 out << indent << std::setw(width) << "Label" << " = \"" << multilayerLabel << "\"\n"
     << indent << std::setw(width) << "Size" << " = "  << std::setprecision(2) << multilayerWidth/mm << " mm x " << std::setprecision(2) << multilayerHeight/mm << " mm\n";

 G4String indent2(indent);
 indent2.prepend("  ");

 const G4int n(multilayerLayersCollection.GetNItems());
 
 for(G4int i(0); i<n; i++)
  GetLayer(i).DumpLayout(out, indent2);
}
