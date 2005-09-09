//
// File name:     RadmonDetectorMultilayerLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerLayout.cc,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMultilayerLayout.hh"



                                                RadmonDetectorMultilayerLayout :: RadmonDetectorMultilayerLayout(const RadmonDetectorMultilayerLayout & copy)
{
}





RadmonDetectorMultilayerLayout &                RadmonDetectorMultilayerLayout :: operator=(const RadmonDetectorMultilayerLayout & copy)
{
}





G4double                                        RadmonDetectorMultilayerLayout :: GetTotalThickness(void) const
{
}





G4int                                           RadmonDetectorMultilayerLayout :: GetNLayers(void) const
{
}



G4bool                                          RadmonDetectorMultilayerLayout :: Empty(void) const
{
}





const RadmonDetectorLayerLayout &               RadmonDetectorMultilayerLayout :: GetLayer(G4int index) const
{
}



RadmonDetectorLayerLayout &                     RadmonDetectorMultilayerLayout :: GetLayer(G4int index)
{
}





G4bool                                          RadmonDetectorMultilayerLayout :: ExistsLayerByLabel(const G4String & layerLabel) const
{
}



G4int                                           RadmonDetectorMultilayerLayout :: MultiplicityLayerByLabel(const G4String & layerLabel) const
{
}





const RadmonDetectorLayerLayout &               RadmonDetectorMultilayerLayout :: FindLayerByLabel(const G4String &layerLabel, G4int count) const
{
}



RadmonDetectorLayerLayout &                     RadmonDetectorMultilayerLayout :: FindLayerByLabel(const G4String & layerLabel, G4int count)
{
}





RadmonDetectorLayerLayout &                     RadmonDetectorMultilayerLayout :: AppendLayer(void)
{
}



RadmonDetectorLayerLayout &                     RadmonDetectorMultilayerLayout :: PrependLayer(void)
{
}





void                                            RadmonDetectorMultilayerLayout :: RemoveLayerByLabel(const G4String & layerName, G4int count)
{
}



void                                            RadmonDetectorMultilayerLayout :: RemoveLayersByLabel(const G4String & layerName)
{
}



void                                            RadmonDetectorMultilayerLayout :: RemoveLayer(G4int index)
{
}



void                                            RadmonDetectorMultilayerLayout :: RemoveLayersByRange(G4int first, G4int last)
{
}



void                                            RadmonDetectorMultilayerLayout :: RemoveAllLayers(void)
{
}





void                                            RadmonDetectorMultilayerLayout :: DumpLayout(std::ostream & out) const
{
}
