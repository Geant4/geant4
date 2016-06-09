//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonDetectorMultilayerLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerLayout.cc,v 1.6.2.2 2006/06/29 16:14:03 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//

// Include files
#include "RadmonDetectorMultilayerLayout.hh"
#include "RadmonDumpStyle.hh"
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
 G4int width(RADMONDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;

 out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << multilayerLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Size";  out.setf(std::ostream::right, std::ostream::adjustfield); out << " = (W) " << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerWidth, "Length") 
                                                                                                                                                                           << " x (H) " << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerHeight, "Length")
                                                                                                                                                                           << " x (T) " << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << G4BestUnit(GetTotalThickness(), "Length") << '\n';

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
