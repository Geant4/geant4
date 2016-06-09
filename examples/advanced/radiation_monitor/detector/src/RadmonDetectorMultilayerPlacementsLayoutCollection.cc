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
// File name:     RadmonDetectorMultilayerPlacementsLayoutCollection.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementsLayoutCollection.cc,v 1.4 2006/06/29 16:14:07 gunter Exp $
// Tag:           $Name: geant4-08-01 $
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
 
 if (n==0)
  out << indent << "No placements defined.\n";

 for(G4int i(0); i<n; i++)
 {
  if (i!=0)
   out << '\n';
   
  out << indent << "Placement # " << i << '\n';
  
  GetPlacement(i).DumpLayout(out, indent2);
 }
}
