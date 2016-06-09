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
// File name:     RadmonDetectorMultilayersLayoutCollection.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayersLayoutCollection.cc,v 1.5 2006/06/29 16:14:09 gunter Exp $
// Tag:           $Name: geant4-09-00 $
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
