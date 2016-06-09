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
// File name:     RadmonGeneratorLayout.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorLayout.cc,v 1.1.2.2 2006/06/29 16:16:25 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//

// Include files
#include "RadmonGeneratorLayout.hh"

void                                            RadmonGeneratorLayout :: InsertSource(const G4String & sourceLabel)
{
 if (labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::InsertSource: Source \"" << sourceLabel << "\" just exists." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceLayout & source(labelledSourcesCollection.AppendItem());
 source.SetLabel(sourceLabel);
 source.SetIntensity(1.);
 NotifyChange();
}



void                                            RadmonGeneratorLayout :: SetRelativeSourceIntensity(const G4String & sourceLabel, G4double relativeIntensity)
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::SetRelativeSourceIntensity: Source \"" << sourceLabel << "\" not found." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));
 
 if (source.GetIntensity()!=relativeIntensity)
 {
  source.SetIntensity(relativeIntensity);
  NotifyChange();
 }
}



G4double                                        RadmonGeneratorLayout :: GetRelativeSourceIntensity(const G4String & sourceLabel) const
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetRelativeSourceIntensity: Source \"" << sourceLabel << "\" not found." << G4endl;
  return 0.;
 }
  
 const RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));
 
 return source.GetIntensity();
}



void                                            RadmonGeneratorLayout :: RemoveSource(const G4String & sourceLabel)
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::RemoveSource: Source \"" << sourceLabel << "\" not found." << G4endl;
  return;
 }
  
 labelledSourcesCollection.RemoveItemByLabel(sourceLabel);
 NotifyChange();
}



G4int                                           RadmonGeneratorLayout :: GetNSources(void) const
{
 return labelledSourcesCollection.GetNItems();
}



const G4String &                                RadmonGeneratorLayout :: GetSourceLabel(G4int index) const
{
 return labelledSourcesCollection.GetItem(index).GetLabel();
}



void                                            RadmonGeneratorLayout :: AppendSourceAlgorithm(const G4String & sourceLabel, const G4String & algorithmLabel)
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::AppendSourceAlgorithm: Source \"" << sourceLabel << "\" not found." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));

 if (source.ExistsAlgorithmByLabel(algorithmLabel))
 {
  G4cout << "RadmonGeneratorLayout::AppendSourceAlgorithm: Source \"" << algorithmLabel << "\" just exists in source \"" << sourceLabel << "\"." << G4endl;
  return;
 }

 RadmonGeneratorSourceAlgorithmLayout & algorithm(source.AppendAlgorithm());
 algorithm.SetLabel(algorithmLabel);
 NotifyChange();
}



void                                            RadmonGeneratorLayout :: SetSourceAlgorithmType(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & typeName)
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::SetSourceAlgorithmType: Source \"" << sourceLabel << "\" not found." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));

 if (! source.ExistsAlgorithmByLabel(algorithmLabel))
 {
  G4cout << "RadmonGeneratorLayout::SetSourceAlgorithmType: Source \"" << algorithmLabel << "\" not found in source \"" << sourceLabel << "\"." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceAlgorithmLayout & algorithm(source.FindAlgorithmByLabel(algorithmLabel));
 if (algorithm.GetType()!=typeName)
 {
  algorithm.SetType(typeName);
  NotifyChange();
 }
}



void                                            RadmonGeneratorLayout :: RemoveSourceAlgorithm(const G4String & sourceLabel, const G4String & algorithmLabel)
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::RemoveSourceAlgorithm: Source \"" << sourceLabel << "\" not found." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));

 if (! source.ExistsAlgorithmByLabel(algorithmLabel))
 {
  G4cout << "RadmonGeneratorLayout::RemoveSourceAlgorithm: Source \"" << algorithmLabel << "\" not found in source \"" << sourceLabel << "\"." << G4endl;
  return;
 }
  
 source.RemoveAlgorithmByLabel(algorithmLabel);
 NotifyChange();
}



G4int                                           RadmonGeneratorLayout :: GetNSourceAlgorithms(const G4String & sourceLabel) const
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetNSourceAlgorithms: Source \"" << sourceLabel << "\" not found." << G4endl;
  return 0;
 }
  
 const RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));
 
 return source.GetNAlgorithms();
}



const G4String &                                RadmonGeneratorLayout :: GetSourceAlgorithmLabel(const G4String & sourceLabel, G4int index) const
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetSourceAlgorithmLabel: Source \"" << sourceLabel << "\" not found." << G4endl;
  return GetNullStr();
 }
  
 const RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));
 
 return source.GetAlgorithm(index).GetLabel();
}



const G4String &                                RadmonGeneratorLayout :: GetSourceAlgorithmType(const G4String & sourceLabel, const G4String & algorithmLabel) const
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetSourceAlgorithmType: Source \"" << sourceLabel << "\" not found." << G4endl;
  return GetNullStr();
 }
  
 const RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));
 
 if (! source.ExistsAlgorithmByLabel(algorithmLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetSourceAlgorithmType: Source \"" << algorithmLabel << "\" not found in source \"" << sourceLabel << "\"." << G4endl;
  return GetNullStr();
 }
  
 const RadmonGeneratorSourceAlgorithmLayout & algorithm(source.FindAlgorithmByLabel(algorithmLabel));
 return algorithm.GetType();
}



void                                            RadmonGeneratorLayout :: SetSourceAlgorithmAttribute(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & attribute, const G4String & value)
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::SetSourceAlgorithmAttribute: Source \"" << sourceLabel << "\" not found." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));

 if (! source.ExistsAlgorithmByLabel(algorithmLabel))
 {
  G4cout << "RadmonGeneratorLayout::SetSourceAlgorithmAttribute: Source \"" << algorithmLabel << "\" not found in source \"" << sourceLabel << "\"." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceAlgorithmLayout & algorithm(source.FindAlgorithmByLabel(algorithmLabel));
 if (algorithm.GetAttribute(attribute, value+'#')!=value)
 {
  algorithm.SetAttribute(attribute, value);
  NotifyChange();
 }
}



void                                            RadmonGeneratorLayout :: ClearSourceAlgorithmAttribute(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & attribute)
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::ClearSourceAlgorithmAttribute: Source \"" << sourceLabel << "\" not found." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));

 if (! source.ExistsAlgorithmByLabel(algorithmLabel))
 {
  G4cout << "RadmonGeneratorLayout::ClearSourceAlgorithmAttribute: Source \"" << algorithmLabel << "\" not found in source \"" << sourceLabel << "\"." << G4endl;
  return;
 }
  
 RadmonGeneratorSourceAlgorithmLayout & algorithm(source.FindAlgorithmByLabel(algorithmLabel));
 if (algorithm.ExistsAttribute(attribute))
 {
  algorithm.ClearAttribute(attribute);
  NotifyChange();
 }
}



G4String                                        RadmonGeneratorLayout :: GetSourceAlgorithmAttribute(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & attribute, const G4String & defaultValue) const
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetSourceAlgorithmAttribute: Source \"" << sourceLabel << "\" not found." << G4endl;
  return G4String();
 }
  
 const RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));

 if (! source.ExistsAlgorithmByLabel(algorithmLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetSourceAlgorithmAttribute: Source \"" << algorithmLabel << "\" not found in source \"" << sourceLabel << "\"." << G4endl;
  return G4String();
 }
  
 const RadmonGeneratorSourceAlgorithmLayout & algorithm(source.FindAlgorithmByLabel(algorithmLabel));
 return algorithm.GetAttribute(attribute, defaultValue);
}



G4int                                            RadmonGeneratorLayout :: GetSourceAlgorithmNAttributes(const G4String & sourceLabel, const G4String & algorithmLabel) const
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetSourceAlgorithmNAttributes: Source \"" << sourceLabel << "\" not found." << G4endl;
  return 0;
 }
  
 const RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));

 if (! source.ExistsAlgorithmByLabel(algorithmLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetSourceAlgorithmNAttributes: Source \"" << algorithmLabel << "\" not found in source \"" << sourceLabel << "\"." << G4endl;
  return 0;
 }
  
 const RadmonGeneratorSourceAlgorithmLayout & algorithm(source.FindAlgorithmByLabel(algorithmLabel));
 return algorithm.GetNAttributes();
}



const G4String &                                 RadmonGeneratorLayout :: GetSourceAlgorithmAttributeName(const G4String & sourceLabel, const G4String & algorithmLabel, G4int index) const
{
 if (! labelledSourcesCollection.ExistsItemByLabel(sourceLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetSourceAlgorithmAttributeName: Source \"" << sourceLabel << "\" not found." << G4endl;
  return GetNullStr();
 }
  
 const RadmonGeneratorSourceLayout & source(labelledSourcesCollection.FindItemByLabel(sourceLabel));

 if (! source.ExistsAlgorithmByLabel(algorithmLabel))
 {
  G4cout << "RadmonGeneratorLayout::GetSourceAlgorithmAttributeName: Source \"" << algorithmLabel << "\" not found in source \"" << sourceLabel << "\"." << G4endl;
  return GetNullStr();
 }
  
 const RadmonGeneratorSourceAlgorithmLayout & algorithm(source.FindAlgorithmByLabel(algorithmLabel));
 return algorithm.GetAttributeName(index);
}



G4bool                                           RadmonGeneratorLayout :: Load(std::istream & /* in */)
{
 // TO BE DONE
 G4cout << "RadmonGeneratorLayout::Load(): PLEASE CHECK" << G4endl;

 return false; 
}



G4bool                                           RadmonGeneratorLayout :: Save(std::ostream & /* out */) const
{
 // TO BE DONE
 G4cout << "RadmonGeneratorLayout::Save(): PLEASE CHECK" << G4endl;

 return false; 
}



void                                             RadmonGeneratorLayout :: DumpLayout(std::ostream & out) const
{
 out << "- Sources:\n";
 const G4int n(labelledSourcesCollection.GetNItems());

 if (n==0)
 {
  out << "  - No sources defined.\n";
  return;
 }

 const G4String indent("    - ");

 for(G4int i(0); i<n; i++)
 {
  out << "  - Source # " << i << '\n';
   
  labelledSourcesCollection.GetItem(i).DumpLayout(out, indent);
 }
}





inline G4String &                                RadmonGeneratorLayout :: GetNullStr() const
{
 static G4String nullStr("");
 
 return nullStr;
}
