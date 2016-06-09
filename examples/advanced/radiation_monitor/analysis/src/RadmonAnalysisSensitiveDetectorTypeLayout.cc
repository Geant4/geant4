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
// File name:     RadmonAnalysisSensitiveDetectorTypeLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysisSensitiveDetectorTypeLayout.cc,v 1.3 2006/06/29 16:07:41 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//

// Include files
#include "RadmonAnalysisSensitiveDetectorTypeLayout.hh"
#include "RadmonDumpStyle.hh"
#include "G4UnitsTable.hh"

#include <iomanip>



                                                RadmonAnalysisSensitiveDetectorTypeLayout :: RadmonAnalysisSensitiveDetectorTypeLayout(const RadmonAnalysisSensitiveDetectorTypeLayout & copy)
:
 sensitiveDetectorLabel(copy.sensitiveDetectorLabel),
 dataAnalysesCollection(copy.dataAnalysesCollection)
{
}





RadmonAnalysisSensitiveDetectorTypeLayout &     RadmonAnalysisSensitiveDetectorTypeLayout :: operator=(const RadmonAnalysisSensitiveDetectorTypeLayout & copy)
{
 sensitiveDetectorLabel=copy.sensitiveDetectorLabel;
 dataAnalysesCollection=copy.dataAnalysesCollection;
 
 return (*this);
}





G4int                                           RadmonAnalysisSensitiveDetectorTypeLayout :: GetNDataAnalyses(void) const
{
 return dataAnalysesCollection.GetNItems();
}



G4bool                                          RadmonAnalysisSensitiveDetectorTypeLayout :: Empty(void) const
{
 return dataAnalysesCollection.Empty();
}





const RadmonDataAnalysisLayout &                RadmonAnalysisSensitiveDetectorTypeLayout :: GetDataAnalysis(G4int index) const
{
 return dataAnalysesCollection.GetItem(index);
}



RadmonDataAnalysisLayout &                      RadmonAnalysisSensitiveDetectorTypeLayout :: GetDataAnalysis(G4int index)
{
 return dataAnalysesCollection.GetItem(index);
}





G4bool                                          RadmonAnalysisSensitiveDetectorTypeLayout :: ExistsDataAnalysisByLabel(const G4String & dataAnalysisLabel) const
{
 return dataAnalysesCollection.ExistsItemByLabel(dataAnalysisLabel);
}



G4int                                           RadmonAnalysisSensitiveDetectorTypeLayout :: MultiplicityDataAnalysisByLabel(const G4String & dataAnalysisLabel) const
{
 return dataAnalysesCollection.MultiplicityItemByLabel(dataAnalysisLabel);
}





const RadmonDataAnalysisLayout &                RadmonAnalysisSensitiveDetectorTypeLayout :: FindDataAnalysisByLabel(const G4String & dataAnalysisLabel, G4int count) const
{
 return dataAnalysesCollection.FindItemByLabel(dataAnalysisLabel, count);
}



RadmonDataAnalysisLayout &                      RadmonAnalysisSensitiveDetectorTypeLayout :: FindDataAnalysisByLabel(const G4String & dataAnalysisLabel, G4int count)
{
 return dataAnalysesCollection.FindItemByLabel(dataAnalysisLabel, count);
}





RadmonDataAnalysisLayout &                      RadmonAnalysisSensitiveDetectorTypeLayout :: AppendDataAnalysis(void)
{
 return dataAnalysesCollection.AppendItem();
}



RadmonDataAnalysisLayout &                      RadmonAnalysisSensitiveDetectorTypeLayout :: PrependDataAnalysis(void)
{
 return dataAnalysesCollection.PrependItem();
}





void                                            RadmonAnalysisSensitiveDetectorTypeLayout :: RemoveDataAnalysisByLabel(const G4String & dataAnalysisLabel, G4int count)
{
 dataAnalysesCollection.RemoveItemByLabel(dataAnalysisLabel, count);
}



void                                            RadmonAnalysisSensitiveDetectorTypeLayout :: RemoveDataAnalysesByLabel(const G4String & dataAnalysisLabel)
{
 dataAnalysesCollection.RemoveItemsByLabel(dataAnalysisLabel);
}



void                                            RadmonAnalysisSensitiveDetectorTypeLayout :: RemoveDataAnalysis(G4int index)
{
 dataAnalysesCollection.RemoveItem(index);
}



void                                            RadmonAnalysisSensitiveDetectorTypeLayout :: RemoveDataAnalysesByRange(G4int first, G4int last)
{
 dataAnalysesCollection.RemoveItemsByRange(first, last);
}



void                                            RadmonAnalysisSensitiveDetectorTypeLayout :: RemoveAllDataAnalyses(void)
{
 dataAnalysesCollection.RemoveAllItems();
}





void                                            RadmonAnalysisSensitiveDetectorTypeLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 G4int width(RADMONDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;

 out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << sensitiveDetectorLabel << "\"\n";

 G4String indent2(indent);
 indent2.prepend("  ");

 const G4int n(dataAnalysesCollection.GetNItems());

 if (n==0)
 {
  out << indent2 << "No data analyses defined.\n";
  return;
 }

 G4String indent3(indent2);
 indent3.prepend("  ");

 for(G4int i(0); i<n; i++)
 {
  out << indent2 << "Data analysis # " << i << '\n';
   
  GetDataAnalysis(i).DumpLayout(out, indent3);
 }
}
