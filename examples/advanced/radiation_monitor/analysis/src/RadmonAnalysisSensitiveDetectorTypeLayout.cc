//
// File name:     RadmonAnalysisSensitiveDetectorTypeLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysisSensitiveDetectorTypeLayout.cc,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
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
