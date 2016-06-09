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
// File name:     RadmonDataAnalysisLayout.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysisLayout.cc,v 1.3 2006/06/29 16:07:37 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//

// Include files
#include "RadmonAnalysisLayout.hh"
#include "RadmonDumpStyle.hh"
#include "G4UnitsTable.hh"

#include <iomanip>


                                                RadmonAnalysisLayout :: RadmonAnalysisLayout()
{
}
 
 
 
                                                RadmonAnalysisLayout :: ~RadmonAnalysisLayout()
{
}





void                                            RadmonAnalysisLayout :: SetOutputFileName(const G4String & outputFileName)
{
 if (fileName==outputFileName)
  return;
  
 fileName=outputFileName;
 NotifyChange();
}



const G4String &                                RadmonAnalysisLayout :: GetOutputFileName(void) const
{
 return fileName;
}



void                                            RadmonAnalysisLayout :: SetOutputFileFormat(const G4String & outputFileFormat)
{
 if (fileFormat==outputFileFormat)
  return;
  
 fileFormat=outputFileFormat;
 NotifyChange();
}



const G4String &                                RadmonAnalysisLayout :: GetOutputFileFormat(void) const
{
 return fileFormat;
}





void                                            RadmonAnalysisLayout :: CreateSensitiveDetector(const G4String & sensitiveDetectorLabel, const G4String & sensitiveDetectorType)
{
 if (ExistsSensitiveDetector(sensitiveDetectorLabel)) 
  return;
  
 if (!FindSensitiveDetectorType(sensitiveDetectorType))
  return;
 
 sensitiveDetectors.push_back(SensitiveDetector(sensitiveDetectorLabel, sensitiveDetectorType));
 NotifyChange();
}



void                                            RadmonAnalysisLayout :: SetSensitiveDetectorType(const G4String & sensitiveDetectorLabel, const G4String & sensitiveDetectorType)
{
 SensitiveDetector * sensitiveDetector(FindSensitiveDetector(sensitiveDetectorLabel));
 
 if (!sensitiveDetector)
  return;
  
 if (sensitiveDetectorType==sensitiveDetector->second)
  return;
  
 if (!FindSensitiveDetectorType(sensitiveDetectorType))
  return;
 
 sensitiveDetector->second=sensitiveDetectorType;
 NotifyChange();
}



const G4String &                                RadmonAnalysisLayout :: GetSensitiveDetectorType(const G4String & sensitiveDetectorLabel) const
{
 const SensitiveDetector * sensitiveDetector(FindSensitiveDetector(sensitiveDetectorLabel));
 
 if (!sensitiveDetector)
  return GetNullStr();
 
 return sensitiveDetector->second;
}



void                                            RadmonAnalysisLayout :: RemoveSensitiveDetector(const G4String & sensitiveDetectorLabel)
{
 SensitiveDetectors::iterator i(sensitiveDetectors.begin());
 const SensitiveDetectors::iterator end(sensitiveDetectors.end());
 
 while (i<end)
 {
  if (i->first==sensitiveDetectorLabel)
  {
   sensitiveDetectors.erase(i);
   return;
  }
  
  i++;
 }  

 G4cout << "RadmonAnalysisLayout::RemoveSensitiveDetector: Sensitive detector \"" << sensitiveDetectorLabel << "\" does not exist." << G4endl;
}



G4int                                           RadmonAnalysisLayout :: GetNSensitiveDetectors(void) const
{
 return sensitiveDetectors.size();
}



const G4String &                                RadmonAnalysisLayout :: GetSensitiveDetectorLabel(G4int index) const
{
 return sensitiveDetectors[index].first;
}





void                                            RadmonAnalysisLayout :: CreateSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel)
{
 if (sensitiveDetectorTypes.ExistsItemByLabel(sensitiveDetectorTypeLabel))
 {
  G4cout << "RadmonAnalysisLayout::CreateSensitiveDetectorType: Sensitive detector \"" << sensitiveDetectorTypeLabel << "\" just exists." << G4endl;
  return;
 }
 
 RadmonAnalysisSensitiveDetectorTypeLayout & sensitiveDetectorType(sensitiveDetectorTypes.AppendItem());
 sensitiveDetectorType.SetLabel(sensitiveDetectorTypeLabel);
}



void                                            RadmonAnalysisLayout :: RemoveSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel)
{
 if (!sensitiveDetectorTypes.ExistsItemByLabel(sensitiveDetectorTypeLabel))
 {
  G4cout << "RadmonAnalysisLayout::RemoveSensitiveDetectorType: Sensitive detector \"" << sensitiveDetectorTypeLabel << "\" does not exist." << G4endl;
  return;
 }
 
 if (SensitiveDetectorTypeInUse(sensitiveDetectorTypeLabel))
 {
  G4cout << "RadmonAnalysisLayout::RemoveSensitiveDetectorType: Sensitive detector \"" << sensitiveDetectorTypeLabel << "\" is in use. Cannot be removed." << G4endl;
  return;
 }

 sensitiveDetectorTypes.RemoveItemByLabel(sensitiveDetectorTypeLabel);
}





void                                            RadmonAnalysisLayout :: AppendDataAnalysisToSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel)
{
 RadmonAnalysisSensitiveDetectorTypeLayout * sensitiveDetectorType(FindSensitiveDetectorType(sensitiveDetectorTypeLabel));
 
 if (! sensitiveDetectorType)
  return;
  
 if (sensitiveDetectorType->ExistsDataAnalysisByLabel(dataAnalysisLabel))
 {
  G4cout << "RadmonAnalysisLayout::AppendDataAnalysisToSensitiveDetectorType: Data analysis label \"" << dataAnalysisLabel << "\" just exits in sensitive detector \"" << sensitiveDetectorTypeLabel << "\"." << G4endl;
  return;
 }
 
 RadmonDataAnalysisLayout & dataAnalysisLayout(sensitiveDetectorType->AppendDataAnalysis());
 dataAnalysisLayout.SetLabel(dataAnalysisLabel);

 if (SensitiveDetectorTypeInUse(sensitiveDetectorTypeLabel))
  NotifyChange();
}



void                                            RadmonAnalysisLayout :: SetDataAnalysisType(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & dataAnalysisType)
{
 RadmonDataAnalysisLayout * dataAnalysis(FindDataAnalysis(sensitiveDetectorTypeLabel, dataAnalysisLabel));
 
 if (!dataAnalysis)
  return;
 
 if (dataAnalysis->GetType()==dataAnalysisType)
  return;
  
 dataAnalysis->SetType(dataAnalysisType);
 
 if (SensitiveDetectorTypeInUse(sensitiveDetectorTypeLabel))
  NotifyChange();
}



void                                            RadmonAnalysisLayout :: RemoveDataAnalysis(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel)
{
 RadmonAnalysisSensitiveDetectorTypeLayout * sensitiveDetectorType(FindSensitiveDetectorType(sensitiveDetectorTypeLabel));
 
 if (! sensitiveDetectorType)
  return;
  
 if (!sensitiveDetectorType->ExistsDataAnalysisByLabel(dataAnalysisLabel))
 {
  G4cout << "RadmonAnalysisLayout::RemoveDataAnalysis: Data analysis label \"" << dataAnalysisLabel << "\" does not exit in sensitive detector \"" << sensitiveDetectorTypeLabel << "\"." << G4endl;
  return;
 }
 
 sensitiveDetectorType->RemoveDataAnalysisByLabel(dataAnalysisLabel);
 
 if (SensitiveDetectorTypeInUse(sensitiveDetectorTypeLabel))
  NotifyChange();
}



G4int                                           RadmonAnalysisLayout :: GetNDataAnalyses(const G4String & sensitiveDetectorTypeLabel) const
{
 const RadmonAnalysisSensitiveDetectorTypeLayout * sensitiveDetectorType(FindSensitiveDetectorType(sensitiveDetectorTypeLabel));
 
 if (! sensitiveDetectorType)
  return 0;
 
 return sensitiveDetectorType->GetNDataAnalyses();
}



const G4String &                                RadmonAnalysisLayout :: GetDataAnalysisLabel(const G4String & sensitiveDetectorTypeLabel, G4int index) const
{
 const RadmonAnalysisSensitiveDetectorTypeLayout * sensitiveDetectorType(FindSensitiveDetectorType(sensitiveDetectorTypeLabel));
 
 if (! sensitiveDetectorType)
  return GetNullStr();
 
 return sensitiveDetectorType->GetDataAnalysis(index).GetLabel();
}



const G4String &                                RadmonAnalysisLayout :: GetDataAnalysisType(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) const
{
 const RadmonDataAnalysisLayout * dataAnalysisType(FindDataAnalysis(sensitiveDetectorTypeLabel, dataAnalysisLabel));
 
 if (!dataAnalysisType)
  return GetNullStr();
 
 return dataAnalysisType->GetType();
}





void                                            RadmonAnalysisLayout :: SetDataAnalysisAttribute(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & attributeName, const G4String & attributeValue)
{
 RadmonDataAnalysisLayout * dataAnalysisType(FindDataAnalysis(sensitiveDetectorTypeLabel, dataAnalysisLabel));
 
 if (!dataAnalysisType)
  return;
  
 dataAnalysisType->SetAttribute(attributeName, attributeValue);
 
 if (SensitiveDetectorTypeInUse(sensitiveDetectorTypeLabel))
  NotifyChange();
}



void                                            RadmonAnalysisLayout :: ClearDataAnalysisAttribute(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & attributeName)
{
 RadmonDataAnalysisLayout * dataAnalysisType(FindDataAnalysis(sensitiveDetectorTypeLabel, dataAnalysisLabel));
 
 if (!dataAnalysisType)
  return;
  
 if (!dataAnalysisType->ExistsAttribute(attributeName))
 {
  G4cout << "RadmonAnalysisLayout::ClearDataAnalysisAttribute: Attribute \"" << attributeName << "\" does not exist in data analysis \"" << dataAnalysisLabel << "\" of sensitive detector \"" << sensitiveDetectorTypeLabel << "\"." << G4endl;
  return;
 } 
 
 dataAnalysisType->ClearAttribute(attributeName);
 
 if (SensitiveDetectorTypeInUse(sensitiveDetectorTypeLabel))
  NotifyChange();
}



G4int                                           RadmonAnalysisLayout :: GetDataAnalysisNAttributes(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) const
{
 const RadmonDataAnalysisLayout * dataAnalysisType(FindDataAnalysis(sensitiveDetectorTypeLabel, dataAnalysisLabel));
 
 if (!dataAnalysisType)
  return 0;
  
 return dataAnalysisType->GetNAttributes();
}



G4String                                        RadmonAnalysisLayout :: GetDataAnalysisAttributeName(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, G4int index) const
{
 const RadmonDataAnalysisLayout * dataAnalysisType(FindDataAnalysis(sensitiveDetectorTypeLabel, dataAnalysisLabel));
 
 if (!dataAnalysisType)
  return GetNullStr();
  
 return dataAnalysisType->GetAttributeName(index);
}



G4String                                        RadmonAnalysisLayout :: GetDataAnalysisAttribute(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel, const G4String & attributeName, const G4String & defaultAttributeValue) const
{
 const RadmonDataAnalysisLayout * dataAnalysisType(FindDataAnalysis(sensitiveDetectorTypeLabel, dataAnalysisLabel));
 
 if (!dataAnalysisType)
  return G4String();
  
 return dataAnalysisType->GetAttribute(attributeName, defaultAttributeValue);
}





void                                            RadmonAnalysisLayout :: DumpLayout(std::ostream & out) const
{
 const G4String indent("  - ");

 out << "- " << std::setw(RADMONDUMP_INDENT_WIDTH-2); out.setf(std::ostream::left, std::ostream::adjustfield); out << "File name";   out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << fileName << "\"\n"
        "  " << std::setw(RADMONDUMP_INDENT_WIDTH-2); out.setf(std::ostream::left, std::ostream::adjustfield); out << "File format"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << fileFormat << "\"\n"
        "\n- Sensitive detectors\n";

 G4String indent2(indent);
 indent2.prepend("  ");

 if (sensitiveDetectors.empty())
  out << indent << "No sensitive detectors defined.\n";
 else
 {
  const G4int width(RADMONDUMP_INDENT_WIDTH-indent2.length());
  const G4int n(sensitiveDetectors.size());
 
  for (G4int i(0); i<n; i++)
  {
   out << indent << "Sensitive Detector # " << i << '\n'
       << indent2 << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << sensitiveDetectors[i].first << "\"\n"
       << indent2 << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Type"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << sensitiveDetectors[i].second << "\"\n";
  }
 }
 
 out << "\n- Sensitive detector types\n";
 G4int n(sensitiveDetectorTypes.GetNItems());
 if (n==0)
 {
  out << indent << "No sensitive detector types defined.\n";
  return;
 }

 for(G4int i(0); i<n; i++)
 {
  out << indent << "Sensitive detector type # " << i << '\n';
   
  sensitiveDetectorTypes.GetItem(i).DumpLayout(out, indent2);
 }
}





G4bool                                          RadmonAnalysisLayout :: Load(std::istream & /* in */)
{
 // TO BE DONE
 G4cout << "RadmonAnalysisLayout::Load(): PLEASE CHECK" << G4endl;

 return false; 
}



G4bool                                          RadmonAnalysisLayout :: Save(std::ostream & /* out */) const
{
 // TO BE DONE
 G4cout << "RadmonAnalysisLayout::Save(): PLEASE CHECK" << G4endl;

 return false; 
}





inline G4bool                                   RadmonAnalysisLayout :: SensitiveDetectorTypeInUse(const G4String & sensitiveDetectorType)
{
 G4int i(0);
 G4int n(sensitiveDetectors.size());
 
 for (; i<n; i++)
  if (sensitiveDetectors[i].second==sensitiveDetectorType)
   return true;

 return false; 
}



inline G4bool                                   RadmonAnalysisLayout :: ExistsSensitiveDetector(const G4String & sensitiveDetectorLabel) const
{
 G4int i(0);
 G4int n(sensitiveDetectors.size());
 
 for (; i<n; i++)
  if (sensitiveDetectors[i].first==sensitiveDetectorLabel)
  {
   G4cout << "RadmonAnalysisLayout::ExistsSensitiveDetector: Sensitive detector \"" << sensitiveDetectorLabel << "\" just exists." << G4endl;
   return true;
  }

 return false; 
}



inline RadmonAnalysisLayout::SensitiveDetector * RadmonAnalysisLayout :: FindSensitiveDetector(const G4String & sensitiveDetectorLabel)
{
 G4int i(0);
 G4int n(sensitiveDetectors.size());
 
 for (; i<n; i++)
  if (sensitiveDetectors[i].first==sensitiveDetectorLabel)
   return &(sensitiveDetectors[i]);

 G4cout << "RadmonAnalysisLayout::FindSensitiveDetector: Sensitive detector \"" << sensitiveDetectorLabel << "\" does not exist." << G4endl;
 return 0;
}



inline const RadmonAnalysisLayout::SensitiveDetector * RadmonAnalysisLayout :: FindSensitiveDetector(const G4String & sensitiveDetectorLabel) const
{
 return const_cast<RadmonAnalysisLayout *>(this)->FindSensitiveDetector(sensitiveDetectorLabel);
}



inline RadmonAnalysisSensitiveDetectorTypeLayout * RadmonAnalysisLayout :: FindSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel)
{
 if (!sensitiveDetectorTypes.ExistsItemByLabel(sensitiveDetectorTypeLabel))
 {
  G4cout << "RadmonAnalysisLayout::FindSensitiveDetectorType: Sensitive detector type \"" << sensitiveDetectorTypeLabel << "\" does not exist." << G4endl;
  return 0;
 }

 return & sensitiveDetectorTypes.FindItemByLabel(sensitiveDetectorTypeLabel);
}



inline const RadmonAnalysisSensitiveDetectorTypeLayout * RadmonAnalysisLayout :: FindSensitiveDetectorType(const G4String & sensitiveDetectorTypeLabel) const
{
 return const_cast<RadmonAnalysisLayout *>(this)->FindSensitiveDetectorType(sensitiveDetectorTypeLabel);
}



inline RadmonDataAnalysisLayout *               RadmonAnalysisLayout :: FindDataAnalysis(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel)
{
 RadmonAnalysisSensitiveDetectorTypeLayout * sensitiveDetectorType(FindSensitiveDetectorType(sensitiveDetectorTypeLabel));
 
 if (sensitiveDetectorType==0)
  return 0;

 if (!sensitiveDetectorType->ExistsDataAnalysisByLabel(dataAnalysisLabel))
 {
  G4cout << "RadmonAnalysisLayout::FindDataAnalysis: Data analysis labelled \"" << dataAnalysisLabel << "\" does not exist in sensitive detector type \"" << sensitiveDetectorTypeLabel << "\"." << G4endl;
  return 0;
 }

 return & sensitiveDetectorType->FindDataAnalysisByLabel(dataAnalysisLabel);
}



inline const RadmonDataAnalysisLayout *         RadmonAnalysisLayout :: FindDataAnalysis(const G4String & sensitiveDetectorTypeLabel, const G4String & dataAnalysisLabel) const
{
 return const_cast<RadmonAnalysisLayout *>(this)->FindDataAnalysis(sensitiveDetectorTypeLabel, dataAnalysisLabel);
}





inline G4String &                               RadmonAnalysisLayout :: GetNullStr() const
{
 static G4String *nullStr(0);
 
 if (nullStr==0)
  nullStr=new G4String("");
  
 return *nullStr;
}
