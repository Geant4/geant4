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
#ifndef DicomFileMgr__HH
#define DicomFileMgr__HH
#include <vector>
#include <map>
#include "globals.hh"

#include "dcmtk/dcmdata/dcfilefo.h"
class DicomVFile;
class DicomFileCT;
class DicomFileStructure;
class DicomFilePlan;
class DicomFilePET;
class DcmDataset;

//typedef std::multimap<OFString,DicomVFile*> msd;
typedef std::map<G4double,DicomFileCT*> mdct;
typedef std::map<G4double,DicomFilePET*> mdpet;
enum VerbLevel {silentVerb = -1, errorVerb = 0, warningVerb = 1, infoVerb=2, debugVerb=3, 
 testVerb=4};

class DicomFileMgr 
{
public:
  static DicomFileMgr* GetInstance();
  ~DicomFileMgr(){};

private:
  DicomFileMgr();

public:
  std::vector<DicomFileStructure*> GetStructFiles() const {
    return theStructFiles;
  }

  void SetCompression( G4String fComp );
  void AddFile( G4String fComp );
  void AddMaterial( std::vector<G4String> data );
  void AddMaterialDensity( std::vector<G4String> data );
  void AddCT2Density( std::vector<G4String> data );

  void Convert( G4String fFileName );
  void CheckNColumns(std::vector<G4String> wl, size_t vsizeTh );
  void ProcessFiles();
  void CheckCTSlices();
  G4double Hounsfield2density(Uint32 Hval);
  size_t GetMaterialIndex( G4double Hval );
  size_t GetMaterialIndexByDensity( G4double density );
  void BuildCTMaterials();
  void MergeCTFiles();
  void CheckPETSlices();
  void BuildPETActivities();
  void MergePETFiles();
  void DumpToTextFile();
  void SetStructureNCheck( G4int nsc ){
    theStructureNCheck = nsc;
  }
  G4int GetStructureNCheck() const {
    return theStructureNCheck;
  }
  void SetStructureNMaxROI( G4int nsc ){
    theStructureNMaxROI = nsc;
  }
  G4int GetStructureNMaxROI() const {
    return theStructureNMaxROI;
  }
  G4int GetCompression() const {
    return fCompression;
  }
  G4String GetFileOutName() const {
    return theFileOutName;
  }

  void SetControlPointMetersets();
  G4bool IsMaterialsDensity() const {
    return bMaterialsDensity;
  }
  
protected:
  G4int fCompression;

private:
  static DicomFileMgr* theInstance;

  G4String theFileOutName;
  //  msd theFiles;
  mdct theCTFiles;
  std::vector<DicomFileStructure*> theStructFiles;
  std::vector<DicomFilePlan*> thePlanFiles;
  mdpet thePETFiles;
  std::map<G4double,G4String> theMaterials;
  std::map<G4double,G4String> theMaterialsDensity;
  std::map<G4int,G4double> theCT2Density;

  DicomFileCT* theCTFileAll;
  DicomFilePET* thePETFileAll;
  G4int theStructureNCheck;
  G4int theStructureNMaxROI;

public:
  static int verbose;
  G4bool bMaterialsDensity;
};

#endif
