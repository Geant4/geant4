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
#ifndef DicomFileCT__HH
#define DicomFileCT__HH

#include "DicomVFileImage.hh"

class DicomFileCT : public DicomVFileImage
{ 
public:
  DicomFileCT();
  DicomFileCT(DcmDataset* dset);
  ~DicomFileCT(){};

public:
  void BuildMaterials();
  void DumpMateIDsToTextFile(std::ofstream& fout);
  void DumpDensitiesToTextFile(std::ofstream& fout);
  void BuildStructureIDs();
  void DumpStructureIDsToTextFile(std::ofstream& fout);

private:
  std::vector<size_t> fMateIDs;
  std::vector<G4double> fDensities;
  std::vector<G4int> fStructure;
  //  G4int* fStructure;
  
};

#endif
