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
#include "DicomVFile.hh"

#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcdeftag.h"
#include "dcmtk/dcmdata/dcpixel.h"
#include "dcmtk/dcmdata/dcpxitem.h"
#include "dcmtk/dcmdata/dcpixseq.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomVFile::DicomVFile(DcmDataset* dset) : theDataset(dset)
{
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4double> DicomVFile::Read1Data( DcmDataset * dset, DcmTagKey tagKey, G4int nData )
{
  std::vector<G4double> dataV;
    
  for(int ii=0; ii<nData; ++ii) {
    G4double data;
    Uint16 datai;
    // see  http://support.dcmtk.org/docs/classDcmItem.html for types
    if (dset->findAndGetFloat64(tagKey, data,ii).good() ) {
      dataV.push_back(data);
    } else if (dset->findAndGetUint16(tagKey, datai,ii).good() ) {
      dataV.push_back(datai);
    } else {
      G4cout <<"ERROR  (" << std::showbase // show the 0x prefix
             << std::internal // fill between the prefix and the number
             << std::setfill('0') << std::hex << std::setw(4) << tagKey.getGroup() 
             << "," << tagKey.getElement() << ") "<< std::dec << ii << std::endl; 
      G4Exception("DicomHandler::ReadData",
                  "",
                  JustWarning,
                  (std::to_string(data) +G4String(" Have not read (")
                   + std::to_string(tagKey.getGroup())+","+std::to_string(tagKey.getElement())
                   +")"+" : "+std::to_string(ii)).c_str());
    }
  }
  
  return dataV;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
OFString DicomVFile::Read1DataStr( DcmDataset * dset, DcmTagKey tagKey )
{
  //  const char* data = "";
  OFString data;  
  // see  http://support.dcmtk.org/docs/classDcmItem.html for types
  if (dset->findAndGetOFString(tagKey, data).good() ) {
  } else {
    G4cout <<"ERROR  (" << std::showbase // show the 0x prefix
           << std::internal // fill between the prefix and the number
           << std::setfill('0') << std::hex << std::setw(4) << tagKey.getGroup() << "," 
           << tagKey.getElement() << ") "<< std::dec << std::endl; 
    G4Exception("DicomHandler::ReadData",
                "",
                JustWarning,
                (" Have not read (" + std::to_string(tagKey.getGroup())+","
                 +std::to_string(tagKey.getElement())+")"+" : ").c_str());
  }

  return data.c_str();
}
