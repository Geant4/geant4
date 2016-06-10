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
// $Id: G4XXXFileViewer.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// John Allison  7th March 2006
// A template for a file-writing graphics driver.
//?? Lines beginning like this require specialisation for your driver.

#ifndef G4XXXFileVIEWER_HH
#define G4XXXFileVIEWER_HH

#include "G4VViewer.hh"

#include <fstream>

class G4XXXFileViewer: public G4VViewer {
public:
  G4XXXFileViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4XXXFileViewer();
  void SetView();
  void ClearView();
  void DrawView();
  void ShowView();
  // A simple class to handle delayed opening, etc.  There are various
  // degrees of sophistication in, for example, the allocation of
  // filenames -- see FukuiRenderer or HepRepFile.
  class FileWriter {
  public:
    FileWriter(): fOpen(false) {}
    G4bool IsOpen() {return fOpen;}
    void WriteItem(const G4String& item);
    void Close();
    // Implement rewind as close and re-open...
    void Rewind() {if (fOpen) {fFile.close(); fFile.open(fFileName.c_str());}}
  private:
    G4String fFileName;
    G4bool fOpen;
    std::ofstream fFile;
  };
  FileWriter& GetFileWriter() {return fFileWriter;}
private:
  FileWriter fFileWriter;
};

#endif
