//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4XXXFileViewer.hh,v 1.1 2006-03-28 17:16:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    FileWriter(const G4String& fileName): fFileName(fileName), fOpen(false) {}
    G4bool IsOpen() {return fOpen;}
    void WriteItem(const G4String& item)
    {
      if (!fOpen) {fFile.open(fFileName); fOpen = true;}
      if (fFile.good()) fFile << item << std::endl;
      else G4cout << "G4XXXFileViewer::FileWriter::WriteItem: ERROR" << G4endl;
    }
    void Close() {if (fOpen) {fFile.close(); fOpen = false;}}
    void Rewind() {if (fOpen) fFile.seekp(0, std::ios::beg);}
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
