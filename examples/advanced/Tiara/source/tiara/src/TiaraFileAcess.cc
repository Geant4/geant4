// $Id: TiaraFileAcess.cc,v 1.2 2003-06-16 17:06:48 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "TiaraFileAcess.hh"
#include <fstream>


void checkFileIsReadable(const G4String &fileName,
			 const G4String &caller) {
  std::ifstream file(fileName);
  if (! file) {
    G4Exception("checkFile: couldn't open file: " + fileName + "\n  called from: " + caller);
  }
  return;
}
