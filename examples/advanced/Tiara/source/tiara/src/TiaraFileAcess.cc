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
