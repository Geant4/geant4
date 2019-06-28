#include "G4HadronicException.hh"

#include <iostream>

G4HadronicException::G4HadronicException(G4String in, G4int at, G4String mess)
    : theMessage{mess}, theName{in}, theLine{at} {

  std::ostringstream os;
  Report(os);
  whatString = os.str();

  G4cout << whatString;

  if (getenv("DumpCoreOnHadronicException")) {
    G4Exception("G4HadronicException", "007", FatalException,
                "Fatal problem in above location");
  }
}

G4HadronicException::~G4HadronicException() throw()
{}

const char* G4HadronicException::what() const noexcept {
  return whatString.c_str();
}

void G4HadronicException::Report(std::ostream &aS) const {
  aS << "In " << theName << ", line " << theLine << ": " << std::endl;
  aS << "===> " << theMessage << std::endl;
}
