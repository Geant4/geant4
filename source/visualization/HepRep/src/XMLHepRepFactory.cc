
#include <iostream>
#include <fstream>

#include "XMLHepRepWriter.h"
#include "XMLHepRepFactory.h"

using namespace std;
using namespace HEPREP;


XMLHepRepFactory::XMLHepRepFactory() {
}

XMLHepRepFactory::~XMLHepRepFactory() {
}

HepRepReader* XMLHepRepFactory::createHepRepReader (istream*) {
    cerr << "XMLHepRepFactory::createHepRepReader not implemented" << endl;
    return NULL;
}

HepRepReader* XMLHepRepFactory::createHepRepReader (std::string) {
    cerr << "XMLHepRepFactory::createHepRepReader not implemented" << endl;
    return NULL;
}

HepRepWriter* XMLHepRepFactory::createHepRepWriter(ostream* out, bool randomAccess, bool compress) {
    return new XMLHepRepWriter(out, randomAccess, compress);
}
