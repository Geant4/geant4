// Copyright FreeHEP, 2005.

#include <iostream>
#include <fstream>

#include "cheprep/XMLHepRepWriter.h"
#include "cheprep/XMLHepRepFactory.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: XMLHepRepFactory.cc,v 1.5 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {


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

} // cheprep

