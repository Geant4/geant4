// Copyright FreeHEP, 2005.
#ifndef XMLHEPREPFACTORY_H
#define XMLHEPREPFACTORY_H 1

#include <string>
#include <iostream>

#include "HEPREP/HepRepReader.h"
#include "HEPREP/HepRepWriter.h"

#include "DefaultHepRepFactory.h"

/**
 * @author Mark Donszelmann
 * @version $Id: XMLHepRepFactory.h,v 1.2 2005-05-25 23:21:59 duns Exp $
 */
namespace cheprep {

class XMLHepRepFactory : public DefaultHepRepFactory {

    public:
        XMLHepRepFactory();
        ~XMLHepRepFactory();

        HEPREP::HepRepReader* createHepRepReader (std::istream* in);
        HEPREP::HepRepReader* createHepRepReader (std::string filename);
        HEPREP::HepRepWriter* createHepRepWriter (std::ostream* out, bool randomAccess, bool compress);
};

} // cheprep


#endif
