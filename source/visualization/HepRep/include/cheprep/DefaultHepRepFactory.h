// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFAULTHEPREPFACTORY_H
#define CHEPREP_DEFAULTHEPREPFACTORY_H 1

#include "cheprep/config.h"

#include <string>
#include <iostream>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepFactory.h"
#include "HEPREP/HepRepReader.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepPoint.h"
#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepInstanceTree.h"
#include "HEPREP/HepRepTreeID.h"
#include "HEPREP/HepRepAction.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepTypeTree.h"

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepFactory.h,v 1.3 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

class DefaultHepRepFactory : public virtual HEPREP::HepRepFactory {

    public:
        DefaultHepRepFactory();
        ~DefaultHepRepFactory();

//        static HEPREP::HepRepFactory* create();
        HEPREP::HepRepReader* createHepRepReader (std::istream* in);
        HEPREP::HepRepReader* createHepRepReader (std::string filename);
        HEPREP::HepRepWriter* createHepRepWriter (std::ostream* out, bool randomAccess, bool compress);
        HEPREP::HepRepPoint* createHepRepPoint (HEPREP::HepRepInstance* instance,
                                   double x, double y, double z);
        HEPREP::HepRepInstance* createHepRepInstance (HEPREP::HepRepInstance* parent, HEPREP::HepRepType* type);
        HEPREP::HepRepInstance* createHepRepInstance (HEPREP::HepRepInstanceTree* parent, HEPREP::HepRepType* type);
        HEPREP::HepRepTreeID* createHepRepTreeID (std::string name, std::string version, std::string qualifier = "top-level");
        HEPREP::HepRepAction* createHepRepAction (std::string name, std::string expression);
        HEPREP::HepRepInstanceTree* createHepRepInstanceTree (std::string name, std::string version,
                                                        HEPREP::HepRepTreeID* typeTreeID);
        HEPREP::HepRepType* createHepRepType (HEPREP::HepRepType* parent, std::string name);
        HEPREP::HepRepType* createHepRepType (HEPREP::HepRepTypeTree* parent, std::string name);
        HEPREP::HepRepTypeTree* createHepRepTypeTree (HEPREP::HepRepTreeID* treeID);
        HEPREP::HepRep* createHepRep ();
};

} // cheprep


#endif
