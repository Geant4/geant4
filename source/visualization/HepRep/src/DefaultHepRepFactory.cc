// Copyright FreeHEP, 2005.

#include <iostream>
#include <fstream>

#include "cheprep/DefaultHepRepFactory.h"
#include "cheprep/DefaultHepRepPoint.h"
#include "cheprep/DefaultHepRepInstance.h"
#include "cheprep/DefaultHepRepInstanceTree.h"
#include "cheprep/DefaultHepRepType.h"
#include "cheprep/DefaultHepRepTypeTree.h"
#include "cheprep/DefaultHepRep.h"
#include "cheprep/DefaultHepRepAction.h"
#include "cheprep/DefaultHepRepTreeID.h"

using namespace std;
using namespace HEPREP;


/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepFactory.cc 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

DefaultHepRepFactory::DefaultHepRepFactory() {
}

DefaultHepRepFactory::~DefaultHepRepFactory() {
}

HepRepReader* DefaultHepRepFactory::createHepRepReader (istream*) {
    cerr << "DefaultHepRepFactory::createHepRepReader not implemented" << endl;
    return NULL;
}

HepRepReader* DefaultHepRepFactory::createHepRepReader (std::string) {
    cerr << "DefaultHepRepFactory::createHepRepReader not implemented" << endl;
    return NULL;
}

HepRepWriter* DefaultHepRepFactory::createHepRepWriter(ostream*, bool, bool) {
    cerr << "DefaultHepRepFactory::createHepRepWriter not implemented" << endl;
    return NULL;
}

HepRepPoint* DefaultHepRepFactory::createHepRepPoint (HepRepInstance* instance,
                               double x, double y, double z) {
    return new DefaultHepRepPoint(instance, x, y, z);
}

HepRepInstance* DefaultHepRepFactory::createHepRepInstance (HepRepInstance* parent, HepRepType* type) {
    return new DefaultHepRepInstance(parent, type);
}

HepRepInstance* DefaultHepRepFactory::createHepRepInstance (HepRepInstanceTree* parent, HepRepType* type) {
    return new DefaultHepRepInstance(parent, type);
}

HepRepTreeID* DefaultHepRepFactory::createHepRepTreeID (string name, string version, string qualifier) {
    return new DefaultHepRepTreeID(name, version, qualifier);
}

HepRepAction* DefaultHepRepFactory::createHepRepAction (string name, string expression) {
    return new DefaultHepRepAction(name, expression);
}

HepRepInstanceTree* DefaultHepRepFactory::createHepRepInstanceTree (string name, string version,
                                                    HepRepTreeID* typeTreeID) {
    return new DefaultHepRepInstanceTree(name, version, typeTreeID);
}

HepRepType* DefaultHepRepFactory::createHepRepType (HepRepType* parent, string name) {
    return new DefaultHepRepType(parent, name);
}

HepRepType* DefaultHepRepFactory::createHepRepType (HepRepTypeTree* parent, string name) {
    return new DefaultHepRepType(parent, name);
}

HepRepTypeTree* DefaultHepRepFactory::createHepRepTypeTree (HepRepTreeID* treeID) {
    return new DefaultHepRepTypeTree(treeID);
}

HepRep* DefaultHepRepFactory::createHepRep () {
    return new DefaultHepRep();
}

} // cheprep

