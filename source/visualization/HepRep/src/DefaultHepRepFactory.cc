
#include <iostream>
#include <fstream>

#include "DefaultHepRepFactory.h"
#include "DefaultHepRepPoint.h"
#include "DefaultHepRepInstance.h"
#include "DefaultHepRepInstanceTree.h"
#include "DefaultHepRepType.h"
#include "DefaultHepRepTypeTree.h"
#include "DefaultHepRep.h"
#include "DefaultHepRepAction.h"
#include "DefaultHepRepTreeID.h"

using namespace std;
using namespace HEPREP;


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

