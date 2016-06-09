
#include "DefaultHepRepType.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepType::DefaultHepRepType(HepRepType* parentType, string typeName)
    : DefaultHepRepDefinition(), parent(parentType), name(typeName) {
    this->description = "No Description";
    this->infoURL = "No Info URL";

    // HepRepTypes are sometimes used without a parent (top-level)
    if (parent != NULL) {
        parent->addType(this);
    }
}

DefaultHepRepType::~DefaultHepRepType() {
    for (vector<HepRepType*>::iterator i1 = types.begin(); i1 != types.end(); i1++) {
        delete (*i1);
    }
}

HepRepType* DefaultHepRepType::getSuperType() {
    return parent;
}

HepRepAttDef* DefaultHepRepType::getAttDef(string defName) {
    HepRepAttDef* def = NULL;
    HepRepType* type = this;
    while ((def == NULL) && (type != NULL)) {
        def = type->getAttDefFromNode(defName);
        type = type->getSuperType();
    }
    if (def == NULL) {
        cerr << "ERROR: No HepRepDefaults, trying to get definition for: " << defName << endl;
        // FIXME, no HepRepDefaults
    }
    return def;
}

/**
 * searched for a value with given name. Search up the type tree if needed.
 */
HepRepAttValue* DefaultHepRepType::getAttValue(string attName) {
    HepRepAttValue* value = NULL;
    HepRepType* type = this;
    while ((value == NULL) && (type != NULL)) {
        value = type->getAttValueFromNode(attName);
        type = type->getSuperType();
    }
    if (value == NULL) {
        cerr << "ERROR: No HepRepDefaults, trying to get value for: " << attName << endl;
        // FIXME, no HepRepDefaults
    }
    return value;
}

HepRepType* DefaultHepRepType::copy(HepRep*, HepRepType*) {
    cerr << "DefaultHepRepType::copy(HepRep*, HepRepType*) not implemented." << endl;
    return NULL;
}

string DefaultHepRepType::getName() {
    return name;
}

string DefaultHepRepType::getFullName() {
    return (getSuperType() == NULL) ? getName() : getSuperType()->getFullName() + "/" + getName();
}

string DefaultHepRepType::getDescription() {
    return description;
}

void DefaultHepRepType::setDescription(string desc) {
    this->description = desc;
}

string DefaultHepRepType::getInfoURL() {
    return infoURL;
}

void DefaultHepRepType::setInfoURL(string info) {
    this->infoURL = info;
}

bool DefaultHepRepType::addType(HepRepType* type) {
    types.push_back(type);
    return true;
}

vector<HepRepType*>* DefaultHepRepType::getTypes() {
    return &types;
}

