
#include "DefaultHepRepAttribute.h"
#include "DefaultHepRepAttValue.h"

#include <iostream>

using namespace std;
using namespace HEPREP;


DefaultHepRepAttribute::DefaultHepRepAttribute() {
}

DefaultHepRepAttribute::~DefaultHepRepAttribute() {
    vector<HepRepAttValue *>* list = getAttValuesFromNode();
    for (vector<HepRepAttValue*>::iterator i1 = list->begin(); i1 != list->end(); i1++) {
        delete (*i1);
    }
}

vector<HepRepAttValue*>* DefaultHepRepAttribute::getAttValuesFromNode() {
    attList.clear();
    for (map<string, HepRepAttValue*>::iterator i = attValues.begin(); i != attValues.end(); i++) {
        if ((*i).first != "layer") attList.push_back((*i).second);
    }
    return &attList;
}

bool DefaultHepRepAttribute::addAttValue(HepRepAttValue* hepRepAttValue) {
    attValues[hepRepAttValue->getLowerCaseName()] = hepRepAttValue;
    return true;
}

bool DefaultHepRepAttribute::addAttValue(string key, char *value, int showLabel) {
    return addAttValue(key, (std::string)value, showLabel);
}

bool DefaultHepRepAttribute::addAttValue(string key, string value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool DefaultHepRepAttribute::addAttValue(string key, int value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool DefaultHepRepAttribute::addAttValue(string key, double value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool DefaultHepRepAttribute::addAttValue(string key, bool value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool DefaultHepRepAttribute::addAttValue(string key, vector<double> value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool DefaultHepRepAttribute::addAttValue(string key, double red, double green, double blue, double alpha, int showLabel) {
    vector<double> color;
    color.push_back(red);
    color.push_back(green);
    color.push_back(blue);
    color.push_back(alpha);
    return addAttValue(new DefaultHepRepAttValue(key, color, showLabel));
}

HepRepAttValue* DefaultHepRepAttribute::getAttValueFromNode(string name) {
    map<string, HepRepAttValue*>::iterator i = attValues.find(name);
    return i != attValues.end() ? (*i).second : NULL;
}

HepRepAttValue* DefaultHepRepAttribute::removeAttValue(string name) {
    HepRepAttValue* attValue = attValues[name];
    attValues.erase(name);
    return attValue;
}

