
#include <iostream>

#include "StreamerHepRepAttribute.h"
#include "DefaultHepRepAttValue.h"

using namespace std;
using namespace HEPREP;


StreamerHepRepAttribute::StreamerHepRepAttribute(HepRepWriter* stream) {
    this->streamer = stream;
}

StreamerHepRepAttribute::~StreamerHepRepAttribute() {
}

vector<HepRepAttValue*>* StreamerHepRepAttribute::getAttValuesFromNode() {
    return NULL;
}

bool StreamerHepRepAttribute::addAttValue(HepRepAttValue* hepRepAttValue) {
    streamer->write(hepRepAttValue);
    delete hepRepAttValue;
    return true;
}

bool StreamerHepRepAttribute::addAttValue(string key, string value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool StreamerHepRepAttribute::addAttValue(string key, int value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool StreamerHepRepAttribute::addAttValue(string key, double value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool StreamerHepRepAttribute::addAttValue(string key, bool value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool StreamerHepRepAttribute::addAttValue(string key, vector<double> value, int showLabel) {
    return addAttValue(new DefaultHepRepAttValue(key, value, showLabel));
}

bool StreamerHepRepAttribute::addAttValue(string key, double red, double green, double blue, double alpha, int showLabel) {
    vector<double> color;
    color.push_back(red);
    color.push_back(green);
    color.push_back(blue);
    color.push_back(alpha);
    return addAttValue(new DefaultHepRepAttValue(key, color, showLabel));
}

HepRepAttValue* StreamerHepRepAttribute::getAttValueFromNode(string lowerCaseName) {
    return NULL;
}

HepRepAttValue* StreamerHepRepAttribute::removeAttValue(string key) {
    return NULL;
}

HepRepAttValue* StreamerHepRepAttribute::getAttValue(string name) {
    return NULL;
}

