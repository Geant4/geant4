//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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

