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

#include "DefaultHepRepAttDef.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepAttDef::DefaultHepRepAttDef(string name, string desc, string category, string extra)
    : name(name), desc(desc), category(category), extra(extra) {
}

DefaultHepRepAttDef::~DefaultHepRepAttDef() {
}

HepRepAttDef* DefaultHepRepAttDef::copy() {
    return NULL;
}

string DefaultHepRepAttDef::getName() {
    return name;
}

string DefaultHepRepAttDef::getLowerCaseName() {
    char* tmp = new char[strlen(name.c_str())];
    strcpy(tmp, name.c_str());
    int i = -1;
    do {
        i++;
        tmp[i] = tolower(tmp[i]);
    } while (tmp[i] != 0);
    return tmp;
}

string DefaultHepRepAttDef::getDescription() {
    return desc;
}

string DefaultHepRepAttDef::getCategory() {
    return category;
}

string DefaultHepRepAttDef::getExtra() {
    return extra;
}

