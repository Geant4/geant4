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
//
#ifndef STREAMERHEPREPATTRIBUTE_H
#define STREAMERHEPREPATTRIBUTE_H 1

#include "FreeHepTypes.h"

#include <string>
#include <vector>

#include "HEPREP/HepRepAttribute.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepWriter.h"

/**
 *
 * @author M.Donszelmann
 */

class StreamerHepRepAttribute : public virtual HEPREP::HepRepAttribute {

    private:
        HEPREP::HepRepWriter* streamer;

    public:
        StreamerHepRepAttribute(HEPREP::HepRepWriter* streamer);
        ~StreamerHepRepAttribute();

        std::vector<HEPREP::HepRepAttValue*>* getAttValuesFromNode();
        bool addAttValue(HEPREP::HepRepAttValue* hepRepAttValue);
        bool addAttValue(std::string key, std::string value, int showLabel);
        bool addAttValue(std::string key, int value, int showLabel);
        bool addAttValue(std::string key, double value, int showLabel);
        bool addAttValue(std::string key, bool value, int showLabel);
        bool addAttValue(std::string key, std::vector<double> value, int showLabel);
        bool addAttValue(std::string key, double red, double green, double blue, double alpha, int showLabel);
        HEPREP::HepRepAttValue* getAttValueFromNode(std::string lowerCaseName);
        HEPREP::HepRepAttValue* removeAttValue(std::string key);
        HEPREP::HepRepAttValue* getAttValue(std::string name);
};

#endif
