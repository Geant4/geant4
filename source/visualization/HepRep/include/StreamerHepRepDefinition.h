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
#ifndef STREAMERHEPREPDEFINITION_H
#define STREAMERHEPREPDEFINITION_H 1

#include "FreeHepTypes.h"

#include <string>
#include <vector>

#include "HEPREP/HepRepDefinition.h"
#include "HEPREP/HepRepWriter.h"

#include "StreamerHepRepAttribute.h"

/**
 *
 * @author M.Donszelmann
 */
class StreamerHepRepDefinition : public StreamerHepRepAttribute, public virtual HEPREP::HepRepDefinition {

    private:
        HEPREP::HepRepWriter* streamer;

    public:
        StreamerHepRepDefinition(HEPREP::HepRepWriter* streamer);
        ~StreamerHepRepDefinition();

        bool addAttDef(HEPREP::HepRepAttDef* hepRepAttDef);
        bool addAttDef(std::string name, std::string desc, std::string type, std::string extra);
        HEPREP::HepRepAttDef* getAttDef(std::string name);
        std::vector<HEPREP::HepRepAttDef *>* getAttDefsFromNode();
        HEPREP::HepRepAttDef* getAttDefFromNode(std::string lowerCaseName);
};

#endif
