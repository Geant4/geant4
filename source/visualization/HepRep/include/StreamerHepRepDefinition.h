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
 * @version $Id: StreamerHepRepDefinition.h,v 1.1 2002-11-13 07:06:22 duns Exp $
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
