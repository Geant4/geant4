#ifndef STREAMERHEPREPTYPE_H
#define STREAMERHEPREPTYPE_H 1

#include "FreeHepTypes.h"

#include <string>
#include <vector>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepAttDef.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepWriter.h"

#include "StreamerHepRepDefinition.h"

/**
 *
 * @author M.Donszelmann
 * @version $Id: StreamerHepRepType.h,v 1.2 2002-11-13 18:39:18 duns Exp $
 */
class StreamerHepRepType : public StreamerHepRepDefinition, public virtual HEPREP::HepRepType {

    private:
        HEPREP::HepRepType* parent;
        std::string name;
        std::string description;
        std::string infoURL;

    public:
        StreamerHepRepType(HEPREP::HepRepWriter* streamer, HEPREP::HepRepType* parent, std::string name);
        ~StreamerHepRepType();

        HEPREP::HepRepType* getSuperType();
        HEPREP::HepRepAttDef* getAttDef(std::string name);
        HEPREP::HepRepAttValue* getAttValue(std::string name);
        HEPREP::HepRepType* copy(HEPREP::HepRep* heprep, HEPREP::HepRepType* parent);
        std::string getName();
        std::string getDescription();
        void setDescription(std::string description);
        std::string getInfoURL();
        void setInfoURL(std::string infoURL);
        bool addType(HEPREP::HepRepType* type);
        std::vector<HEPREP::HepRepType*>* getTypes();
};

#endif
