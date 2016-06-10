// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFAULTHEPREPTYPE_H
#define CHEPREP_DEFAULTHEPREPTYPE_H 1

#include "cheprep/config.h"

#include <string>
#include <vector>
#include <set>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepAttDef.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepWriter.h"

#include "DefaultHepRepDefinition.h"

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepType.h 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

class DefaultHepRepType : public DefaultHepRepDefinition, public virtual HEPREP::HepRepType {

    private:
        HEPREP::HepRepType* parent;
        std::vector<HEPREP::HepRepType*> types;
        std::string name;
        std::string description;
        std::string infoURL;

    public:
        DefaultHepRepType(HEPREP::HepRepType* parent, std::string name);
        DefaultHepRepType(HEPREP::HepRepTypeTree* parent, std::string name);
        ~DefaultHepRepType();

        HEPREP::HepRepType* getSuperType();
        HEPREP::HepRepAttDef* getAttDef(std::string name);
        HEPREP::HepRepAttValue* getAttValue(std::string name);
        HEPREP::HepRepType* copy(HEPREP::HepRepType* parent);
        std::string getName();
        std::string getFullName();
        std::string getDescription();
        void setDescription(std::string description);
        std::string getInfoURL();
        void setInfoURL(std::string infoURL);
        void addType(HEPREP::HepRepType* type);
        std::vector<HEPREP::HepRepType*> getTypeList();
};

} // cheprep


#endif
