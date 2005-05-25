// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFAULTHEPREPDEFINITION_H
#define CHEPREP_DEFAULTHEPREPDEFINITION_H 1

#include "cheprep/config.h"

#include <string>
#include <vector>
#include <set>

#include "HEPREP/HepRepDefinition.h"
#include "HEPREP/HepRepWriter.h"

#include "DefaultHepRepAttribute.h"

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepDefinition.h,v 1.2 2005-05-25 23:21:59 duns Exp $
 */
namespace cheprep {

class DefaultHepRepDefinition : public DefaultHepRepAttribute, public virtual HEPREP::HepRepDefinition {

    private:
        std::map<std::string, HEPREP::HepRepAttDef*> attDefs;

    public:
        DefaultHepRepDefinition();
        ~DefaultHepRepDefinition();

        void addAttDef(HEPREP::HepRepAttDef* hepRepAttDef);
        void addAttDef(std::string name, std::string desc, std::string type, std::string extra);
        std::set<HEPREP::HepRepAttDef *> getAttDefsFromNode();
        HEPREP::HepRepAttDef* getAttDefFromNode(std::string lowerCaseName);

        HEPREP::HepRepAttDef* getAttDef(std::string name) = 0;
        HEPREP::HepRepAttValue* getAttValue(std::string name) = 0;
};

} // cheprep


#endif
