#ifndef DEFAULTHEPREPDEFINITION_H
#define DEFAULTHEPREPDEFINITION_H 1

#include "FreeHepTypes.h"

#include <string>
#include <vector>

#include "HEPREP/HepRepDefinition.h"
#include "HEPREP/HepRepWriter.h"

#include "DefaultHepRepAttribute.h"

/**
 *
 * @author M.Donszelmann
 */
class DefaultHepRepDefinition : public DefaultHepRepAttribute, public virtual HEPREP::HepRepDefinition {

    private:
        std::map<std::string, HEPREP::HepRepAttDef*> attDefs;
        std::vector<HEPREP::HepRepAttDef*> attList;

    public:
        DefaultHepRepDefinition();
        ~DefaultHepRepDefinition();

        bool addAttDef(HEPREP::HepRepAttDef* hepRepAttDef);
        bool addAttDef(std::string name, std::string desc, std::string type, std::string extra);
        std::vector<HEPREP::HepRepAttDef *>* getAttDefsFromNode();
        HEPREP::HepRepAttDef* getAttDefFromNode(std::string lowerCaseName);

        HEPREP::HepRepAttDef* getAttDef(std::string name) = 0;
        HEPREP::HepRepAttValue* getAttValue(std::string name) = 0;
};

#endif
