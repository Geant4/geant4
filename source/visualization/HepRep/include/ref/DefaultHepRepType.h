#ifndef DEFAULTHEPREPTYPE_H
#define DEFAULTHEPREPTYPE_H 1

#include "FreeHepTypes.h"

#include <string>
#include <vector>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepAttDef.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepWriter.h"

#include "DefaultHepRepDefinition.h"

/**
 *
 * @author M.Donszelmann
 */
class DefaultHepRepType : public DefaultHepRepDefinition, public virtual HEPREP::HepRepType {

    private:
        HEPREP::HepRepType* parent;
        std::vector<HEPREP::HepRepType*> types;
        std::string name;
        std::string description;
        std::string infoURL;

    public:
        DefaultHepRepType(HEPREP::HepRepType* parent, std::string name);
        ~DefaultHepRepType();

        HEPREP::HepRepType* getSuperType();
        HEPREP::HepRepAttDef* getAttDef(std::string name);
        HEPREP::HepRepAttValue* getAttValue(std::string name);
        HEPREP::HepRepType* copy(HEPREP::HepRep* heprep, HEPREP::HepRepType* parent);
        std::string getName();
        std::string getFullName();
        std::string getDescription();
        void setDescription(std::string description);
        std::string getInfoURL();
        void setInfoURL(std::string infoURL);
        bool addType(HEPREP::HepRepType* type);
        std::vector<HEPREP::HepRepType*>* getTypes();
};

#endif
