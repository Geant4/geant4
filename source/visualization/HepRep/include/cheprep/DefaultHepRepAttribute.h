// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFAULTHEPREPATTRIBUTE_H
#define CHEPREP_DEFAULTHEPREPATTRIBUTE_H 1

#include "cheprep/config.h"

#include <string>
#include <map>
#include <set>
#include <vector>

#include "HEPREP/HepRepAttribute.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepConstants.h"
#include "HEPREP/HepRepWriter.h"

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepAttribute.h,v 1.3 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

class DefaultHepRepAttribute : public virtual HEPREP::HepRepAttribute {

    private:
        std::map<std::string, HEPREP::HepRepAttValue*> attValues;

    public:
        DefaultHepRepAttribute();
        ~DefaultHepRepAttribute();

        std::set<HEPREP::HepRepAttValue*> getAttValuesFromNode();
        void addAttValue(HEPREP::HepRepAttValue* hepRepAttValue);
        void addAttValue(std::string key, char *value, int showLabel);
        void addAttValue(std::string key, std::string value, int showLabel);
        void addAttValue(std::string key, int value, int showLabel);
        void addAttValue(std::string key, int64 value, int showLabel);
        void addAttValue(std::string key, double value, int showLabel);
        void addAttValue(std::string key, bool value, int showLabel);
        void addAttValue(std::string key, std::vector<double> value, int showLabel);
        void addAttValue(std::string key, double red, double green, double blue, double alpha, int showLabel);
        HEPREP::HepRepAttValue* getAttValueFromNode(std::string lowerCaseName);
        HEPREP::HepRepAttValue* removeAttValue(std::string key);

        HEPREP::HepRepAttValue* getAttValue(std::string name) = 0;
};

} // cheprep


#endif
