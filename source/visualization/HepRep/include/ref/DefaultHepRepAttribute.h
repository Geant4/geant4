#ifndef DEFAULTHEPREPATTRIBUTE_H
#define DEFAULTHEPREPATTRIBUTE_H 1

#include "FreeHepTypes.h"

#include <string>
#include <map>
#include <vector>

#include "HEPREP/HepRepAttribute.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepConstants.h"
#include "HEPREP/HepRepWriter.h"

/**
 *
 * @author M.Donszelmann
 */
class DefaultHepRepAttribute : public virtual HEPREP::HepRepAttribute {

    private:
        std::map<std::string, HEPREP::HepRepAttValue*> attValues;
        std::vector<HEPREP::HepRepAttValue*> attList;

    public:
        DefaultHepRepAttribute();
        ~DefaultHepRepAttribute();

        std::vector<HEPREP::HepRepAttValue*>* getAttValuesFromNode();
        bool addAttValue(HEPREP::HepRepAttValue* hepRepAttValue);
        bool addAttValue(std::string key, char *value, int showLabel);
        bool addAttValue(std::string key, std::string value, int showLabel);
        bool addAttValue(std::string key, int value, int showLabel);
        bool addAttValue(std::string key, double value, int showLabel);
        bool addAttValue(std::string key, bool value, int showLabel);
        bool addAttValue(std::string key, std::vector<double> value, int showLabel);
        bool addAttValue(std::string key, double red, double green, double blue, double alpha, int showLabel);
        HEPREP::HepRepAttValue* getAttValueFromNode(std::string lowerCaseName);
        HEPREP::HepRepAttValue* removeAttValue(std::string key);

        HEPREP::HepRepAttValue* getAttValue(std::string name) = 0;
};

#endif
