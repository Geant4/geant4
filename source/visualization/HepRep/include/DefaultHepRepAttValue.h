#ifndef DEFAULTHEPREPATTVALUE_H
#define DEFAULTHEPREPATTVALUE_H 1

#include "FreeHepTypes.h"

#include <string>

#include "HEPREP/HepRepAttValue.h"

/**
 *
 * @author M.Donszelmann
 * @version $Id: DefaultHepRepAttValue.h,v 1.7 2002-11-19 21:53:37 duns Exp $
 */

class DefaultHepRepAttValue : public virtual HEPREP::HepRepAttValue {

    private:
        enum { LABELSTRINGS_LEN = 4 };
        std::string name;
        int type;

        // values implemented as separate items, so that they do not take up unnecessary space for an Object
        // only ONE of these is filled
        std::string stringValue;
        long longValue;
        double doubleValue;
        bool booleanValue;
        std::vector<double> colorValue;

        int showLabelValue;
        std::string labelStrings[LABELSTRINGS_LEN];

        void init();

    public:
        DefaultHepRepAttValue(std::string name, std::string value, int showLabel);
        DefaultHepRepAttValue(std::string name, long value, int showLabel);
        DefaultHepRepAttValue(std::string name, int value, int showLabel);
        DefaultHepRepAttValue(std::string name, double value, int showLabel);
        DefaultHepRepAttValue(std::string name, bool value, int showLabel);
        DefaultHepRepAttValue(std::string name, std::vector<double> value, int showLabel);
        ~DefaultHepRepAttValue();

        HEPREP::HepRepAttValue* copy();

        std::string getName();
        std::string getLowerCaseName();
        int getType();
        std::string getTypeName();
        int showLabel();
        std::string getString();
        long getLong();
        int getInteger();
        double getDouble();
        bool getBoolean();
        std::vector<double> getColor();
        std::string getAsString();

        std::string toShowLabel();
};

#endif

