// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFAULTHEPREPATTVALUE_H
#define CHEPREP_DEFAULTHEPREPATTVALUE_H 1

#include "cheprep/config.h"

#include <string>

#include "HEPREP/HepRepAttValue.h"

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepAttValue.h,v 1.3 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

class DefaultHepRepAttValue : public virtual HEPREP::HepRepAttValue {

    private:
        enum { LABELSTRINGS_LEN = 4 };
        std::string name;
        int type;

        // values implemented as separate items, so that they do not take up unnecessary space for an Object
        // only ONE of these is filled
        std::string stringValue;
        int64 longValue;
        double doubleValue;
        bool booleanValue;
        std::vector<double> colorValue;

        int showLabelValue;
        static std::string labelStrings[LABELSTRINGS_LEN];

        void init();

    public:
        DefaultHepRepAttValue(std::string name, std::string value, int showLabel);
        DefaultHepRepAttValue(std::string name, int64 value, int showLabel);
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
        std::string getLowerCaseString();
        int64 getLong();
        int getInteger();
        double getDouble();
        bool getBoolean();
        std::vector<double> getColor();

        std::string getAsString();
        static std::string getAsString(std::vector<double> c);
        static std::string getAsString(int i);
        static std::string getAsString(int64 i);
        static std::string getAsString(double d);
        static std::string getAsString(bool b);
        
        std::string toShowLabel();
        static std::string toShowLabel(int showLabel);
};

} // cheprep


#endif

