#ifndef DEFAULTHEPREPATTDEF_H
#define DEFAULTHEPREPATTDEF_H 1

#include "FreeHepTypes.h"

#include <string>

#include "HEPREP/HepRepAttDef.h"

/**
 *
 * @author M.Donszelmann
 * @version $Id: DefaultHepRepAttDef.h,v 1.1 2002-11-13 07:05:45 duns Exp $
 */

class DefaultHepRepAttDef : public virtual HEPREP::HepRepAttDef {

    private:
        std::string name, desc, category, extra;

    public:
        DefaultHepRepAttDef(std::string name, std::string desc, std::string category, std::string extra);
        ~DefaultHepRepAttDef();

        HEPREP::HepRepAttDef* copy();
        std::string getName();
        std::string getLowerCaseName();
        std::string getDescription();
        std::string getCategory();
        std::string getExtra();
};

#endif
