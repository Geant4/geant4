// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFAULTHEPREPTYPETREE_H
#define CHEPREP_DEFAULTHEPREPTYPETREE_H 1

#include "cheprep/config.h"

#include <string>
#include <vector>
#include <set>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepTypeTree.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepTreeID.h"

#include "DefaultHepRepTreeID.h"

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepTypeTree.h,v 1.3 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

class DefaultHepRepTypeTree : public DefaultHepRepTreeID, public virtual HEPREP::HepRepTypeTree {

    private:
        std::vector<HEPREP::HepRepType*> types;

    public:
        DefaultHepRepTypeTree(HEPREP::HepRepTreeID* typeTree);
        ~DefaultHepRepTypeTree();

        HEPREP::HepRepTypeTree* copy();
        void addType(HEPREP::HepRepType* type);
        std::vector<HEPREP::HepRepType* > getTypeList();
        HEPREP::HepRepType* getType(std::string name);
};

} // cheprep


#endif
