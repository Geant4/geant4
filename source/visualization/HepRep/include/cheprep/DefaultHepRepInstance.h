// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFAULTHEPREPINSTANCE_H
#define CHEPREP_DEFAULTHEPREPINSTANCE_H 1

#include "cheprep/config.h"

#include <string>
#include <vector>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepSelectFilter.h"
#include "HEPREP/HepRepInstanceTree.h"
#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepPoint.h"
#include "HEPREP/HepRepAttValue.h"

#include "DefaultHepRepAttribute.h"

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepInstance.h 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

class DefaultHepRepInstance : public DefaultHepRepAttribute, public virtual HEPREP::HepRepInstance {

    private:
        HEPREP::HepRepInstance* parent;
        HEPREP::HepRepType* type;
        std::vector<HEPREP::HepRepPoint*> points;
        std::vector<HEPREP::HepRepInstance*> instances;

    public:
        DefaultHepRepInstance(HEPREP::HepRepInstance* parent, HEPREP::HepRepType* type);
        DefaultHepRepInstance(HEPREP::HepRepInstanceTree* parent, HEPREP::HepRepType* type);
        ~DefaultHepRepInstance();

        void overlay(HEPREP::HepRepInstance * instance);
        HEPREP::HepRepInstance* copy(HEPREP::HepRepTypeTree* typeTree, HEPREP::HepRepInstance* parent, HEPREP::HepRepSelectFilter* filter);
        HEPREP::HepRepInstance* copy(HEPREP::HepRepTypeTree* typeTree, HEPREP::HepRepInstanceTree* parent, HEPREP::HepRepSelectFilter* filter);
        HEPREP::HepRepType* getType();
        void addPoint(HEPREP::HepRepPoint* point);
        std::vector<HEPREP::HepRepPoint *> getPoints();
        HEPREP::HepRepInstance* getSuperInstance();
        void addInstance(HEPREP::HepRepInstance* instance);
        void removeInstance(HEPREP::HepRepInstance* instance);
        std::vector<HEPREP::HepRepInstance *> getInstances();

        HEPREP::HepRepAttValue* getAttValue(std::string name);

        void *getParent() { return parent; }
};

} // cheprep


#endif
