#ifndef DEFAULTHEPREPINSTANCE_H
#define DEFAULTHEPREPINSTANCE_H 1

#include "FreeHepTypes.h"

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
 *
 * @author M.Donszelmann
 */
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

#endif
