#ifndef STREAMERHEPREPINSTANCETREE_H
#define STREAMERHEPREPINSTANCETREE_H 1

#include "FreeHepTypes.h"

#include <string>
#include <vector>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepSelectFilter.h"
#include "HEPREP/HepRepInstanceTree.h"
#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepTreeID.h"

#include "DefaultHepRepTreeID.h"

/**
 *
 * @author M.Donszelmann
 * @version $Id: StreamerHepRepInstanceTree.h,v 1.6 2002-11-14 18:36:21 duns Exp $
 */
class StreamerHepRepInstanceTree : public DefaultHepRepTreeID, public virtual HEPREP::HepRepInstanceTree {

    private:
        HEPREP::HepRepWriter* streamer;
        HEPREP::HepRepTreeID* typeTree;

    public:
        StreamerHepRepInstanceTree(HEPREP::HepRepWriter* streamer, std::string name, std::string version, HEPREP::HepRepTreeID* typeTree);
        ~StreamerHepRepInstanceTree();

        HEPREP::HepRepInstanceTree* copy(HEPREP::HepRep* heprep, HEPREP::HepRepSelectFilter* filter);
        HEPREP::HepRepTreeID* copy();
        bool addInstance(HEPREP::HepRepInstance* instance);
        void removeInstance(HEPREP::HepRepInstance* instance);
        std::vector<HEPREP::HepRepInstance*>* getInstances();
        bool addInstanceTree(HEPREP::HepRepTreeID* treeID);
        HEPREP::HepRepTreeID* getTypeTree();
        std::vector<HEPREP::HepRepInstanceTree*>* getInstanceTrees();
};

#endif
