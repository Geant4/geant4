#ifndef STREAMERHEPREP_H
#define STREAMERHEPREP_H 1

#include "FreeHepTypes.h"

#include <string>
#include <vector>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepSelectFilter.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepTypeTree.h"
#include "HEPREP/HepRepInstanceTree.h"

/**
 *
 * @author M.Donszelmann
 * @version $Id: StreamerHepRep.h,v 1.3 2002-11-13 18:49:47 duns Exp $
 */
class StreamerHepRep : public virtual HEPREP::HepRep {

    private:
        std::vector<std::string> layers;

    public:
        StreamerHepRep(HEPREP::HepRepWriter* streamer);
        ~StreamerHepRep();

        HEPREP::HepRep* copy(HEPREP::HepRepSelectFilter* filter);
        std::vector<std::string>* getLayerOrder();
        bool addLayer(std::string layer);
        bool addTypeTree(HEPREP::HepRepTypeTree* typeTree);
        void removeTypeTree(HEPREP::HepRepTypeTree* typeTree);
        HEPREP::HepRepTypeTree* getTypeTree(std::string name, std::string version);
        std::vector<HEPREP::HepRepTypeTree*>* getTypeTrees();
        HEPREP::HepRepType* getType(std::string name);
        bool addInstanceTree(HEPREP::HepRepInstanceTree* instanceTree);
        void removeInstanceTree(HEPREP::HepRepInstanceTree* instanceTree);
        HEPREP::HepRepInstanceTree* getInstanceTreeTop(std::string name, std::string version);
        HEPREP::HepRepInstanceTree* getInstances(std::string name, std::string version,
                                                 std::vector<std::string> typeNames);
        HEPREP::HepRepInstanceTree* getInstancesAfterAction(
                                            std::string name,
                                            std::string version,
                                            std::vector<std::string> typeNames,
                                            std::vector<HEPREP::HepRepAction*> actions,
                                            bool getPoints,
                                            bool getDrawAtts,
                                            bool getNonDrawAtts,
                                            std::vector<std::string> invertAtts);
        std::string checkForException();
        std::vector<HEPREP::HepRepInstanceTree*>* getInstanceTrees();
};

#endif
