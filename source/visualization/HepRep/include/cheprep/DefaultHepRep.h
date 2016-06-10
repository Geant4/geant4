// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFAULTHEPREP_H
#define CHEPREP_DEFAULTHEPREP_H 1

#include "cheprep/config.h"

#include <string>
#include <vector>
#include <set>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepSelectFilter.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepTypeTree.h"
#include "HEPREP/HepRepInstanceTree.h"

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRep.h 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

class DefaultHepRep : public virtual HEPREP::HepRep {

    private:
        std::vector<std::string> layers;
        std::vector<HEPREP::HepRepTypeTree*> typeTrees;
        std::vector<HEPREP::HepRepInstanceTree*> instanceTrees;

    public:
        DefaultHepRep();
        ~DefaultHepRep();

        HEPREP::HepRep* copy(HEPREP::HepRepSelectFilter* filter);
        std::vector<std::string> getLayerOrder();
        void addLayer(std::string layer);
        void addTypeTree(HEPREP::HepRepTypeTree* typeTree);
        void removeTypeTree(HEPREP::HepRepTypeTree* typeTree);
        HEPREP::HepRepTypeTree* getTypeTree(std::string name, std::string version);
        std::vector<HEPREP::HepRepTypeTree*> getTypeTreeList();
        void addInstanceTree(HEPREP::HepRepInstanceTree* instanceTree);
        void overlayInstanceTree(HEPREP::HepRepInstanceTree * instanceTree);
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
        std::vector<HEPREP::HepRepInstanceTree*> getInstanceTreeList();
};

} // cheprep

#endif
