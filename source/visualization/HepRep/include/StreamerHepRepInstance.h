#ifndef STREAMERHEPREPINSTANCE_H
#define STREAMERHEPREPINSTANCE_H 1

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

#include "StreamerHepRepAttribute.h"

/**
 *
 * @author M.Donszelmann
 * @version $Id: StreamerHepRepInstance.h,v 1.7 2002-11-19 21:54:10 duns Exp $
 */
class StreamerHepRepInstance : public StreamerHepRepAttribute, public virtual HEPREP::HepRepInstance {

    private:
        void* parent;
        HEPREP::HepRepType* type;

    public:
        StreamerHepRepInstance(HEPREP::HepRepWriter* streamer, HEPREP::HepRepInstance* parent, HEPREP::HepRepType* type);
        StreamerHepRepInstance(HEPREP::HepRepWriter* streamer, HEPREP::HepRepInstanceTree* parent, HEPREP::HepRepType* type);
        ~StreamerHepRepInstance();

        HEPREP::HepRepInstance* copy(HEPREP::HepRep* heprep, HEPREP::HepRepInstance* parent, HEPREP::HepRepSelectFilter* filter);
        HEPREP::HepRepInstance* copy(HEPREP::HepRep* heprep, HEPREP::HepRepInstanceTree* parent, HEPREP::HepRepSelectFilter* filter);
        HEPREP::HepRepType* getType();
        bool addPoint(HEPREP::HepRepPoint* point);
        std::vector<HEPREP::HepRepPoint *>* getPoints();
        bool addInstance(HEPREP::HepRepInstance* instance);
        void removeInstance(HEPREP::HepRepInstance* instance);
        std::vector<HEPREP::HepRepInstance *>* getInstances();
        HEPREP::HepRepAttValue* getAttValue(std::string name);

        void *getParent() { return parent; }
};

#endif
