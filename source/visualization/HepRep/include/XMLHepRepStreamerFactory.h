#ifndef XMLHEPREPSTREAMERFACTORY_H
#define XMLHEPREPSTREAMERFACTORY_H 1

#include "FreeHepTypes.h"

#include <string>
#include <iostream>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepFactory.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepPoint.h"
#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepInstanceTree.h"
#include "HEPREP/HepRepTreeID.h"
#include "HEPREP/HepRepAction.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepTypeTree.h"

/**
 *
 * @author M.Donszelmann
 * @version $Id: XMLHepRepStreamerFactory.h,v 1.6 2002-11-14 18:36:34 duns Exp $
 */
class XMLHepRepStreamerFactory : public virtual HEPREP::HepRepFactory {

    private:
        HEPREP::HepRepWriter* streamer;

    public:
        XMLHepRepStreamerFactory();
        ~XMLHepRepStreamerFactory();

//        static HEPREP::HepRepFactory* create();
        HEPREP::HepRepWriter* createHepRepWriter (std::ostream* out);
        HEPREP::HepRepPoint* createHepRepPoint (HEPREP::HepRepInstance* instance,
                                   double x, double y, double z);
        HEPREP::HepRepInstance* createHepRepInstance (HEPREP::HepRepInstance* parent, HEPREP::HepRepType* type);
        HEPREP::HepRepInstance* createHepRepInstance (HEPREP::HepRepInstanceTree* parent, HEPREP::HepRepType* type);
        HEPREP::HepRepTreeID* createHepRepTreeID (std::string name, std::string version, std::string qualifier = "top-level");
        HEPREP::HepRepAction* createHepRepAction (std::string name, std::string expression);
        HEPREP::HepRepInstanceTree* createHepRepInstanceTree (std::string name, std::string version,
                                                        HEPREP::HepRepTreeID* typeTreeID);
        HEPREP::HepRepType* createHepRepType (HEPREP::HepRepType* parent, std::string name);
        HEPREP::HepRepTypeTree* createHepRepTypeTree (HEPREP::HepRepTreeID* treeID);
        HEPREP::HepRep* createHepRep ();
};

#endif
