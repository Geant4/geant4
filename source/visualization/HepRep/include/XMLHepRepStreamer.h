#ifndef XMLHEPREPSTREAMER_H
#define XMLHEPREPSTREAMER_H 1

#include "FreeHepTypes.h"

#include <iostream>
#include <string>
#include <stack>

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
#include "HEPREP/HepRepAttDef.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepAttribute.h"

#include "XMLWriter.h"

/**
 *
 * @author M.Donszelmann
 * @version $Id: XMLHepRepStreamer.h,v 1.2 2002-11-13 18:39:26 duns Exp $
 */
class XMLHepRepStreamer : public XMLWriter, public virtual HEPREP::HepRepWriter {

    private:
        std::string nameSpace;
        bool writtenLayers;
        std::stack<void *> hepreps;
        HEPREP::HepRep *heprep;
        HEPREP::HepRepTypeTree *heprepTypeTree;

    public:
        XMLHepRepStreamer(std::ostream* out);
        ~XMLHepRepStreamer();

        bool close();
        bool write(HEPREP::HepRep* heprep);
        bool write(HEPREP::HepRepTypeTree* typeTree);
        bool write(HEPREP::HepRepType* type);
        bool write(HEPREP::HepRepTreeID* treeID);
        bool write(HEPREP::HepRepAction* action);
        bool write(HEPREP::HepRepInstanceTree* instanceTree);
        bool write(HEPREP::HepRepInstance* instance);
        bool write(HEPREP::HepRepPoint* point);
        bool write(HEPREP::HepRepAttribute* attribute);
        bool write(HEPREP::HepRepDefinition* definition);
        bool write(HEPREP::HepRepAttValue* attValue);
        bool write(HEPREP::HepRepAttDef* attDef);
};

#endif
