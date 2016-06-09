#ifndef XMLHEPREPWRITER_H
#define XMLHEPREPWRITER_H 1

#include "FreeHepTypes.h"

#include <iostream>
#include <string>
#include <stack>
#include <set>
#include <vector>
#include <map>

#include "zipios++/zipoutputstream.h"
#include "zipios++/gzipoutputstream.h"

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
 */
class XMLHepRepWriter : public virtual HEPREP::HepRepWriter {

    private:
        std::string nameSpace;
        XMLWriter *xml;
        zipios::ZipOutputStream *zip;
        zipios::GZIPOutputStream *gz;
        std::map<std::string, std::string> properties;
        
    public:
        XMLHepRepWriter(std::ostream* out, bool randomAccess, bool compress);
        ~XMLHepRepWriter();

        bool addProperty(std::string key, std::string value);
        bool close();
        bool write(HEPREP::HepRep* heprep, std::string name);
        bool write(std::vector<std::string> layers);
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
