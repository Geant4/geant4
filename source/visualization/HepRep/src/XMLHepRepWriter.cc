// Copyright FreeHEP, 2005.

#include "cheprep/XMLHepRepWriter.h"
#include "cheprep/XMLWriter.h"
#include "cheprep/BHepRepWriter.h"

#include "cheprep/DefaultHepRepInstance.h"
#include "cheprep/DefaultHepRepAttValue.h"

#define NAMESPACE "heprep"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: XMLHepRepWriter.cc,v 1.15 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

XMLHepRepWriter::XMLHepRepWriter(ostream* os, bool randomAccess, bool useCompression) 
        : out(os),
          compress(useCompression) {
            
    this->nameSpace = NAMESPACE;

    if (randomAccess) {
        zip = new ZipOutputStream(*os);
        out = zip;
        gz = NULL;
    } else {
        zip = NULL;
        if (useCompression) {
#ifndef CHEPREP_NO_ZLIB
            gz = new GZIPOutputStream(*os);
            out = gz;
#else
            cerr << "WARNING: the .gz output stream you are creating will be a plain file," << endl;
            cerr << "since compression support (ZLIB) was not compiled into the library." << endl;
            cerr << "To add ZLIB support, you need to undefine CHEPREP_NO_ZLIB." << endl;  
            gz = NULL;            
#endif
        } else {
            gz = NULL;
        }
    }
}

XMLHepRepWriter::~XMLHepRepWriter() {
    delete gz;
    delete zip;
}

bool XMLHepRepWriter::addProperty(std::string key, std::string value) {
    properties[key] = value;
    return true;
}

bool XMLHepRepWriter::close() {
    if (zip != NULL) {
        zip->putNextEntry("heprep.properties", true);
        
        map<string, string>::iterator i = properties.begin();
        while (i != properties.end()) {
            *zip << (*i).first << "=" << (*i).second << endl;
            i++;
        }
        zip->closeEntry();
        zip->close();
    }

    if (gz != NULL) {
        gz->close();
    }
    return true;
}

bool XMLHepRepWriter::write(HepRep* heprep, string name) {
    if (zip != NULL) {
        zip->putNextEntry(name, compress);
    }
    
    if (name.rfind(".bheprep") == name.length()-8) {
       xml = new BHepRepWriter(*out);
    } else { 
       xml = new XMLWriter(out, "  ", NAMESPACE);
    }
            
    xml->openDoc();
    xml->setAttribute("version", (string)"2.0");
    xml->setAttribute("xmlns", (string)"http://java.freehep.org/schemas/heprep/2.0");
    xml->setAttribute("xmlns", "xsi", "http://www.w3.org/2001/XMLSchema-instance");
    xml->setAttribute("xsi", "schemaLocation", "http://java.freehep.org/schemas/heprep/2.0 http://java.freehep.org/schemas/heprep/2.0/HepRep.xsd");
    xml->openTag(nameSpace, "heprep");
    write(heprep->getLayerOrder());
    vector<HepRepTypeTree*> typeTreeSet = heprep->getTypeTreeList();
    for (vector<HepRepTypeTree*>::iterator i1=typeTreeSet.begin(); i1 != typeTreeSet.end(); i1++) {
        write(*i1);
    }
    vector<HepRepInstanceTree*> instanceTreeSet = heprep->getInstanceTreeList();
    for (vector<HepRepInstanceTree*>::iterator i2=instanceTreeSet.begin(); i2 != instanceTreeSet.end(); i2++) {
        write(*i2);
    }
    xml->closeTag();
    xml->closeDoc();
//    xml->close();
    delete xml;
   
    if (zip != NULL) {
        zip->closeEntry();
    }

    return true;
}

bool XMLHepRepWriter::write(vector<string> layers) {
    string layerOrder = "";
    bool comma = false;
    for (vector<string>::iterator i=layers.begin(); i != layers.end(); i++) {
        if (comma) {
            layerOrder.append(", ");
        }
        layerOrder.append(*i);
        comma = true;
    }
    xml->setAttribute("order", layerOrder);
    xml->printTag(nameSpace, "layer");
    return true;
}

bool XMLHepRepWriter::write(HepRepTypeTree* typeTree) {
    xml->setAttribute("name", typeTree->getName());
    xml->setAttribute("version", typeTree->getVersion());
    xml->openTag(nameSpace, "typetree");

    vector<HepRepType*> types = typeTree->getTypeList();
    for (vector<HepRepType*>::iterator i=types.begin(); i != types.end(); i++) {
        write(*i);
    }
    
    xml->closeTag();
    return true;
}

bool XMLHepRepWriter::write(HepRepType* type) {
    xml->setAttribute("name", type->getName());
    xml->openTag(nameSpace, "type");
    write((HepRepDefinition*)type);
    write((HepRepAttribute*)type);
    
    vector<HepRepType*> types = type->getTypeList();
    for (vector<HepRepType*>::iterator i=types.begin(); i != types.end(); i++) {
        write(*i);
    }
    xml->closeTag();
    return true;
}

bool XMLHepRepWriter::write(HepRepTreeID* treeID) {
    xml->setAttribute("qualifier", treeID->getQualifier());
    xml->setAttribute("name", treeID->getName());
    xml->setAttribute("version", treeID->getVersion());
    xml->printTag(nameSpace, "treeid");
    return true;
}

bool XMLHepRepWriter::write(HepRepAction* action) {
    xml->setAttribute("name", action->getName());
    xml->setAttribute("expression", action->getExpression());
    xml->printTag(nameSpace, "action");
    return true;
}

bool XMLHepRepWriter::write(HepRepInstanceTree* instanceTree) {
    xml->setAttribute("name", instanceTree->getName());
    xml->setAttribute("version", instanceTree->getVersion());
    xml->setAttribute("typetreename", instanceTree->getTypeTree()->getName());
    xml->setAttribute("typetreeversion", instanceTree->getTypeTree()->getVersion());
    xml->openTag(nameSpace, "instancetree");
    // refs
    vector<HepRepTreeID*> instanceTreeSet = instanceTree->getInstanceTreeList();
    for (vector<HepRepTreeID*>::iterator i1=instanceTreeSet.begin(); i1 != instanceTreeSet.end(); i1++) {
        write(*i1);
    }

    // instances
    vector<HepRepInstance*> instanceList = instanceTree->getInstances();
    for (vector<HepRepInstance*>::iterator i2=instanceList.begin(); i2 != instanceList.end(); i2++) {
        write(*i2);
    }
    xml->closeTag();
    return true;
}

bool XMLHepRepWriter::write(HepRepInstance* instance) {
    // FIXME FREEHEP-356
    xml->setAttribute("type", instance->getType()->getFullName());
    xml->openTag(nameSpace, "instance");
    write((HepRepAttribute*)instance);

    vector<HepRepPoint*> pointList = instance->getPoints();
    for (vector<HepRepPoint*>::iterator i1=pointList.begin(); i1 != pointList.end(); i1++) {
        write(*i1);
    }

    vector<HepRepInstance*> instanceList = instance->getInstances();
    for (vector<HepRepInstance*>::iterator i2=instanceList.begin(); i2 != instanceList.end(); i2++) {
        write(*i2);
    }
    xml->closeTag();
    return true;
}

bool XMLHepRepWriter::write(HepRepPoint* point) {
    xml->setAttribute("x", point->getX());
    xml->setAttribute("y", point->getY());
    xml->setAttribute("z", point->getZ());
    if (point->getAttValuesFromNode().size() != 0) {
        xml->openTag(nameSpace, "point");
        write((HepRepAttribute*)point);
        xml->closeTag();
    } else {
        xml->printTag(nameSpace, "point");
    }
    return true;
}

bool XMLHepRepWriter::write(HepRepAttribute* attribute) {
    // BUG FIX.  Do something special for layers, because these do not end
    // up in the normal iteration.
    HepRepAttValue* layerAtt = attribute->getAttValueFromNode("layer");
    if (layerAtt != NULL) write(layerAtt);

    set<HepRepAttValue*> attSet = attribute->getAttValuesFromNode();
    for (set<HepRepAttValue*>::iterator i=attSet.begin(); i != attSet.end(); i++) {
        write(*i);
    }
    return true;
}

bool XMLHepRepWriter::write(HepRepDefinition* definition) {
    set<HepRepAttDef*> list = definition->getAttDefsFromNode();
    for (set<HepRepAttDef*>::iterator i=list.begin(); i != list.end(); i++) {
        write(*i);
    }
    return true;
}

bool XMLHepRepWriter::write(HepRepAttValue* attValue) {
    string name = attValue->getName();

    xml->setAttribute("name", name);

    switch(attValue->getType()) {
        default:                            xml->setAttribute("value", attValue->getAsString());
                                            break;                           
        case HepRepConstants::TYPE_STRING:  xml->setAttribute("value", attValue->getString());
                                            break;
        case HepRepConstants::TYPE_LONG:    xml->setAttribute("value", attValue->getLong());
                                            break;
        case HepRepConstants::TYPE_INT:     xml->setAttribute("value", attValue->getInteger());
                                            break;
        case HepRepConstants::TYPE_DOUBLE:  xml->setAttribute("value", attValue->getDouble());
                                            break;
        case HepRepConstants::TYPE_BOOLEAN: xml->setAttribute("value", attValue->getBoolean());
                                            break;
        case HepRepConstants::TYPE_COLOR:   xml->setAttribute("value", attValue->getColor());
    }

    if (attValue->showLabel() != HepRepConstants::SHOW_NONE) {
        xml->setAttribute("showlabel", attValue->showLabel());
    }

    xml->printTag(nameSpace, "attvalue");
    return true;
}

bool XMLHepRepWriter::write(HepRepAttDef* attDef) {
    xml->setAttribute("name", attDef->getName());
    xml->setAttribute("desc", attDef->getDescription());
    xml->setAttribute("category", attDef->getCategory());
    xml->setAttribute("extra", attDef->getExtra());
    xml->printTag(nameSpace, "attdef");
    return true;
}

} // cheprep

