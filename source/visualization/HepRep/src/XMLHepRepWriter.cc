
#include "XMLHepRepWriter.h"

#include "DefaultHepRepInstance.h"
#include "DefaultHepRepAttValue.h"

#define NAMESPACE "heprep"

using namespace zipios;
using namespace std;
using namespace HEPREP;

XMLHepRepWriter::XMLHepRepWriter(ostream* out, bool randomAccess, bool compress) {

    this->nameSpace = NAMESPACE;

    if (randomAccess) {
        zip = new ZipOutputStream(*out);
        zip->setLevel(compress ? 6 : 0);
        out = zip;
        gz = NULL;
    } else {
        zip = NULL;
        if (compress) {
            gz = new GZIPOutputStream(*out);
            out = gz;
        } else {
            gz = NULL;
        }
    }
    xml = new XMLWriter(out, "  ", NAMESPACE);
}

XMLHepRepWriter::~XMLHepRepWriter() {
    delete gz;
    delete zip;
    delete xml;
}

bool XMLHepRepWriter::addProperty(std::string key, std::string value) {
    properties[key] = value;
    return true;
}

bool XMLHepRepWriter::close() {
    xml->closeDoc(true);
    if (zip != NULL) {
        zip->putNextEntry(ZipCDirEntry("heprep.properties"));
        
        map<string, string>::iterator i = properties.begin();
        while (i != properties.end()) {
            *zip << (*i).first << "=" << (*i).second << endl;
            i++;
        }
        zip->closeEntry();
        zip->finish();
        zip->close();
    }
    if (gz != NULL) {
        gz->finish();
        gz->close();
    }
    xml->close();
    return true;
}

bool XMLHepRepWriter::write(HepRep* heprep, string name) {
    if (zip != NULL) {
        zip->putNextEntry(ZipCDirEntry(name));
    }
    xml->openDoc();
    xml->setAttribute("xmlns", "http://java.freehep.org/schemas/heprep/2.0");
    xml->setAttribute("xmlns", "xsi", "http://www.w3.org/2001/XMLSchema-instance");
    xml->setAttribute("xsi", "schemaLocation", "http://java.freehep.org/schemas/heprep/2.0 http://java.freehep.org/schemas/heprep/2.0/HepRep.xsd");
    xml->openTag(nameSpace, "heprep");
    write(heprep->getLayerOrder());
    set<HepRepTypeTree*> typeTreeSet = heprep->getTypeTrees();
    for (set<HepRepTypeTree*>::iterator i1=typeTreeSet.begin(); i1 != typeTreeSet.end(); i1++) {
        write(*i1);
    }
    set<HepRepInstanceTree*> instanceTreeSet = heprep->getInstanceTrees();
    for (set<HepRepInstanceTree*>::iterator i2=instanceTreeSet.begin(); i2 != instanceTreeSet.end(); i2++) {
        write(*i2);
    }
    xml->closeTag();
    xml->closeDoc();
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

    set<HepRepType*> types = typeTree->getTypes();
    for (set<HepRepType*>::iterator i=types.begin(); i != types.end(); i++) {
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
    
    set<HepRepType*> types = type->getTypes();
    for (set<HepRepType*>::iterator i=types.begin(); i != types.end(); i++) {
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
    set<HepRepTreeID*> instanceTreeSet = instanceTree->getInstanceTrees();
    for (set<HepRepTreeID*>::iterator i1=instanceTreeSet.begin(); i1 != instanceTreeSet.end(); i1++) {
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
    xml->setAttribute("value", attValue->getAsString());

    if (attValue->getType() != HepRepConstants::TYPE_STRING) {
        xml->setAttribute("type", attValue->getTypeName());
    }
        
    if (attValue->showLabel() != HepRepConstants::SHOW_NONE) {
        string label = dynamic_cast<DefaultHepRepAttValue*>(attValue)->toShowLabel();
        xml->setAttribute("showlabel", label);
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

