
#include "XMLHepRepStreamer.h"

#include "StreamerHepRepInstance.h"
#include "DefaultHepRepAttValue.h"

#define NAMESPACE "heprep"

using namespace std;
using namespace HEPREP;

XMLHepRepStreamer::XMLHepRepStreamer(ostream* out)
    : XMLWriter(out, "  ", NAMESPACE) {

    this->nameSpace = NAMESPACE;
    openDoc();
}

XMLHepRepStreamer::~XMLHepRepStreamer() {
}

bool XMLHepRepStreamer::close() {
    while (!hepreps.empty()) {
        hepreps.pop();
        closeTag();
    }

    closeDoc();
    XMLWriter::close();
    return true;
}

bool XMLHepRepStreamer::write(HepRep* root) {
    setAttribute("xmlns", "http://www.freehep.org/HepRep");
    setAttribute("xmlns", "xsi", "http://www.w3.org/2001/XMLSchema-instance");
    setAttribute("xsi", "schemaLocation", "HepRep.xsd");
    openTag(nameSpace, "heprep");
    writtenLayers = false;
    hepreps.push(root);
    this->heprep = root;
    return true;
}

bool XMLHepRepStreamer::write(HepRepTypeTree* typeTree) {
    while (hepreps.top() != heprep) {
        if (hepreps.empty()) {
            cerr << "XMLHepRepStreamer: cannot write typeTree without HepRep" << endl;
            return false;
        }
        hepreps.pop();
        closeTag();
    }
    if (!writtenLayers) {
        writtenLayers = true;
        string layerOrder = "";
        bool comma = false;
        vector<string> *layers = heprep->getLayerOrder();
        for (vector<string>::iterator i=layers->begin(); i != layers->end(); i++) {
            if (comma) {
                layerOrder.append(", ");
            }
            layerOrder.append(*i);
            comma = true;
        }
        setAttribute("order", layerOrder);
        printTag(nameSpace, "layer");
    }

    setAttribute("name", typeTree->getName());
    setAttribute("version", typeTree->getVersion());
    openTag(nameSpace, "typetree");
    hepreps.push(typeTree);
    heprepTypeTree = typeTree;
    return true;
}

bool XMLHepRepStreamer::write(HepRepType* type) {
    HepRepType* parent = type->getSuperType();
    if (parent != NULL) {
        while (hepreps.top() != parent) {
            if (hepreps.empty()) {
                cerr << "XMLHepRepStreamer: cannot write type, parent not anymore on stack..." << endl;
                return false;
            }
            hepreps.pop();
            closeTag();
        }
    } else {
        while (hepreps.top() != heprepTypeTree) {
            if (hepreps.empty()) {
                cerr << "XMLHepRepStreamer: cannot write type without parent typeTree." << endl;
                return false;
            }
            hepreps.pop();
            closeTag();
        }
    }

    setAttribute("name", type->getName());
    openTag(nameSpace, "type");
    hepreps.push(type);
    return true;
}

bool XMLHepRepStreamer::write(HepRepTreeID* treeID) {
    setAttribute("qualifier", treeID->getQualifier());
    setAttribute("name", treeID->getName());
    setAttribute("version", treeID->getVersion());
    printTag(nameSpace, "treeid");
    return true;
}

bool XMLHepRepStreamer::write(HepRepAction* action) {
    setAttribute("name", action->getName());
    setAttribute("expression", action->getExpression());
    printTag(nameSpace, "action");
    return true;
}

bool XMLHepRepStreamer::write(HepRepInstanceTree* instanceTree) {
    while (hepreps.top() != heprep) {
        if (hepreps.empty()) {
            cerr << "XMLHepRepStreamer cannot write instanceTree without HepRep." << endl;
            return false;
        }
        hepreps.pop();
        closeTag();
    }

    setAttribute("name", instanceTree->getName());
    setAttribute("version", instanceTree->getVersion());
    setAttribute("typetreename", instanceTree->getTypeTree()->getName());
    setAttribute("typetreeversion", instanceTree->getTypeTree()->getVersion());
    openTag(nameSpace, "instancetree");
    hepreps.push(instanceTree);
    return true;
}

bool XMLHepRepStreamer::write(HepRepInstance* instance) {
    void* parent = dynamic_cast<StreamerHepRepInstance*>(instance)->getParent();
    while (hepreps.top() != parent) {
        if (hepreps.empty()) {
            cerr << "XMLHepRepStreamer: cannot write instance, parent (instance/instanceTree) not anymore on Stack..." << endl;
            return false;
        }
        hepreps.pop();
        closeTag();
    }

    setAttribute("type", instance->getType()->getName());
    openTag(nameSpace, "instance");
    hepreps.push(instance);
    return true;
}

bool XMLHepRepStreamer::write(HepRepPoint* point) {
    while (hepreps.top() != point->getInstance()) {
        if (hepreps.empty()) {
            cerr << "XMLHepRepStreamer: cannot write point, parent instance not anymore on Stack..." << endl;
            return false;
        }
        hepreps.pop();
        closeTag();
    }

    setAttribute("x", point->getX());
    setAttribute("y", point->getY());
    setAttribute("z", point->getZ());
    openTag(nameSpace, "point");
    hepreps.push(point);
    return true;
}

bool XMLHepRepStreamer::write(HepRepAttribute* attribute) {
    return true;
}

bool XMLHepRepStreamer::write(HepRepDefinition* definition) {
    return true;
}

bool XMLHepRepStreamer::write(HepRepAttValue* attValue) {
    string name = attValue->getName();

    setAttribute("name", name);
    setAttribute("value", attValue->getAsString());
    setAttribute("type", attValue->getTypeName());

    string label = dynamic_cast<DefaultHepRepAttValue*>(attValue)->toShowLabel();
    setAttribute("showLabel", label);
    printTag(nameSpace, "attvalue");
    return true;
}

bool XMLHepRepStreamer::write(HepRepAttDef* attDef) {
    setAttribute("name", attDef->getName());
    setAttribute("desc", attDef->getDescription());
    setAttribute("category", attDef->getCategory());
    setAttribute("extra", attDef->getExtra());
    printTag(nameSpace, "attdef");
    return true;
}

