
#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepAttribute.h"
#include "HEPREP/HepRepConstants.h"
#include "HEPREP/HepRepFactory.h"
#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepInstanceTree.h"
#include "HEPREP/HepRepPoint.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepTypeTree.h"

#include "XMLHepRepStreamerFactory.h"

using namespace std;
using namespace HEPREP;

int main(int argc, char** argv) {

    char* fname = "HepRepStreamerTest.xml";

    printf("Creating Factory...\n");
	HepRepFactory* factory = new XMLHepRepStreamerFactory();

    printf("Creating Writer...\n");
    ostream* out = new ofstream(fname);
	HepRepWriter* writer = factory->createHepRepWriter(out);

    printf("Creating HepRep...\n");
	HepRep* heprep = factory->createHepRep();

    printf("Creating Layer...\n");
	string layer = "Detector";
	heprep->addLayer(layer);

    printf("Creating HepRepTreeID...\n");
    HepRepTreeID* treeID = factory->createHepRepTreeID("CylinderType", "1.0");

//    <heprep:typetree name="CylinderType" version="1.0">
    printf("Creating HepRepTypeTree...\n");
    HepRepTypeTree* typeTree = factory->createHepRepTypeTree(treeID);
    heprep->addTypeTree(typeTree);

    printf("Creating HepRepType...\n");
    HepRepType* type = factory->createHepRepType(NULL, "Cylinder");
    typeTree->addType(type);
//    <heprep:attvalue name="banner" value="true" />
//    <heprep:attvalue name="framed" value="true" />
//    <heprep:attvalue name="drawAs" value="Cylinder" />
//    <heprep:attvalue name="radius1" value="2" />
//    <heprep:attvalue name="radius2" value="2" />
    type->addAttValue("banner", true);
    type->addAttValue("framed", true);
    type->addAttValue("drawAs", string("Cylinder"));
    type->addAttValue("radius1", 2.0);
    type->addAttValue("radius2", 2.0, HepRepConstants::SHOW_VALUE + HepRepConstants::SHOW_NAME + 0x0400);

    printf("Creating HepRepInstanceTree...\n");
//    <heprep:instancetree name="TestCylinder" version="CPlusPlus Generated" typename="CylinderType" typeversion="1.0">
    HepRepInstanceTree* instanceTree = factory->createHepRepInstanceTree("TestCylinder", "CPlusPlus Generated", typeTree);
    heprep->addInstanceTree(instanceTree);

    printf("Creating HepRepInstance...\n");
//    <heprep:instance type="Cylinder">
    HepRepInstance* instance1 = factory->createHepRepInstance(instanceTree, type);

//    <heprep:attvalue name="color" value="red" />
//    <heprep:attvalue name="label" value="x-axis" showLabel="VALUE" />
    vector<double> red;
    red.push_back(1.0);
    red.push_back(0.0);
    red.push_back(0.0);
    instance1->addAttValue("color", red);
    instance1->addAttValue("label", string("x-axis"), HepRepConstants::SHOW_VALUE);

//    <heprep:point x="0.0" y="0.0" z="0.0" />
//    <heprep:point x="4.0" y="0.0" z="0.0" />
    factory->createHepRepPoint(instance1, 0, 0, 0);
    factory->createHepRepPoint(instance1, 4, 0, 0);

//    <heprep:instance type="Cylinder">
    HepRepInstance* instance2 = factory->createHepRepInstance(instanceTree, type);

//    <heprep:attvalue name="color" value="green" />
//    <heprep:attvalue name="label" value="y-axis" showLabel="VALUE" />
    vector<double> green;
    green.push_back(0.0);
    green.push_back(1.0);
    green.push_back(0.0);
    instance2->addAttValue("color", green, 0);
    instance2->addAttValue("label", string("y-axis"), HepRepConstants::SHOW_VALUE);

//    <heprep:point x="0.0" y="0.0" z="0.0" />
//    <heprep:point x="0.0" y="4.0" z="0.0" />
    factory->createHepRepPoint(instance2, 0, 0, 0);
    factory->createHepRepPoint(instance2, 0, 4, 0);

//    <heprep:instance type="Cylinder">
    HepRepInstance* instance3 = factory->createHepRepInstance(instanceTree, type);

//    <heprep:attvalue name="color" value="blue" />
//    <heprep:attvalue name="label" value="z-axis" showLabel="VALUE" />
    vector<double> blue;
    blue.push_back(0.0);
    blue.push_back(1.0);
    blue.push_back(0.0);
    instance3->addAttValue("color", blue, 0);
    instance3->addAttValue("label", string("z-axis"), HepRepConstants::SHOW_VALUE);

//    <heprep:point x="0.0" y="0.0" z="0.0" />
//    <heprep:point x="0.0" y="0.0" z="4.0" />
    factory->createHepRepPoint(instance3, 0, 0, 0);
    factory->createHepRepPoint(instance3, 0, 0, 4);

//    <heprep:instance type="Cylinder">
    HepRepInstance* instance4 = factory->createHepRepInstance(instanceTree, type);

//    <heprep:attvalue name="Phi1" value="0.2" showLabel="false" />
//    <heprep:attvalue name="Phi2" value="0.3" showLabel="false" />
    instance4->addAttValue("Phi1", 0.2, 0);
    instance4->addAttValue("Phi2", 0.3, 0);

//    <heprep:point x="10.0" y="10.0" z="10.0" />
//    <heprep:point x="12.3" y="12.3" z="12.3" />
    factory->createHepRepPoint(instance4, 10, 10, 10);
    factory->createHepRepPoint(instance4, 12.3, 12.3, 12.3);

    printf("Closing Writer...\n");
    writer->close();

    printf("Test finished ok.\n");

    return 0;
}
