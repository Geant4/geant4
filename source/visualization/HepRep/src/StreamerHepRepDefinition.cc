
#include "StreamerHepRepDefinition.h"
#include "DefaultHepRepAttDef.h"

using namespace std;
using namespace HEPREP;

StreamerHepRepDefinition::StreamerHepRepDefinition(HepRepWriter* streamer)
    : StreamerHepRepAttribute(streamer), streamer(streamer) {
}

StreamerHepRepDefinition::~StreamerHepRepDefinition() {
}

vector<HepRepAttDef *>* StreamerHepRepDefinition::getAttDefsFromNode() {
    return NULL;
}

bool StreamerHepRepDefinition::addAttDef(HepRepAttDef* hepRepAttDef) {
    streamer->write(hepRepAttDef);
    delete hepRepAttDef;
    return true;
}

bool StreamerHepRepDefinition::addAttDef(string name, string desc, string type, string extra) {
    addAttDef(new DefaultHepRepAttDef(name, desc, type, extra));
    return true;
}

HepRepAttDef* StreamerHepRepDefinition::getAttDefFromNode(string lowerCaseName) {
    return NULL;
}

HepRepAttDef* StreamerHepRepDefinition::getAttDef(string name) {
    return NULL;
}

