
#include "DefaultHepRepAction.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepAction::DefaultHepRepAction(string name, string expression)
    : name(name), expression(expression) {
}

DefaultHepRepAction::~DefaultHepRepAction() {
}

string DefaultHepRepAction::getName() {
    return name;
}

string DefaultHepRepAction::getExpression() {
    return expression;
}

HepRepAction* DefaultHepRepAction::copy() {
    return NULL;
}

