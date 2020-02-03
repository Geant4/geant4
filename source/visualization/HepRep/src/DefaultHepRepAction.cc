// Copyright FreeHEP, 2005.

#include <iostream>

#include "cheprep/DefaultHepRepAction.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 */
namespace cheprep {

DefaultHepRepAction::DefaultHepRepAction(string aName, string anExpression)
    : name(aName), expression(anExpression) {
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
    return new DefaultHepRepAction(name, expression);
}

} // cheprep

