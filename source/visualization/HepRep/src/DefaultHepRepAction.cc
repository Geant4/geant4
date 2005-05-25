// Copyright FreeHEP, 2005.

#include <iostream>

#include "cheprep/DefaultHepRepAction.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepAction.cc,v 1.7 2005-05-25 23:22:25 duns Exp $
 */
namespace cheprep {

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
    return new DefaultHepRepAction(name, expression);
}

} // cheprep

