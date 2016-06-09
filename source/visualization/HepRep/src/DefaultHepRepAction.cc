// Copyright FreeHEP, 2005.

#include <iostream>

#include "cheprep/DefaultHepRepAction.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepAction.cc,v 1.8 2005-06-02 21:28:45 duns Exp $
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

