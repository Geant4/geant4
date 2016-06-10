// Copyright FreeHEP, 2005.

#include <string>
#include <iostream>
#include <cmath>

#include "cheprep/DefaultHepRepPoint.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepPoint.cc 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

DefaultHepRepPoint::DefaultHepRepPoint(HepRepInstance* inst, double xx, double yy, double zz)
    : DefaultHepRepAttribute(), instance(inst), x(xx), y(yy), z(zz) {

    if (instance == NULL) {
        cerr << "HepRepPoints cannot be created without a HepRepInstance." << endl;
    } else {
        instance->addPoint(this);
    }
}

DefaultHepRepPoint::~DefaultHepRepPoint() {
}

HepRepInstance* DefaultHepRepPoint::getInstance() {
    return instance;
}

HepRepAttValue* DefaultHepRepPoint::getAttValue(string lowerCaseName) {
    HepRepAttValue* value = getAttValueFromNode(lowerCaseName);
    return (value != NULL) ? value : instance->getAttValue(lowerCaseName);
}

HepRepPoint* DefaultHepRepPoint::copy(HepRepInstance* inst) {
    return new DefaultHepRepPoint(inst, x, y, z);
}

double DefaultHepRepPoint::getX() {
    return x;
}

double DefaultHepRepPoint::getY() {
    return y;
}

double DefaultHepRepPoint::getZ() {
    return z;
}

vector<double>* DefaultHepRepPoint::getXYZ(vector<double>* xyz) {
    (*xyz)[0] = x;
    (*xyz)[1] = y;
    (*xyz)[2] = z;
    return xyz;
}

double DefaultHepRepPoint::getRho() {
    return sqrt(x*x + y*y);
}

double DefaultHepRepPoint::getPhi() {
    return atan2(y, x);
}

double DefaultHepRepPoint::getTheta() {
    return atan2(getRho(), z);
}

double DefaultHepRepPoint::getR() {
    double r = getRho();
    return sqrt(r*r + z*z);
}

double DefaultHepRepPoint::getEta() {
    double ct = cos(getTheta());
    return -0.5*log((1.-ct)/(1.+ct));
}

double DefaultHepRepPoint::getX(double xVertex, double, double) {
    return x - xVertex;
}

double DefaultHepRepPoint::getY(double, double yVertex, double) {
    return y - yVertex;
}

double DefaultHepRepPoint::getZ(double, double, double zVertex) {
    return z - zVertex;
}

double DefaultHepRepPoint::getRho(double xVertex, double yVertex, double zVertex) {
    double dx = getX(xVertex, yVertex, zVertex);
    double dy = getY(xVertex, yVertex, zVertex);
    return sqrt(dx*dx + dy*dy);
}

double DefaultHepRepPoint::getPhi(double xVertex, double yVertex, double zVertex) {
    return atan2(getY(xVertex, yVertex, zVertex), getX(xVertex, yVertex, zVertex));
}

double DefaultHepRepPoint::getTheta(double xVertex, double yVertex, double zVertex) {
    return atan2(getRho(xVertex, yVertex, zVertex), getZ(xVertex, yVertex, zVertex));
}

double DefaultHepRepPoint::getR(double xVertex, double yVertex, double zVertex) {
    double dr = getRho(xVertex, yVertex, zVertex);
    double dz = getZ(xVertex, yVertex, zVertex);
    return sqrt(dr*dr + dz*dz);
}

double DefaultHepRepPoint::getEta(double xVertex, double yVertex, double zVertex) {
    double ct = cos(getTheta(xVertex, yVertex, zVertex));
    return -0.5*log((1.-ct)/(1.+ct));
}


} // cheprep

