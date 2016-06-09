
#include <string>
#include <iostream>
#include <cmath>

#include "DefaultHepRepPoint.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepPoint::DefaultHepRepPoint(HepRepInstance* inst, double x, double y, double z)
    : DefaultHepRepAttribute(), instance(inst), x(x), y(y), z(z) {

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
    return std::sqrt(x*x + y*y);
}

double DefaultHepRepPoint::getPhi() {
    return std::atan2(y, x);
}

double DefaultHepRepPoint::getTheta() {
    return std::atan2(getRho(), z);
}

double DefaultHepRepPoint::getR() {
    double r = getRho();
    return std::sqrt(r*r + z*z);
}

double DefaultHepRepPoint::getEta() {
    double ct = std::cos(getTheta());
    return -0.5*std::log((1.-ct)/(1.+ct));
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
    return std::sqrt(dx*dx + dy*dy);
}

double DefaultHepRepPoint::getPhi(double xVertex, double yVertex, double zVertex) {
    return std::atan2(getY(xVertex, yVertex, zVertex), getX(xVertex, yVertex, zVertex));
}

double DefaultHepRepPoint::getTheta(double xVertex, double yVertex, double zVertex) {
    return std::atan2(getRho(xVertex, yVertex, zVertex), getZ(xVertex, yVertex, zVertex));
}

double DefaultHepRepPoint::getR(double xVertex, double yVertex, double zVertex) {
    double dr = getRho(xVertex, yVertex, zVertex);
    double dz = getZ(xVertex, yVertex, zVertex);
    return std::sqrt(dr*dr + dz*dz);
}

double DefaultHepRepPoint::getEta(double xVertex, double yVertex, double zVertex) {
    double ct = std::cos(getTheta(xVertex, yVertex, zVertex));
    return -0.5*std::log((1.-ct)/(1.+ct));
}


