//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#include <cmath>

#include "StreamerHepRepPoint.h"

using namespace std;
using namespace HEPREP;

StreamerHepRepPoint::StreamerHepRepPoint(HepRepWriter* stream, HepRepInstance* inst, double x, double y, double z)
    : StreamerHepRepAttribute(stream), instance(inst), x(x), y(y), z(z) {

    if (instance == NULL) {
        cerr << "HepRepPoints cannot be created without a HepRepInstance." << endl;
    }

    stream->write(this);
}

StreamerHepRepPoint::~StreamerHepRepPoint() {
}

HepRepInstance* StreamerHepRepPoint::getInstance() {
    return instance;
}

HepRepAttValue* StreamerHepRepPoint::getAttValue(string lowerCaseName) {
    HepRepAttValue* value = getAttValueFromNode(lowerCaseName);
    return (value != NULL) ? value : instance->getAttValue(lowerCaseName);
}

HepRepPoint* StreamerHepRepPoint::copy(HepRepInstance* parent) {
    return NULL;
}

double StreamerHepRepPoint::getX() {
    return x;
}

double StreamerHepRepPoint::getY() {
    return y;
}

double StreamerHepRepPoint::getZ() {
    return z;
}

vector<double>* StreamerHepRepPoint::getXYZ(vector<double>* xyz) {
    (*xyz)[0] = x;
    (*xyz)[1] = y;
    (*xyz)[2] = z;
    return xyz;
}

double StreamerHepRepPoint::getRho() {
    return sqrt(x*x + y*y);
}

double StreamerHepRepPoint::getPhi() {
    return atan2(y, x);
}

double StreamerHepRepPoint::getTheta() {
    return atan2(getRho(), z);
}

double StreamerHepRepPoint::getR() {
    double r = getRho();
    return sqrt(r*r + z*z);
}

double StreamerHepRepPoint::getEta() {
    double ct = cos(getTheta());
    return -0.5*log((1.-ct)/(1.+ct));
}

double StreamerHepRepPoint::getX(double xVertex, double yVertex, double zVertex) {
    return x - xVertex;
}

double StreamerHepRepPoint::getY(double xVertex, double yVertex, double zVertex) {
    return y - yVertex;
}

double StreamerHepRepPoint::getZ(double xVertex, double yVertex, double zVertex) {
    return z - zVertex;
}

double StreamerHepRepPoint::getRho(double xVertex, double yVertex, double zVertex) {
    double dx = getX(xVertex, yVertex, zVertex);
    double dy = getY(xVertex, yVertex, zVertex);
    return sqrt(dx*dx + dy*dy);
}

double StreamerHepRepPoint::getPhi(double xVertex, double yVertex, double zVertex) {
    return atan2(getY(xVertex, yVertex, zVertex), getX(xVertex, yVertex, zVertex));
}

double StreamerHepRepPoint::getTheta(double xVertex, double yVertex, double zVertex) {
    return atan2(getRho(xVertex, yVertex, zVertex), getZ(xVertex, yVertex, zVertex));
}

double StreamerHepRepPoint::getR(double xVertex, double yVertex, double zVertex) {
    double dr = getRho(xVertex, yVertex, zVertex);
    double dz = getZ(xVertex, yVertex, zVertex);
    return sqrt(dr*dr + dz*dz);
}

double StreamerHepRepPoint::getEta(double xVertex, double yVertex, double zVertex) {
    double ct = cos(getTheta(xVertex, yVertex, zVertex));
    return -0.5*log((1.-ct)/(1.+ct));
}


