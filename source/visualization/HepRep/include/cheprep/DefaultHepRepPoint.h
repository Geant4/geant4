// Copyright FreeHEP, 2005.
#ifndef CHEPREP_DEFAULTHEPREPPOINT_H
#define CHEPREP_DEFAULTHEPREPPOINT_H 1

#include "cheprep/config.h"

#include <string>
#include <vector>

#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepPoint.h"

#include "DefaultHepRepAttribute.h"

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepPoint.h 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

class DefaultHepRepPoint : public DefaultHepRepAttribute, public virtual HEPREP::HepRepPoint {

    private:
        HEPREP::HepRepInstance* instance;

    protected:
        double x, y, z;

    public:
        DefaultHepRepPoint(HEPREP::HepRepInstance* instance, double x, double y, double z);
        ~DefaultHepRepPoint();

        HEPREP::HepRepInstance* getInstance();

        HEPREP::HepRepAttValue* getAttValue(std::string lowerCaseName);

        HEPREP::HepRepPoint* copy(HEPREP::HepRepInstance* parent);
        double getX();
        double getY();
        double getZ();
        std::vector<double>* getXYZ(std::vector<double>* xyz);
        double getRho();
        double getPhi();
        double getTheta();
        double getR();
        double getEta();
        double getX(double xVertex, double yVertex, double zVertex);
        double getY(double xVertex, double yVertex, double zVertex);
        double getZ(double xVertex, double yVertex, double zVertex);
        double getRho(double xVertex, double yVertex, double zVertex);
        double getPhi(double xVertex, double yVertex, double zVertex);
        double getTheta(double xVertex, double yVertex, double zVertex);
        double getR(double xVertex, double yVertex, double zVertex);
        double getEta(double xVertex, double yVertex, double zVertex);
};

} // cheprep


#endif
