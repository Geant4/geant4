#ifndef STREAMERHEPREPPOINT_H
#define STREAMERHEPREPPOINT_H 1

#include "FreeHepTypes.h"

#include <string>
#include <vector>

#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepPoint.h"

#include "StreamerHepRepAttribute.h"

/**
 *
 * @author M.Donszelmann
 * @version $Id: StreamerHepRepPoint.h,v 1.1 2002-11-13 07:06:32 duns Exp $
 */
class StreamerHepRepPoint : public StreamerHepRepAttribute, public virtual HEPREP::HepRepPoint {

    private:
        HEPREP::HepRepInstance* instance;

    protected:
        double x, y, z;

    public:
        StreamerHepRepPoint(HEPREP::HepRepWriter* streamer, HEPREP::HepRepInstance* instance, double x, double y, double z);
        ~StreamerHepRepPoint();

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

#endif
