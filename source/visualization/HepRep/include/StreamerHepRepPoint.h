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
