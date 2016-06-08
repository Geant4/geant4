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
//
#ifndef STREAMERHEPREPINSTANCE_H
#define STREAMERHEPREPINSTANCE_H 1

#include "FreeHepTypes.h"

#include <string>
#include <vector>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepSelectFilter.h"
#include "HEPREP/HepRepInstanceTree.h"
#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepPoint.h"
#include "HEPREP/HepRepAttValue.h"

#include "StreamerHepRepAttribute.h"

/**
 *
 * @author M.Donszelmann
 */
class StreamerHepRepInstance : public StreamerHepRepAttribute, public virtual HEPREP::HepRepInstance {

    private:
        void* parent;
        HEPREP::HepRepType* type;

    public:
        StreamerHepRepInstance(HEPREP::HepRepWriter* streamer, HEPREP::HepRepInstance* parent, HEPREP::HepRepType* type);
        StreamerHepRepInstance(HEPREP::HepRepWriter* streamer, HEPREP::HepRepInstanceTree* parent, HEPREP::HepRepType* type);
        ~StreamerHepRepInstance();

        HEPREP::HepRepInstance* copy(HEPREP::HepRep* heprep, HEPREP::HepRepInstance* parent, HEPREP::HepRepSelectFilter* filter);
        HEPREP::HepRepInstance* copy(HEPREP::HepRep* heprep, HEPREP::HepRepInstanceTree* parent, HEPREP::HepRepSelectFilter* filter);
        HEPREP::HepRepType* getType();
        bool addPoint(HEPREP::HepRepPoint* point);
        std::vector<HEPREP::HepRepPoint *>* getPoints();
        bool addInstance(HEPREP::HepRepInstance* instance);
        void removeInstance(HEPREP::HepRepInstance* instance);
        std::vector<HEPREP::HepRepInstance *>* getInstances();
        HEPREP::HepRepAttValue* getAttValue(std::string name);

        void *getParent() { return parent; }
};

#endif
