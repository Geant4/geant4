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
#ifndef XMLHEPREPSTREAMERFACTORY_H
#define XMLHEPREPSTREAMERFACTORY_H 1

#include "FreeHepTypes.h"

#include <string>
#include <iostream>

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepFactory.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepPoint.h"
#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepInstanceTree.h"
#include "HEPREP/HepRepTreeID.h"
#include "HEPREP/HepRepAction.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepTypeTree.h"

/**
 *
 * @author M.Donszelmann
 */
class XMLHepRepStreamerFactory : public virtual HEPREP::HepRepFactory {

    private:
        HEPREP::HepRepWriter* streamer;

    public:
        XMLHepRepStreamerFactory();
        ~XMLHepRepStreamerFactory();

//        static HEPREP::HepRepFactory* create();
        HEPREP::HepRepWriter* createHepRepWriter (std::ostream* out);
        HEPREP::HepRepPoint* createHepRepPoint (HEPREP::HepRepInstance* instance,
                                   double x, double y, double z);
        HEPREP::HepRepInstance* createHepRepInstance (HEPREP::HepRepInstance* parent, HEPREP::HepRepType* type);
        HEPREP::HepRepInstance* createHepRepInstance (HEPREP::HepRepInstanceTree* parent, HEPREP::HepRepType* type);
        HEPREP::HepRepTreeID* createHepRepTreeID (std::string name, std::string version, std::string qualifier = "top-level");
        HEPREP::HepRepAction* createHepRepAction (std::string name, std::string expression);
        HEPREP::HepRepInstanceTree* createHepRepInstanceTree (std::string name, std::string version,
                                                        HEPREP::HepRepTreeID* typeTreeID);
        HEPREP::HepRepType* createHepRepType (HEPREP::HepRepType* parent, std::string name);
        HEPREP::HepRepTypeTree* createHepRepTypeTree (HEPREP::HepRepTreeID* treeID);
        HEPREP::HepRep* createHepRep ();
};

#endif
