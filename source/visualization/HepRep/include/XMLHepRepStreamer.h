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

#include <iostream>
#include <string>
#include <stack>

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
#include "HEPREP/HepRepAttDef.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepAttribute.h"

#include "XMLWriter.h"

/**
 *
 * @author M.Donszelmann
 */
class XMLHepRepStreamer : public XMLWriter, public virtual HEPREP::HepRepWriter {

    private:
        std::string nameSpace;
        bool writtenLayers;
        std::stack<void *> hepreps;
        HEPREP::HepRep *heprep;
        HEPREP::HepRepTypeTree *heprepTypeTree;

    public:
        XMLHepRepStreamer(std::ostream* out);
        ~XMLHepRepStreamer();

        bool close();
        bool write(HEPREP::HepRep* heprep);
        bool write(HEPREP::HepRepTypeTree* typeTree);
        bool write(HEPREP::HepRepType* type);
        bool write(HEPREP::HepRepTreeID* treeID);
        bool write(HEPREP::HepRepAction* action);
        bool write(HEPREP::HepRepInstanceTree* instanceTree);
        bool write(HEPREP::HepRepInstance* instance);
        bool write(HEPREP::HepRepPoint* point);
        bool write(HEPREP::HepRepAttribute* attribute);
        bool write(HEPREP::HepRepDefinition* definition);
        bool write(HEPREP::HepRepAttValue* attValue);
        bool write(HEPREP::HepRepAttDef* attDef);
};

#endif
