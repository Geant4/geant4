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

#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepTypeTree.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepWriter.h"
#include "HEPREP/HepRepTreeID.h"

#include "DefaultHepRepTreeID.h"

/**
 *
 * @author M.Donszelmann
 */
class StreamerHepRepTypeTree : public DefaultHepRepTreeID, public virtual HEPREP::HepRepTypeTree {

    public:
        StreamerHepRepTypeTree(HEPREP::HepRepWriter* streamer, HEPREP::HepRepTreeID* typeTree);
        ~StreamerHepRepTypeTree();

        HEPREP::HepRepTypeTree* copy(HEPREP::HepRep* heprep);
        HEPREP::HepRepTreeID* copy();
        bool addType(HEPREP::HepRepType* type);
        std::vector<HEPREP::HepRepType* >* getTypes();
};

#endif
