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
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepAttDef.h"
#include "HEPREP/HepRepAttValue.h"
#include "HEPREP/HepRepWriter.h"

#include "StreamerHepRepDefinition.h"

/**
 *
 * @author M.Donszelmann
 */
class StreamerHepRepType : public StreamerHepRepDefinition, public virtual HEPREP::HepRepType {

    private:
        HEPREP::HepRepType* parent;
        std::string name;
        std::string description;
        std::string infoURL;

    public:
        StreamerHepRepType(HEPREP::HepRepWriter* streamer, HEPREP::HepRepType* parent, std::string name);
        ~StreamerHepRepType();

        HEPREP::HepRepType* getSuperType();
        HEPREP::HepRepAttDef* getAttDef(std::string name);
        HEPREP::HepRepAttValue* getAttValue(std::string name);
        HEPREP::HepRepType* copy(HEPREP::HepRep* heprep, HEPREP::HepRepType* parent);
        std::string getName();
        std::string getDescription();
        void setDescription(std::string description);
        std::string getInfoURL();
        void setInfoURL(std::string infoURL);
        bool addType(HEPREP::HepRepType* type);
        std::vector<HEPREP::HepRepType*>* getTypes();
};

#endif
