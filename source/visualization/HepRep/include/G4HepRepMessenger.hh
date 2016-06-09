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
// $Id: G4HepRepMessenger.hh,v 1.5 2004/05/27 05:55:19 duns Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//
#ifndef G4HepRepMessenger_HH
#define G4HepRepMessenger_HH 1

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"

class G4HepRepMessenger : public G4UImessenger {
    
    public:
        G4HepRepMessenger();
        virtual ~G4HepRepMessenger();

        virtual G4String GetCurrentValue(G4UIcommand * command);
        virtual void SetNewValue(G4UIcommand * command, G4String newValue);
        
        virtual G4String getEventNumberSuffix();
        virtual G4bool appendGeometry();
        virtual G4bool addPointAttributes();
        virtual G4String getCoordinateSystem();


    private:            
        G4UIdirectory* heprepDirectory;
        
        G4String suffix;
        G4UIcmdWithAString* setEventNumberSuffixCommand;
        
        G4bool geometry;
        G4UIcmdWithABool* appendGeometryCommand;

        G4bool pointAttributes;
        G4UIcmdWithABool* addPointAttributesCommand;

        G4String coordinateSystem;
        G4UIcmdWithAString* setCoordinateSystemCommand;
};

#endif
