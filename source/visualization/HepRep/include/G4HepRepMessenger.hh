//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4HepRepMessenger.hh,v 1.9 2010-06-05 06:25:03 perl Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
		static G4HepRepMessenger* GetInstance();  // Singleton constructor.
        virtual ~G4HepRepMessenger();

        virtual G4String GetCurrentValue(G4UIcommand * command);
        virtual void SetNewValue(G4UIcommand * command, G4String newValue);
        
	// Used by HepRepFile
        virtual G4String getFileDir();
        virtual G4String getFileName();
        virtual G4bool getOverwrite();
        virtual G4bool getCullInvisibles();
		virtual G4bool renderCylAsPolygons();
	
	// Used by HepRepXML
        virtual G4String getEventNumberSuffix();
        virtual G4bool appendGeometry();
        virtual G4bool addPointAttributes();
        virtual G4bool useSolids();
        virtual G4bool writeInvisibles();

    private:            
		G4HepRepMessenger();  // Private constructor.
		static G4HepRepMessenger* fpInstance;
        
        G4UIdirectory* heprepDirectory;
        
		G4String fileDir;
		G4UIcmdWithAString* setFileDirCommand;
        
		G4String fileName;
		G4UIcmdWithAString* setFileNameCommand;
        
		G4bool overwrite;
		G4UIcmdWithABool* setOverwriteCommand;
        
		G4bool cullInvisibles;
		G4UIcmdWithABool* setCullInvisiblesCommand;
	
		G4bool cylAsPolygons;
		G4UIcmdWithABool* renderCylAsPolygonsCommand;
	
        G4String suffix;
        G4UIcmdWithAString* setEventNumberSuffixCommand;
        
        G4bool geometry;
        G4UIcmdWithABool* appendGeometryCommand;

        G4bool pointAttributes;
        G4UIcmdWithABool* addPointAttributesCommand;

        G4bool solids;
        G4UIcmdWithABool* useSolidsCommand;
        
        G4bool invisibles;
        G4UIcmdWithABool* writeInvisiblesCommand;
};

#endif
