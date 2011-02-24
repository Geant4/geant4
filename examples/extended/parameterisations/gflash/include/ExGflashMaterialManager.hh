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
#ifndef MATERIALMANAGER_H
#define MATERIALMANAGER_H

#include "G4Element.hh"
#include "G4Material.hh"

#include <map>
#include <iostream>
#include <functional>

typedef std::map< G4String, G4Element*, std::less<G4String> > ElementList;
typedef std::map< G4String, G4Material*, std::less<G4String> > MaterialList;




class ExGflashMaterialManager {
private:
	ElementList elist;
	MaterialList mlist;
	ExGflashMaterialManager() {initialize();}
	static ExGflashMaterialManager *mpointer;
public:
	static ExGflashMaterialManager* GetMaterialManager()
	{
		if (!mpointer)
			mpointer=new ExGflashMaterialManager;
		return mpointer;
	}
	void storeElement(G4String, G4String, double, double);
	G4Element* getElement(G4String);
	void storeMaterial(G4String, double, double, double);
	void storeMaterial(G4String, double, int);
	G4Material* getMaterial(G4String);
	void addMaterial(G4String,G4String,double);
	void addElement(G4String,G4String,double);
	void addElement(G4String,G4String,int);
	void printElementTable();
	void printMaterialTable();
	void initialize();
};
#endif
