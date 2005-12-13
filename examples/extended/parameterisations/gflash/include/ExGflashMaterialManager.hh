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
#ifndef MATERIALMANAGER_H
#define MATERIALMANAGER_H

using namespace std;
#include "G4Element.hh"
#include "G4Material.hh"

#include <map>
#include <iostream>


typedef map< G4String, G4Element*, less<G4String> > ElementList;
typedef map< G4String, G4Material*, less<G4String> > MaterialList;




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
