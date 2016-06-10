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
// $Id: pyTestEm0.cc 66241 2012-12-13 18:34:42Z gunter $
// ====================================================================
//   pyTestEm0.cc
//
//   python wrapper for user application
//                                         2007 Q
// ====================================================================
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ParticleTable.hh"

#include <boost/python.hpp>
#include <boost/python/list.hpp>

#include <vector>
#include <string>

using namespace boost::python;
// ========================================================================================
//   Wrap Gean4 methods  which are not accessible from python because they return a pointer. 
// ========================================================================================
boost::python::list  getParticleTable()
{
	// Create a list on heap which will be return to python
	boost::python::list *particleList = new boost::python::list();

	// Get particle list fron G4ParticleTable
	G4ParticleTable *g4ParticleList = G4ParticleTable::GetParticleTable();

	// Fill python list from g4ParticleList
	for ( int index = 0 ; index <= g4ParticleList->size() ; index++ ) {
		 particleList->append ( (std::string) g4ParticleList->GetParticleName(index) );
	}

	return *particleList;
}
// ====================================================================
boost::python::list getMaterialTable()
	{
		// Create a list on heap which will be return to python
		boost::python::list *materialTableList = new boost::python::list();

		// Get material list fron G4Material
		G4MaterialTable materialList = *G4Material::GetMaterialTable();
		
		std::vector<G4Material*>::iterator itVectorData;
		for(itVectorData = materialList.begin(); itVectorData != materialList.end(); itVectorData++) {
			materialTableList->append   (   (std::string)(*(itVectorData))->GetName())    ;
		
		}
	return *materialTableList;
}


// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(testem0) {

  def ("getMaterialTable", getMaterialTable);

  def ("getParticleTable", getParticleTable);

  class_<DetectorConstruction, DetectorConstruction*,
    bases<G4VUserDetectorConstruction> >
    ("DetectorConstruction", "testEm0 detector")
	.def("SetMaterial",&DetectorConstruction::SetMaterial)
    ;

  class_<PrimaryGeneratorAction, PrimaryGeneratorAction*,
    bases<G4VUserPrimaryGeneratorAction> >
    ("PrimaryGeneratorAction", init<DetectorConstruction*>())
    ;

  class_<RunAction, RunAction*,
	 bases<G4UserRunAction> >
	 ("RunAction", init<DetectorConstruction*, PrimaryGeneratorAction*>())
    ;

  class_<PhysicsList, PhysicsList*,
    bases<G4VUserPhysicsList> >
    ("PhysicsList", "testEm0 physics list")
    ;
}
