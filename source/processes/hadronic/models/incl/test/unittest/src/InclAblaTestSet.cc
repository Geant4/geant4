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
// $Id: InclAblaTestSet.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>
#include "InclAblaTestSet.hh"

#include <iostream>

ClassImp(InclAblaTestSet)

//////////////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML
<h1>Container for tests</h1>

<h2>Overview</h2>

<p>
InclAblaTestSet provides a container for tests. All tests must be placed into a
test set. Test sets should contain tests belonging somehow logically together
(e.g. tests for init_incl42 module should be in the same test set). The test
sets provide metadata (currently test set name) to the log writers. This allows
the log writer to e.g. group certain tests together and generate a more
readable log.
</p>

<h2>Examples of usage</h2>
<pre>
// Create a new test set
InclAblaTestSet *aTestSet = new InclAblaTestSet("TestSetName");

// Create new test and place it into a test set
TestRibm *theRibmTest = new TestRibm();
aTestSet->addTest(theRibmTest);

// Add the test set to the test manager
theTestManager->addTestSet(aTestSet);
</pre>
<p>
</p>


  END_HTML*/

InclAblaTestSet::InclAblaTestSet()
{

}

InclAblaTestSet::InclAblaTestSet(TString name)
{
	theTestSetName = name;
}
InclAblaTestSet::~InclAblaTestSet()
{

}

void InclAblaTestSet::addTest(InclAblaTest *aTest)
{
	// Adds a new test to the test set.

	theListOfTests.push_back(aTest);
}

std::list<InclAblaTest*> InclAblaTestSet::getListOfTests()
{
	// Get the list of tests contained by this test set.
	return theListOfTests;
}

TString InclAblaTestSet::getName()
{
	// Returns the name of the software unit tested by this class.
	return theTestSetName;
}

void InclAblaTestSet::setName(TString name)
{
	// Set the name of the code unit to be tested.
 	theTestSetName = name;
}

void InclAblaTestSet::runTests()
{
	// Runs the tests in this set.
	
	std::cout <<"Running tests..." << std::endl;

	std::list<InclAblaTest*>::iterator aTest;

	for(aTest = theListOfTests.begin(); aTest != theListOfTests.end(); aTest++) {
		(*aTest)->runMe();
	}
}
