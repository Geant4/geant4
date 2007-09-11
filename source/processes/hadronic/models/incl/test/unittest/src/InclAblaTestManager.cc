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
// $Id: InclAblaTestManager.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>
#include "InclAblaTestManager.hh"

ClassImp(InclAblaTestManager)

//////////////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML
<h1>Manager class for controlling all tests</h1>

<h2>Overview</h2>

<p>
This class manages all tests and log writers. The Tests themselves are added to
test sets. A test set is simply a group of tests that somehow belong together
(e.g. tests for all incl42.cc routines form one test set). All test sets and
one log writer must be registered with the test manager. This arrangement
allows us to have one central point where we can control testing and log
writing tasks.
</p>

<h2>Examples of usage</h2>

<p>
<pre>
// Create a test manager
InclAblaTestManager *theTestManager = new InclAblaTestManager();

// Register new log writer
theTestManager->registerLogWriter(new InclAblaTestConsoleLogWriter());

// Create a test set
InclAblaTestSet *inclTestSet = new InclAblaTestSet("incl42.cc");

// Add tests to the set
inclTestSet->addTest(new TestRibm());

// Add test set to the test manager
theTestManager->addTestSet(inclTestSet);

// Run tests
theTestManager->runTests();

// Write log
theTestManager->writeLog();

</pre>
</p>


  END_HTML*/

InclAblaTestManager::InclAblaTestManager()
{

}

InclAblaTestManager::~InclAblaTestManager()
{

}

void InclAblaTestManager::addTestSet(InclAblaTestSet *aTestSet)
{
	// Add a test set.

	theListOfTestSets.push_back(aTestSet);
}

void InclAblaTestManager::registerLogWriter(InclAblaTestLogWriter *aLogWriter)
{
	// Register the log writer you wish to use in order to get suitable log
	// output. Only one log writer can be in use. If you wish to generate the
	// log in multiple formats you have to first register one log writer, then
	// call writeLog(), register another writer, call writeLog() again, etc.

	theLogWriter = aLogWriter;
}

std::list<InclAblaTestSet*> InclAblaTestManager::getListOfTestSets()
{
	// Returns the list of test sets registered with the testing manager.

	return theListOfTestSets;
}

void InclAblaTestManager::runTests()
{
	// Runs the tests in registered test sets.
	
	std::list<InclAblaTestSet*>::iterator aTestSet;
	for(aTestSet = theListOfTestSets.begin(); aTestSet != theListOfTestSets.end(); aTestSet++) {
		(*aTestSet)->runTests();
	}
}

void InclAblaTestManager::writeLog()
{
	// Calls the registered InclAblaTestLogWriter that writes the log file.
	
	theLogWriter->createLog(theListOfTestSets);
}
