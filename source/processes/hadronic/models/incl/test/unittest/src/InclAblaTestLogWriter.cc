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
// $Id: InclAblaTestLogWriter.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>
#include "InclAblaTestLogWriter.hh"

ClassImp(InclAblaTestLogWriter)

//////////////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML
<h1>Base class for INCL4/ABLA unit test log writer classes</h1>

<h2>Overview</h2>

<p>
This class implements the common properties of the INCL4/ABLA test log writer
classes. It also defines an interface that can be used by an instance
of the <a href="InclAblaTestManager.html">InclAblaTestManager</a> class to
write the test log. This abstraction of the interface allows us to cleanly
separate the testing and log writing facilities. One important advantage of
this technique is that it allows us to write log in several formats (e.g. HTML,
plain text, console output, etc).
</p>

<h2>Connecting log writer to the testing framework</h2>

<p>
Log writers must be registered with an instance of the <a
href="InclAblaTestManager.html">InclAblaTestManager</a> class using method <a
href="InclAblaTestManager.html#InclAblaTestManager:registerLogWriter">InclAblaTestManager::registerLogWriter()</a>.
</p>

<h2>Examples of usage</h2>

<p>
In order to use a log writer module the writer must first be instantiated and
later registered to the <a href="InclAblaTestManager.html">test manager</a>.
This is demonstrated by the code below:
</p>
<p>
<pre>
// Here we write the log to the console
InclAblaTestConsoleLogWriter *theLogWriter = new InclAblaConsoleLogWriter();
theTestManager-><a
href="InclAblaTestManager.html#InclAblaTestManager:registerLogWriter">registerLogWriter</a>(theLogWriter);

theTestManager-><a href="InclAblaTestManager.html#InclAblaTestManager:writeLog">writeLog</a>();
</pre>
</p>
  END_HTML*/

InclAblaTestLogWriter::InclAblaTestLogWriter()
{

}

InclAblaTestLogWriter::~InclAblaTestLogWriter()
{

}

void InclAblaTestLogWriter::createLog(std::list<InclAblaTestSet*> aListOfTestSets)
{
	// Generates and writes the log. The log is generated by going through the
	// list of test sets and tests. The list processing is done by this base
	// class and the actual log writing is done using abstract methods so this
	// class does not care what the actual output format is! It only provides
	// the basic infrastructure needed for generating the log!
	
	beginLog();

	std::list<InclAblaTestSet*>::iterator aTestSet;
	for(aTestSet = aListOfTestSets.begin(); aTestSet != aListOfTestSets.end(); aTestSet++) {
		processTestSet(*aTestSet);
	}

	endLog();
}

void InclAblaTestLogWriter::processTestSet(InclAblaTestSet *aTestSet) 
{
	// Process one test set. 
	
	beginTestSet(aTestSet->getName());

	std::list<InclAblaTest*> aListOfTests = aTestSet->getListOfTests();
	std::list<InclAblaTest*>::iterator aTest;

	for(aTest = aListOfTests.begin(); aTest != aListOfTests.end(); aTest++) {
		logTest(*aTest);
	}

	endTestSet();
}

