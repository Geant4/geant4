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
// $Id: runtest.C,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

{
// Main script for running the tests.
	
InclAblaTestManager *theTestManager = new InclAblaTestManager();

InclAblaTestSet *theInclTestSet = new InclAblaTestSet("Cascade: G4Incl");
InclAblaTestSet *theAblaTestSet = new InclAblaTestSet("Evaporation/Fission: G4Abla");

// Add tests
theInclTestSet->addTest(new TestInitIncl());
theInclTestSet->addTest(new TestInitMat());
theInclTestSet->addTest(new TestDensDeut());
theInclTestSet->addTest(new TestInteg());
theInclTestSet->addTest(new TestRibm());
theInclTestSet->addTest(new TestRgaus());
theInclTestSet->addTest(new TestWsax());
theInclTestSet->addTest(new TestDerivWsax());
theInclTestSet->addTest(new TestDmho());
theInclTestSet->addTest(new TestDerivMho());
theInclTestSet->addTest(new TestDerivGaus());
theInclTestSet->addTest(new TestDeutv());
theInclTestSet->addTest(new TestFm2());
theInclTestSet->addTest(new TestDens());
theInclTestSet->addTest(new TestFlin());
theInclTestSet->addTest(new TestFlin2());
theInclTestSet->addTest(new TestTexp());
theInclTestSet->addTest(new TestSigReac());
theInclTestSet->addTest(new TestRadius());
theInclTestSet->addTest(new TestSpl2Ab());
theInclTestSet->addTest(new TestSplineab());
theInclTestSet->addTest(new TestClmb1());
theInclTestSet->addTest(new TestClmb2());
theInclTestSet->addTest(new TestCoulombTransm());
theInclTestSet->addTest(new TestXabs2());
theInclTestSet->addTest(new TestForceAbs());

theTestManager->addTestSet(theInclTestSet);
theTestManager->addTestSet(theAblaTestSet);


InclAblaTestHtmlLogWriter *theHtmlLogWriter = new InclAblaTestHtmlLogWriter();

// Register new log writer:
// Just console output:
//theTestManager->registerLogWriter(new InclAblaTestConsoleLogWriter());
// Console and HTML log output:
theTestManager->registerLogWriter(new InclAblaTestHtmlLogWriter());

theTestManager->runTests();
theTestManager->writeLog();
}
