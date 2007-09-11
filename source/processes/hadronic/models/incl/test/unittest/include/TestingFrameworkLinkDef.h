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
// $Id: TestingFrameworkLinkDef.h,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Manager classes
#pragma link C++ class InclAblaTest+;
#pragma link C++ class InclAblaTestManager+;
#pragma link C++ class InclAblaTestSet+;

// Logging functionality
#pragma link C++ class InclAblaTestLogWriter+;
#pragma link C++ class InclAblaTestConsoleLogWriter+;
#pragma link C++ class InclAblaTestHtmlLogWriter+;

// Util classes
#pragma link C++ class InclAblaTestDifferencePlotter+;

// Testing classes

// INCL4
#pragma link C++ class TestWsax+;
#pragma link C++ class TestDerivWsax+;
#pragma link C++ class TestDmho+;
#pragma link C++ class TestDerivMho+;
#pragma link C++ class TestDerivGaus+;
#pragma link C++ class TestFm2+;
#pragma link C++ class TestDeutv+;
#pragma link C++ class TestDens+;
#pragma link C++ class TestSpl2Ab+;
#pragma link C++ class TestSplineab+;
#pragma link C++ class TestFlin+;
#pragma link C++ class TestFlin2+;
#pragma link C++ class TestInteg+;
#pragma link C++ class TestDensDeut+;
#pragma link C++ class TestInitMat+;
#pragma link C++ class TestInitIncl+;

#pragma link C++ class TestRibm+;
#pragma link C++ class TestRgaus+;
#pragma link C++ class TestTexp+;
#pragma link C++ class TestSigReac+;
#pragma link C++ class TestRadius+;
#pragma link C++ class TestForceAbs+;
#pragma link C++ class TestXabs2+;
#pragma link C++ class TestCoulombTransm+;
#pragma link C++ class TestClmb1+;
#pragma link C++ class TestClmb2+;


#endif
