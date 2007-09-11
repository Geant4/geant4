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
// $Id: InclAblaTestConsoleLogWriter.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>
#include "InclAblaTestConsoleLogWriter.hh"

#include <iostream>

ClassImp(InclAblaTestConsoleLogWriter)

//////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test log output in the console window</h1>

     <h2>Overview</h2>

     <p>
     This logger module writes log to the console window. It also serves as a simple
     example of a log writer. 
     </p>

     END_HTML*/

  InclAblaTestConsoleLogWriter::InclAblaTestConsoleLogWriter()
{

}

InclAblaTestConsoleLogWriter::~InclAblaTestConsoleLogWriter()
{

}

void InclAblaTestConsoleLogWriter::beginLog()
{
  // Write the log header.

  std::cout <<"INCL/ABLA unit test log" << std::endl;
  std::cout <<"--------------------------------------------------" << std::endl;
  std::cout << std::endl;
}

void InclAblaTestConsoleLogWriter::endLog()
{
  // Print end of log message.

  std::cout << std::endl;
  std::cout <<"End of test log." << std::endl;
}

void InclAblaTestConsoleLogWriter::beginTestSet(TString name)
{
  // Print some header information at the beginning of each test set output.
	
  std::cout <<"Test set: " << name << std::endl;
}

void InclAblaTestConsoleLogWriter::endTestSet()
{
  // Information at the end of a test set.

  std::cout <<"End of test set" << std::endl;
}

void InclAblaTestConsoleLogWriter::logTest(InclAblaTest *aTest)
{
  TString *testStatus;
  
  // Writes the information of a test to the log.
  if(aTest->getTestStatus()) {
    testStatus = new TString("pass");
  }
  else {
    testStatus = new TString("fail");
  }

  std::cout << aTest->getUnitName() << " "<< (*testStatus) << std::endl;

  delete testStatus;
}
