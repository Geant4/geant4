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
// $Id: InclAblaTestHtmlLogWriter.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>
#include "InclAblaTestHtmlLogWriter.hh"

#include "TDatime.h"

#include <iostream>

ClassImp(InclAblaTestHtmlLogWriter)

//////////////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML
   <h1>Log writer for creating test log in HTML format (NOT FULLY IMPLEMENTED YET!!!)</h1>

   <h2>Overview</h2>

   <p>
   This log writer generates the test log in HTML format.
   </p>

   <h2>Examples of usage</h2>

   <p>
   The HTML log writer implements the interface defined in its superclass <a
   href="InclAblaTestLogWriter.html#overview">InclAblaTestLogWriter</a>. The <a
   href="InclAblaTestLogWriter.html">InclAblaTestLogWriter</a> documentation
   provides an example on how to use log writers.
   </p>
   END_HTML*/

  InclAblaTestHtmlLogWriter::InclAblaTestHtmlLogWriter()
{

}

InclAblaTestHtmlLogWriter::~InclAblaTestHtmlLogWriter()
{

}

void InclAblaTestHtmlLogWriter::beginLog()
{
  // Write the log header.

  outFile.open("htmldoc/testlog.html");

  // Test run version information:
  TDatime *dt = new TDatime();
  Int_t year = dt->GetYear();
  Int_t month = dt->GetMonth();
  Int_t day = dt->GetDay();
  Int_t hour = dt->GetHour();
  Int_t minute = dt->GetMinute();
  char versionstring[400];
  char monthstr[10];
  char daystr[10];
  char minutestr[10];

  if(minute < 10) {
    sprintf(minutestr, "0%i", minute);
  }
  else {
    sprintf(minutestr, "%i", minute);
  }

  if(month < 10) {
    sprintf(monthstr, "0%i", month);
  }
  else {
    sprintf(monthstr, "%i", month);
  }
        
  if(day < 10) {
    sprintf(daystr, "0%i", day);
  }
  else {
    sprintf(daystr, "%i", day);
  }

  sprintf(versionstring, "Test log generated (YYYY-MM-DD hh:mm): %i-%s-%s %i:%s", year, monthstr, daystr, hour, minutestr);

  outFile << "<html>" << std::endl;
  outFile << "<title>INCL/ABLA Unit Test Log</title>" << std::endl;
  outFile << "<body>" << std::endl;
  outFile << "<h1>INCL/ABLA Unit Test Log</h1>" << std::endl;
  outFile << "<br><br>" << std::endl;
  outFile << versionstring << std::endl;
  outFile << "<br><br>" << std::endl;

  std::cout <<"INCL/ABLA unit test log" << std::endl;
  std::cout <<"--------------------------------------------------" << std::endl;
  std::cout << std::endl;
}

void InclAblaTestHtmlLogWriter::endLog()
{
  // Print end of log message.

  outFile << "</body>" << std::endl;
  outFile << "</html>" << std::endl;
  outFile.close();
  
  std::cout << std::endl;
  std::cout <<"End of test log." << std::endl;
}

void InclAblaTestHtmlLogWriter::beginTestSet(TString name)
{
  // Print some header information at the beginning of each test set output.

  outFile << "<br> <br>" << std::endl;
  outFile << "<h2>" << name << "</h2>" << std::endl;

  outFile << "<table border>" << std:: endl;
  outFile << "<caption>" << name << "</caption>" << std::endl;
  outFile <<"<tr><td><b>C++ Name</b></td><td><b>FORTRAN Name</b></td><td><b>Test Status</b></td></tr>" << std::endl;
  
  std::cout <<"Test set: " << name << std::endl;
}

void InclAblaTestHtmlLogWriter::endTestSet()
{
  // Information at the end of a test set.

  outFile << "</table>" << std::endl;
  
  std::cout <<"End of test set" << std::endl;
}

void InclAblaTestHtmlLogWriter::logTest(InclAblaTest *aTest)
{
  // Writes the information of a test to the log.

  TString *statusString;
  
  if(aTest->getTestStatus()) {
    statusString = new TString("Pass");
  }
  else {
    statusString = new TString("Fail");
  }

  outFile << "<tr><td>" << aTest->getUnitName() << "</td><td>" << aTest->getOriginalUnitName() << "</td><td><a href=\"./" << aTest->GetName() <<".html#" << aTest->GetName() << ":description" <<"\">" << (*statusString) << "</a></td></tr>" << std::endl;
  
  std::cout << aTest->getUnitName() << " \t " << (*statusString) << std::endl;
}
