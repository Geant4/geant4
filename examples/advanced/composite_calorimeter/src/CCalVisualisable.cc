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
///////////////////////////////////////////////////////////////////////////////
// File: CCalVisualisable.cc
// Description: Sets visualisable attributes from .vis files
///////////////////////////////////////////////////////////////////////////////
#include "CCalVisualisable.hh"

#include <fstream>
#include "CCalutils.hh"

//Comment/Uncomment next line to unset/set debug information printing
//#define debug

const char* visEnvName = "CCAL_VISPATH";
const char* CCalVisualisable::pathName=0;

CCalVisualisable::CCalVisualisable(G4String file):visFile(file) {
  if (!pathName)
    setPath();
  setDefault();
  readFile();
}

G4bool CCalVisualisable::readFile(G4String file) { 
  visFile=file;
  return readFile();
}

void CCalVisualisable::setDefault(){
  theParameters[Sensitive] = visParameters(true);
  theParameters[Electronics] = visParameters(false);
  theParameters[Support] = visParameters(false);
  theParameters[Cable] = visParameters(false);
  theParameters[OtherServices] = visParameters(false);
  theParameters[PseudoVolumes] = visParameters(false);
}

void CCalVisualisable::setColor(visType v, G4double r, G4double g, G4double b){
  theParameters[v].rColor=r;
  theParameters[v].gColor=g;
  theParameters[v].bColor=b;
}

void CCalVisualisable::setPath() {
  pathName = std::getenv(visEnvName);
  if (!pathName) {
     G4ExceptionDescription ed;
     ed << "ERROR: " << visEnvName << " environmental variable not set!" 
        << G4endl;
     ed << "       Set it and restart." << G4endl;
     G4Exception("CCalVisualisable::setPath()","ccal007",
                 FatalException,ed);
  }
}

G4bool CCalVisualisable::readFile() {
  if (visFile=="") {
    G4cerr << "ERROR: No file was specified from which to read Visualisation parameters" 
         << G4endl;
    return false;
  }

  ///////////////////////////////////////////////////////////////
  //Let's open the file
  G4cout << " ==> Opening file " << visFile 
       << " to read visualisation parameters..."
       << G4endl;

  G4String pathname(pathName);
  std::ifstream is;
#ifdef debug
  G4cout << "Viualisable : Path " << pathname << " FIle " << visFile << G4endl;
#endif
  G4bool ok = openGeomFile(is, pathname, visFile);
  if (!ok) {
    G4cout << "WARNING: Could not read " << visFile << G4endl;
    G4cout << "         Default visualization parameters will be used." << G4endl;
    return false;
  }
  
  while (is) {
    G4String viewvol;
    readName(is,viewvol);

    visType vt = Undefined;
    if (viewvol=="Sensitive")
      vt=Sensitive;
    else if (viewvol=="Electronics")
      vt=Electronics;
    else if (viewvol=="Support")
      vt=Support;
    else if (viewvol=="Absorber")
      vt=Absorber;
    else if (viewvol=="Cable")
      vt=Cable;
    else if (viewvol=="OtherServices")
      vt=OtherServices;
    else if (viewvol=="PseudoVolumes")
      vt=PseudoVolumes;
    else {
      vt=Undefined;
      G4cerr << "WARNING: Unknown type of visualisable object \"" << viewvol
           << "\"." << G4endl;
    }


    G4int isvisible, wireframe;
    G4double r, g, b;
    
    is >> isvisible >> r >> g >> b >> wireframe >> jump;

    r=checkColorRange(r,'R');
    g=checkColorRange(g,'G');
    b=checkColorRange(b,'B');

    if (vt!=Undefined) {
#ifdef debug
      G4cout << tab << viewvol << tab << isvisible << tab 
           << r << " " << g << " "<< b << tab
           << wireframe << G4endl;
#endif
      theParameters[vt]=visParameters(isvisible, r, g, b, wireframe);
    }

  }

  ///////////////////////////////////////////////////////////////
  // Close the file
  G4cout << " ==> Closing file " << visFile << G4endl;
  is.close();

  return true;
}

G4double CCalVisualisable::checkColorRange(G4double cvalue, char ctype) const {
  if (cvalue>1) {
    G4cerr << "ERROR: In " << visFile << ". Color " << ctype << "=" 
         << cvalue << " > 1" << G4endl;
    G4cerr << "       It will be reset to 1." << G4endl;
    return 1.;
  }
  if (cvalue<0) {
    G4cerr << "ERROR: In " << visFile << ". Color " << ctype << "=" 
         << cvalue << " < 0" << G4endl;
    G4cerr << "       It will be reset to 0." << G4endl;
    return 0.;
  }
  return cvalue;
}
