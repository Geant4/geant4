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
//
// $Id: G4ComptonTest.hh,v 1.1 2001-10-29 09:28:53 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Test of electromagnetic physics processes
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4COMPTONTEST_HH
#define G4COMPTONTEST_HH 1

#include "globals.hh"
#include "G4ProcessTest.hh"

class G4VProcess;

class G4ComptonTest : public G4ProcessTest 
{
public:

  G4ComptonTest(const G4String& category, G4bool isPolarised);
  virtual ~G4ComptonTest(); 

  protected:

  virtual G4VProcess* createProcess();
  virtual G4VProcess* createBremsstrahlung();
  virtual G4VProcess* createElectronIonisation();

private:
  
  // Hide copy constructor and assignment operator
  G4ComptonTest(const G4ComptonTest&);
  G4ComptonTest & operator=(const G4ComptonTest &right);

  const G4String& type ;
  G4bool polarised;

};
 
#endif



