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
// $Id: G4ProcessTest.hh,v 1.1 2001-10-15 12:33:20 pia Exp $
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
// Test DoIt method of physics processes
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4PROCESSTEST_HH
#define G4PROCESSTEST_HH 1

#include "globals.hh"

class G4VProcess;
class G4Track;
class G4Step;

class G4ProcessTest {
 
public:

  G4ProcessTest();

  virtual ~G4ProcessTest();
 
  void postStepTest(G4VProcess* process, 
		    const G4Track& track,
		    const G4Step& step) const;

  void alongStepTest(const G4VProcess* process, 
		     const G4Track& track,
		     const G4Step& step) const;

private:
  
  // Hide copy constructor and assignment operator
  G4ProcessTest(const G4ProcessTest&);
  G4ProcessTest & operator=(const G4ProcessTest &right);

};
 
#endif
