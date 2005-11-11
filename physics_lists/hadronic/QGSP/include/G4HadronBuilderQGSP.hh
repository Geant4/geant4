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
// $Id: G4HadronBuilderQGSP.hh,v 1.1 2005-11-11 22:57:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronBuilderQGSP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard 
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronBuilderQGSP_h
#define G4HadronBuilderQGSP_h 1

#include "globals.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4HadronQEDBuilder.hh"
#include "G4StoppingHadronBuilder.hh"
#include "G4MiscLHEPBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4LEPPiKBuilder.hh"
#include "G4QGSPPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4LEPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4LEPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"

class G4NeutronBuilder;
class G4LEPNeutronBuilder;
class G4QGSPNeutronBuilder;
class G4PiKBuilder;
class G4LEPPiKBuilder;
class G4QGSPPiKBuilder;
class G4ProtonBuilder;
class G4LEPProtonBuilder;
class G4QGSPProtonBuilder; 
class G4MiscLHEPBuilder;
class G4StoppingHadronBuilder;

class G4HadronBuilderQGSP : public G4VPhysicsConstructor
{
public: 
  G4HadronBuilderQGSP(const G4String& name ="QGSP");
  virtual ~G4HadronBuilderQGSP();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:

  G4bool wasActivated;

  G4NeutronBuilder* theNeutrons;
  G4LEPNeutronBuilder* theLEPNeutron;
  G4QGSPNeutronBuilder* theQGSPNeutron;
  
  G4PiKBuilder* thePiK;
  G4LEPPiKBuilder* theLEPPiK;
  G4QGSPPiKBuilder* theQGSPPiK;
    
  G4ProtonBuilder* thePro;
  G4LEPProtonBuilder* theLEPPro;
  G4QGSPProtonBuilder* theQGSPPro;    
    
  G4MiscLHEPBuilder* theMiscLHEP;
  G4StoppingHadronBuilder* theStoppingHadron;
};

#endif

