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
// $Id: G4KaonBuilder.hh,v 1.1 2010-11-19 16:21:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4KaonBuilder
//
// Author: 2010 G.Folger
//  devired from G4PiKBuilder
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
// 13.06.2006 G.Folger: (re)move elastic scatterring 
//
//----------------------------------------------------------------------------
//
#ifndef G4KaonBuilder_h
#define G4KaonBuilder_h 1

#include "globals.hh"

#include "G4ProtonInelasticProcess.hh"
#include "G4VKaonBuilder.hh"

#include <vector>

class G4KaonBuilder
{
  public: 
    G4KaonBuilder();
    virtual ~G4KaonBuilder();

  public: 
    void Build();
    void RegisterMe(G4VKaonBuilder * aB) {theModelCollections.push_back(aB);}

  private:
    G4KaonPlusInelasticProcess*  theKaonPlusInelastic;
    G4KaonMinusInelasticProcess* theKaonMinusInelastic;
    G4KaonZeroLInelasticProcess* theKaonZeroLInelastic;
    G4KaonZeroSInelasticProcess* theKaonZeroSInelastic;
     
    std::vector<G4VKaonBuilder *> theModelCollections;

    G4bool wasActivated;
};
#endif

