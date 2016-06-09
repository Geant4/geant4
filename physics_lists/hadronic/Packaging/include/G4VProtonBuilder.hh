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
// $Id: G4VProtonBuilder.hh,v 1.2 2005/11/25 15:38:50 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4VProtonBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 21.11.2005 G.Folger: don't  keep processes as data members, but new these
//
//----------------------------------------------------------------------------
//
#ifndef G4VProtonBuilder_h
#define G4VProtonBuilder_h

class G4ProtonInelasticProcess;
class G4HadronElasticProcess;

class G4VProtonBuilder
{
  public:
    G4VProtonBuilder() {}
    virtual ~G4VProtonBuilder() {}
    virtual void Build(G4HadronElasticProcess * aP) = 0;
    virtual void Build(G4ProtonInelasticProcess * aP) = 0;
};
// 2002 by J.P. Wellisch

#endif
