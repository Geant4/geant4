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
// $Id: G4ShortLivedConstructor.hh,v 1.4 2001-07-11 10:02:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
#ifndef G4ShortLivedConstructor_h
#define G4ShortLivedConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4ShortLivedConstructor
{
  //This class is a utility class for construction 
  //short lived particles

  public:
    G4ShortLivedConstructor();
    ~G4ShortLivedConstructor();
  
  public:
    void ConstructParticle();
 
  protected:
    void ConstructResonances();
    void ConstructBaryons();
    void ConstructMesons();
    void ConstructQuarks();

  private:
    static G4bool isConstructed;
    // flag for checking whether resonces exist or not
};

#endif











