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
// $Id: Tst33IStoreBuilder.hh,v 1.3 2002-11-20 09:38:25 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33IStoreBuilder
//
// Class description:
// Creats the importance store and sets the importance
// values. 
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33IStoreBuilder_hh
#define Tst33IStoreBuilder_hh Tst33IStoreBuilder_hh

class G4VIStore;
class Tst33VGeometry;

class Tst33IStoreBuilder {
public:
  Tst33IStoreBuilder();
  ~Tst33IStoreBuilder();
  G4VIStore *CreateIStore(Tst33VGeometry *samplegeo);  
};


#endif
