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
// $Id: Tst33WeightWindowStoreBuilder.hh,v 1.1 2003-08-15 15:35:54 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33WeightWindowStoreBuilder
//
// Class description:
// Creats the importance store and sets the importance
// values. 
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33WeightWindowStoreBuilder_hh
#define Tst33WeightWindowStoreBuilder_hh Tst33WeightWindowStoreBuilder_hh

class G4VWeightWindowStore;
class Tst33VGeometry;

class Tst33WeightWindowStoreBuilder {
public:
  Tst33WeightWindowStoreBuilder();
  ~Tst33WeightWindowStoreBuilder();
  G4VWeightWindowStore *CreateWeightWindowStore(Tst33VGeometry *samplegeo);  
};


#endif
