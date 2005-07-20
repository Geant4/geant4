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
// $Id: G4DNAGenericIonsManager.cc,v 1.1 2005-07-20 10:01:00 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4DNAGenericIonsManager.hh"
#include "G4Alpha.hh"
#include "G4VIon.hh"

G4DNAGenericIonsManager *                G4DNAGenericIonsManager :: Instance(void)
{
 if (!theInstance)
  theInstance=new G4DNAGenericIonsManager;
 
 return theInstance;
}

G4ParticleDefinition *                   G4DNAGenericIonsManager :: GetIon(const G4String & name)
{
 IonsMap::const_iterator i(map.find(name));
 
 if (i==map.end())
  return 0;
  
// return map[name];

 return i->second;
}

                                         G4DNAGenericIonsManager :: G4DNAGenericIonsManager()
{
 //               name             mass          width         charge
 //             2*spin           parity  C-conjugation
 //          2*Isospin       2*Isospin3       G-parity
 //               type    lepton number  baryon number   PDG encoding
 //             stable         lifetime    decay table


 G4VIon *helium;
 G4VIon *hydrogen;
 G4VIon *alphaPlus;

 helium=     new G4VIon("helium",    3.727417*GeV,       0.0*MeV,  +0.0*eplus,
                               0,              +1,             0,
                               0,               0,             0,
                       "nucleus",              +2,            +4,           0,
                            true,            -1.0,           NULL);

 alphaPlus=  new G4VIon("alpha+",    3.727417*GeV,       0.0*MeV,  +1.0*eplus,
                               1,              +1,             0,
                               0,               0,             0,
                       "nucleus",              +1,            +4,           0,
                            true,            -1.0,           NULL);

 hydrogen= new G4VIon("hydrogen",   0.9382723*GeV,       0.0*MeV,  +0.0*eplus,
                               0,              +1,             0,
                               0,               0,             0,
                       "nucleus",              +1,            +1,           0,
                            true,            -1.0,           NULL);

 map["helium"]=helium;
 map["hydrogen"]=hydrogen;
 map["alpha+"]=alphaPlus;
 map["alpha++"]=G4Alpha::Alpha();
}


G4DNAGenericIonsManager *                G4DNAGenericIonsManager::theInstance(0);
   
