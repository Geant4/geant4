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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              Nucleus.cc
//
// Version:             0.b.3
// Date:                29/02/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
//
#include "Nucleus.hh"
////////////////////////////////////////////////////////////////////////////////
//
Nucleus::Nucleus ()
  : a(24), z(11), e(0.0)
{;}
///////////////////////////////////////////////////////////////////////////////
//
Nucleus::Nucleus (G4int a1, G4int z1, G4double e1)
{
  //
  //
  a = a1;
  z = z1;
  e = e1;
}
///////////////////////////////////////////////////////////////////////////////
//
Nucleus::~Nucleus ()
{;}
///////////////////////////////////////////////////////////////////////////////
//
G4std::ostream &operator << (G4std::ostream &s, const Nucleus &q)
//
//
// Definition of the insertion operator << to provide the nucleus limits to
// ostream.
//
{
  s <<"Atomic weight: " <<q.GetA()
    <<"Atomic number: " <<q.GetZ()
    <<"Excitation energy: "<<q.GetE();
  return s;
}
///////////////////////////////////////////////////////////////////////////////






