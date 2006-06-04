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
// $Id: pymodG4particles.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4particles.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4ParticleDefinition();
void export_G4DynamicParticle();
void export_G4ParticleTable();
void export_G4DecayTable();
void export_G4PrimaryParticle();
void export_G4PrimaryVertex();
void export_PyG4ParticleList();


BOOST_PYTHON_MODULE(G4particles)
{
  export_G4ParticleDefinition();
  export_G4DynamicParticle();
  export_G4ParticleTable();
  export_G4DecayTable();
  export_G4PrimaryParticle();
  export_G4PrimaryVertex();
  export_PyG4ParticleList();
}

