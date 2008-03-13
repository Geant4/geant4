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
// $Id: pyG4EmCalculator.cc,v 1.7 2008-03-13 07:32:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4EmCalculator.cc
//
//                                         2006 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4EmCalculator.hh"
#include "G4Region.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4EmCalculator {

// GetDEDX
G4double (G4EmCalculator::*f1_GetDEDX)
  (G4double, const G4ParticleDefinition*, const G4Material*, const G4Region*)
  = &G4EmCalculator::GetDEDX;

G4double (G4EmCalculator::*f2_GetDEDX)
  (G4double, const G4String&, const G4String&, const G4String&)
  = &G4EmCalculator::GetDEDX;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetDEDX, GetDEDX, 3, 4);

// GetRange
G4double (G4EmCalculator::*f1_GetRange)
  (G4double, const G4ParticleDefinition*, const G4Material*, const G4Region*)
  = &G4EmCalculator::GetRange;

G4double (G4EmCalculator::*f2_GetRange)
  (G4double, const G4String&, const G4String&, const G4String&)
  = &G4EmCalculator::GetRange;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetRange, GetRange, 3, 4);

// GetKinEnergy
G4double (G4EmCalculator::*f1_GetKinEnergy)
  (G4double, const G4ParticleDefinition*, const G4Material*, const G4Region*)
  = &G4EmCalculator::GetKinEnergy;

G4double (G4EmCalculator::*f2_GetKinEnergy)
  (G4double, const G4String&, const G4String&, const G4String&)
  = &G4EmCalculator::GetKinEnergy;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetKinEnergy, GetKinEnergy, 3, 4);

// GetCrossSectionPerVolume
G4double (G4EmCalculator::*f1_GetCrossSectionPerVolume)
  (G4double, const G4ParticleDefinition*, 
   const G4String&, const G4Material*, const G4Region*)
  = &G4EmCalculator::GetCrossSectionPerVolume;

G4double (G4EmCalculator::*f2_GetCrossSectionPerVolume)
  (G4double, const G4String&, const G4String&, 
   const G4String&, const G4String&)
  = &G4EmCalculator::GetCrossSectionPerVolume;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetCrossSectionPerVolume, 
				       GetCrossSectionPerVolume, 4, 5);

// GetCrossSectionPerAtom
#if G4VERSION_NUMBER <= 801
G4double (G4EmCalculator::*f1_GetCrossSectionPerAtom)
  (G4double, const G4ParticleDefinition*, 
   const G4String&, const G4Material*, const G4Region*)
  = &G4EmCalculator::GetCrossSectionPerAtom;

G4double (G4EmCalculator::*f2_GetCrossSectionPerAtom)
  (G4double, const G4String&, const G4String&, 
   const G4String&, const G4String&)
  = &G4EmCalculator::GetCrossSectionPerAtom;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetCrossSectionPerAtom, 
				       GetCrossSectionPerAtom, 4, 5);
#endif

// GetMeanFreePath
G4double (G4EmCalculator::*f1_GetMeanFreePath)
  (G4double, const G4ParticleDefinition*, 
   const G4String&, const G4Material*, const G4Region*)
  = &G4EmCalculator::GetMeanFreePath;

G4double (G4EmCalculator::*f2_GetMeanFreePath)
  (G4double, const G4String&, const G4String&, 
   const G4String&, const G4String&)
  = &G4EmCalculator::GetMeanFreePath;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetMeanFreePath, 
				       GetMeanFreePath, 4, 5);

// ComputeDEDX
G4double (G4EmCalculator::*f1_ComputeDEDX)
  (G4double, const G4ParticleDefinition*, 
   const G4String&, const G4Material*, G4double)
  = &G4EmCalculator::ComputeDEDX;

G4double (G4EmCalculator::*f2_ComputeDEDX)
  (G4double, const G4String&, const G4String&, const G4String&, G4double)
  = &G4EmCalculator::ComputeDEDX;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ComputeDEDX, ComputeDEDX, 4, 5);

// ComputeNuclearDEDX
#if G4VERSION_NUMBER >=710
G4double (G4EmCalculator::*f1_ComputeNuclearDEDX)
  (G4double, const G4ParticleDefinition*, const G4Material*)
  = &G4EmCalculator::ComputeNuclearDEDX;

G4double (G4EmCalculator::*f2_ComputeNuclearDEDX)
  (G4double, const G4String&, const G4String&)
  = &G4EmCalculator::ComputeNuclearDEDX;
#endif

#if G4VERSION_NUMBER >= 810
// ComputeElectronicDEDX
G4double (G4EmCalculator::*f1_ComputeElectronicDEDX)
  (G4double, const G4ParticleDefinition*, const G4Material*, G4double)
  = &G4EmCalculator::ComputeElectronicDEDX;

G4double (G4EmCalculator::*f2_ComputeElectronicDEDX)
  (G4double, const G4String&, const G4String&, G4double)
  = &G4EmCalculator::ComputeElectronicDEDX;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ComputeElectronicDEDX, 
                                       ComputeElectronicDEDX, 3, 4);

// ComputeTotalDEDX
G4double (G4EmCalculator::*f1_ComputeTotalDEDX)
  (G4double, const G4ParticleDefinition*, const G4Material*, G4double)
  = &G4EmCalculator::ComputeTotalDEDX;

G4double (G4EmCalculator::*f2_ComputeTotalDEDX)
  (G4double, const G4String&, const G4String&, G4double)
  = &G4EmCalculator::ComputeTotalDEDX;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ComputeTotalDEDX, ComputeTotalDEDX, 
                                       3, 4);
#endif


// ComputeCrossSectionPerVolume
G4double (G4EmCalculator::*f1_ComputeCrossSectionPerVolume)
  (G4double, const G4ParticleDefinition*, 
   const G4String&, const G4Material*, G4double)
  = &G4EmCalculator::ComputeCrossSectionPerVolume;

G4double (G4EmCalculator::*f2_ComputeCrossSectionPerVolume)
  (G4double, const G4String&, const G4String&, const G4String&, G4double)
  = &G4EmCalculator::ComputeCrossSectionPerVolume;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ComputeCrossSectionPerVolume, 
				       ComputeCrossSectionPerVolume, 4, 5);

// ComputeCrossSectionPerAtom
G4double (G4EmCalculator::*f1_ComputeCrossSectionPerAtom)
  (G4double, const G4ParticleDefinition*, const G4String&, 
   G4double, G4double, G4double)
  = &G4EmCalculator::ComputeCrossSectionPerAtom;

G4double (G4EmCalculator::*f2_ComputeCrossSectionPerAtom)
  (G4double, const G4String&, const G4String&, const G4Element*, G4double)
  = &G4EmCalculator::ComputeCrossSectionPerAtom;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ComputeCrossSectionPerAtom, 
				       ComputeCrossSectionPerAtom, 5, 6);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(g_ComputeCrossSectionPerAtom, 
				       ComputeCrossSectionPerAtom, 4, 5);

#if G4VERSION_NUMBER >= 830
// ComputeEnergyCutFromRangeCut
G4double (G4EmCalculator::*f1_ComputeEnergyCutFromRangeCut)
  (G4double, const G4ParticleDefinition*, const G4Material*)
  = &G4EmCalculator::ComputeEnergyCutFromRangeCut;

G4double (G4EmCalculator::*f2_ComputeEnergyCutFromRangeCut)
  (G4double range, const G4String&, const G4String&)
  = &G4EmCalculator::ComputeEnergyCutFromRangeCut;
#endif


// ComputeMeanFreePath
G4double (G4EmCalculator::*f1_ComputeMeanFreePath)
  (G4double, const G4ParticleDefinition*, 
   const G4String&, const G4Material*, G4double)
  = &G4EmCalculator::ComputeMeanFreePath;

G4double (G4EmCalculator::*f2_ComputeMeanFreePath)
  (G4double, const G4String&, const G4String&, const G4String&, G4double)
  = &G4EmCalculator::ComputeMeanFreePath;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ComputeMeanFreePath, 
				       ComputeMeanFreePath, 4, 5);

// FindCouple
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_FindCouple, FindCouple, 1, 2);

};

using namespace pyG4EmCalculator;

// ====================================================================
// module definition
// ====================================================================
void export_G4EmCalculator()
{
  class_<G4EmCalculator, boost::noncopyable>
    ("G4EmCalculator", "Provide access to dE/dx and cross section")
    // ---
    .def("GetDEDX",       f1_GetDEDX,       f_GetDEDX())
    .def("GetDEDX",       f2_GetDEDX,       f_GetDEDX())
    .def("GetRange",      f1_GetRange,      f_GetRange())
    .def("GetRange",      f2_GetDEDX,       f_GetRange())
    .def("GetKinEnergy",  f1_GetKinEnergy,  f_GetKinEnergy())
    .def("GetKinEnergy",  f2_GetKinEnergy,  f_GetKinEnergy())
    .def("GetCrossSectionPerVolume",  
	 f1_GetCrossSectionPerVolume, f_GetCrossSectionPerVolume())
    .def("GetCrossSectionPerVolume",  
	 f2_GetCrossSectionPerVolume, f_GetCrossSectionPerVolume())
#if G4VERSION_NUMBER <= 801
    .def("GetCrossSectionPerAtom",  
         f1_GetCrossSectionPerAtom, f_GetCrossSectionPerAtom())
    .def("GetCrossSectionPerAtom",  
         f2_GetCrossSectionPerAtom, f_GetCrossSectionPerAtom())
#endif
    .def("GetMeanFreePath",  f1_GetMeanFreePath,  f_GetMeanFreePath())
    .def("GetMeanFreePath",  f2_GetMeanFreePath,  f_GetMeanFreePath())
    // ---
    .def("PrintDEDXTable",         &G4EmCalculator::PrintDEDXTable)
    .def("PrintRangeTable",        &G4EmCalculator::PrintRangeTable)
    .def("PrintInverseRangeTable", &G4EmCalculator::PrintInverseRangeTable)
    // ---
    .def("ComputeDEDX",            f1_ComputeDEDX,  f_ComputeDEDX())
    .def("ComputeDEDX",            f2_ComputeDEDX,  f_ComputeDEDX())
#if G4VERSION_NUMBER >=710
    .def("ComputeNuclearDEDX",     f1_ComputeNuclearDEDX)
    .def("ComputeNuclearDEDX",     f2_ComputeNuclearDEDX)
#endif
#if G4VERSION_NUMBER >= 810
    .def("ComputeElectronicDEDX",  f1_ComputeElectronicDEDX,  
         f_ComputeElectronicDEDX())
    .def("ComputeDEDX",            f2_ComputeElectronicDEDX,  
         f_ComputeElectronicDEDX())
    .def("ComputeTotalDEDX",       f1_ComputeTotalDEDX,  f_ComputeTotalDEDX())
    .def("ComputeTotalDEDX",       f2_ComputeTotalDEDX,  f_ComputeTotalDEDX())
#endif
    // ---
    .def("ComputeCrossSectionPerVolume",
	 f1_ComputeCrossSectionPerVolume, f_ComputeCrossSectionPerVolume())
    .def("ComputeCrossSectionPerVolume",
	 f2_ComputeCrossSectionPerVolume, f_ComputeCrossSectionPerVolume())
    .def("ComputeCrossSectionPerAtom",
         f1_ComputeCrossSectionPerAtom, f_ComputeCrossSectionPerAtom())
    .def("ComputeCrossSectionPerAtom",
         f2_ComputeCrossSectionPerAtom, g_ComputeCrossSectionPerAtom())
#if G4VERSION_NUMBER >= 830
    .def("ComputeEnergyCutFromRangeCut", f1_ComputeEnergyCutFromRangeCut)
    .def("ComputeEnergyCutFromRangeCut", f2_ComputeEnergyCutFromRangeCut)
#endif
    // ---
    .def("ComputeMeanFreePath", 
	 f1_ComputeMeanFreePath, f_ComputeMeanFreePath())
    .def("ComputeMeanFreePath",
	 f2_ComputeMeanFreePath, f_ComputeMeanFreePath())
    // ---
    .def("FindParticle",  &G4EmCalculator::FindParticle,
         return_value_policy<reference_existing_object>())
    .def("FindMaterial",  &G4EmCalculator::FindMaterial,
         return_value_policy<reference_existing_object>())
    .def("FindRegion",    &G4EmCalculator::FindRegion,
         return_value_policy<reference_existing_object>())
    .def("FindCouple",    &G4EmCalculator::FindCouple,
	 f_FindCouple()[return_value_policy<reference_existing_object>()])
    // ---
    .def("SetVerbose",    &G4EmCalculator::SetVerbose)
    ;
}

