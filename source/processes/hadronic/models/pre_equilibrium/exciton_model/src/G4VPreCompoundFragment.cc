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
//
// $Id: G4VPreCompoundFragment.cc,v 1.6 2008-05-08 10:42:35 quesada Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara
//
//J. M. Quesada (May. 2008) . Coulomb barrier has been kept  (NOW implicitely included through cross sections)
// because of inheritance (just in case..) by HETC classes. They are not activeted by default.
#include "G4VPreCompoundFragment.hh"
#include "G4PreCompoundParameters.hh"

G4VPreCompoundFragment::
G4VPreCompoundFragment(const G4VPreCompoundFragment & right)
{
  theA = right.theA;
  theZ = right.theZ;
  theRestNucleusA = right.theRestNucleusA;
  theRestNucleusZ = right.theRestNucleusZ;
  theCoulombBarrier = right.theCoulombBarrier;
  theCoulombBarrierPtr = right.theCoulombBarrierPtr;
  theMaximalKineticEnergy = right.theMaximalKineticEnergy;
  theEmissionProbability = right.theEmissionProbability;
  theMomentum = right.theMomentum;
  theFragmentName = right.theFragmentName;
  theStage = right.theStage;
}

G4VPreCompoundFragment::
G4VPreCompoundFragment(const G4double anA,
		       const G4double aZ, 
		       G4VCoulombBarrier* aCoulombBarrier,
		       const G4String & aName):
  theA(anA),theZ(aZ), 
  theRestNucleusA(0.0),theRestNucleusZ(0.0),theCoulombBarrier(0.0),
  theCoulombBarrierPtr(aCoulombBarrier),
  theBindingEnergy(0.0), theMaximalKineticEnergy(-1.0),
  theEmissionProbability(0.0), theMomentum(0.0,0.0,0.0,0.0),
  theFragmentName(aName),theStage(0),G4Prec(false)
{}

//JMQ (May 08): new  constructor to be used by new G4PrecompoundFragment class (without coulomb barrier)


G4VPreCompoundFragment::
G4VPreCompoundFragment(const G4double anA,
		       const G4double aZ, 		       
		       const G4String & aName):
  theA(anA),theZ(aZ), 
  theRestNucleusA(0.0),theRestNucleusZ(0.0),
  theBindingEnergy(0.0), theMaximalKineticEnergy(-1.0),
  theEmissionProbability(0.0), theMomentum(0.0,0.0,0.0,0.0),
  theFragmentName(aName),theStage(0),G4Prec(true)
{}

//JMQ end of new constructor 

G4VPreCompoundFragment::~G4VPreCompoundFragment()
{
}


const G4VPreCompoundFragment & G4VPreCompoundFragment::
operator= (const G4VPreCompoundFragment & right)
{
  if (this != &right) {
    theA = right.theA;
    theZ = right.theZ;
    theRestNucleusA = right.theRestNucleusA;
    theRestNucleusZ = right.theRestNucleusZ;
    theCoulombBarrier = right.theCoulombBarrier;
    theCoulombBarrierPtr = right.theCoulombBarrierPtr;
    theMaximalKineticEnergy = right.theMaximalKineticEnergy;
    theEmissionProbability = right.theEmissionProbability;
    theMomentum = right.theMomentum;
    theFragmentName = right.theFragmentName;
    theStage = right.theStage;
  }
  return *this;
}

G4int G4VPreCompoundFragment::operator==(const G4VPreCompoundFragment & right) const
{
  return (this == (G4VPreCompoundFragment *) &right);
}

G4int G4VPreCompoundFragment::operator!=(const G4VPreCompoundFragment & right) const
{
  return (this != (G4VPreCompoundFragment *) &right);
}


std::ostream& 
operator << (std::ostream &out, const G4VPreCompoundFragment &theFragment)
{
  out << &theFragment;
  return out; 
}


std::ostream& 
operator << (std::ostream &out, const G4VPreCompoundFragment *theFragment)
{
  std::ios::fmtflags old_floatfield = out.flags();
  out.setf(std::ios::floatfield);
    
  out 
    << "PreCompound Model Emitted Fragment: A = " 
    << std::setprecision(3) << theFragment->theA 
    << ", Z = " << std::setprecision(3) << theFragment->theZ;
    out.setf(std::ios::scientific,std::ios::floatfield);
    //   out
    //     << ", U = " << theFragment->theExcitationEnergy/MeV 
    //     << " MeV" << endl
    //     << "          P = (" 
    //     << theFragment->theMomentum.x()/MeV << ","
    //     << theFragment->theMomentum.y()/MeV << ","
    //     << theFragment->theMomentum.z()/MeV 
    //     << ") MeV   E = " 
    //     << theFragment->theMomentum.t()/MeV << " MeV";
    
    out.setf(old_floatfield,std::ios::floatfield);
    
    return out;
}


void G4VPreCompoundFragment::
Initialize(const G4Fragment & aFragment)
{
  theRestNucleusA = aFragment.GetA() - theA;
  theRestNucleusZ = aFragment.GetZ() - theZ;

  if ((theRestNucleusA < theRestNucleusZ) ||
      (theRestNucleusA < theA) ||
      (theRestNucleusZ < theZ)) 
    {
      // In order to be sure that emission probability will be 0.
      theMaximalKineticEnergy = 0.0;
      return;
    }
  
 //JMQ New "if" control structure to conform  to the case when no Coumlomb barrier is considered (ours).  
  if (!G4Prec){
  // Calculate Coulomb barrier
    theCoulombBarrier = theCoulombBarrierPtr->
    GetCoulombBarrier(static_cast<G4int>(theRestNucleusA),static_cast<G4int>(theRestNucleusZ),
    aFragment.GetExcitationEnergy());}
  
  // Compute Binding Energies for fragments 
  // (needed to separate a fragment from the nucleus)
  
  theBindingEnergy = G4NucleiProperties::GetMassExcess(static_cast<G4int>(theA),static_cast<G4int>(theZ)) +
    G4NucleiProperties::GetMassExcess(static_cast<G4int>(theRestNucleusA),static_cast<G4int>(theRestNucleusZ)) -
    G4NucleiProperties::GetMassExcess(static_cast<G4int>(aFragment.GetA()),static_cast<G4int>(aFragment.GetZ()));
  
  // Compute Maximal Kinetic Energy which can be carried by fragments after separation
  G4double m = aFragment.GetMomentum().m();
  G4double rm = GetRestNuclearMass();
  G4double em = GetNuclearMass();
  theMaximalKineticEnergy = ((m - rm)*(m + rm) + em*em)/(2.0*m) - em;
//JMQ : this (previously) commented line expresses doubts of previous coder. Trivial matter. Warning! 
  // - theCoulombBarrier;
  
  return;
}




