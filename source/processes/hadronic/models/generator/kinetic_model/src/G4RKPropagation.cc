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
// -------------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4RKPropagation.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 6 June 2000
// -------------------------------------------------------------------
#include "G4RKPropagation.hh"
// nuclear fields
#include "G4VNuclearField.hh"
#include "G4ProtonField.hh"
#include "G4NeutronField.hh"
#include "G4AntiProtonField.hh"
#include "G4KaonPlusField.hh"
#include "G4KaonMinusField.hh"
#include "G4KaonZeroField.hh"
#include "G4PionPlusField.hh"
#include "G4PionMinusField.hh"
#include "G4PionZeroField.hh"
#include "G4SigmaPlusField.hh"
#include "G4SigmaMinusField.hh"
#include "G4SigmaZeroField.hh"
// particles properties
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4AntiProton.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaMinus.hh"
#include "G4SigmaZero.hh"

#include "globals.hh"

#include "G4KM_OpticalEqRhs.hh"
#include "G4KM_NucleonEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorDriver.hh"

// unsigned EncodingHashFun(const G4int& aEncoding);

G4RKPropagation::G4RKPropagation() : theNucleus(0), 
				     theFieldMap(0), theEquationMap(0),
				     theField(0)
{ }


G4RKPropagation::G4RKPropagation(const  G4RKPropagation &right)
{ }


G4RKPropagation::~G4RKPropagation()
{
// free theFieldMap memory
  if(theFieldMap) delete_FieldsAndMap(theFieldMap);

// free theEquationMap memory
  if(theEquationMap) delete_EquationsAndMap(theEquationMap);

  if (theField) delete theField;
}



const G4RKPropagation & G4RKPropagation::operator=(const G4RKPropagation & right)
{
  G4Exception("G4RKPropagation::operator= meant not to be accessible");
  return *this;
}

G4int G4RKPropagation::operator==(const G4RKPropagation & right) const
{
  G4Exception("G4RKPropagation::operator== meant not to be accessible");
  return 0;
}

G4int G4RKPropagation::operator!=(const G4RKPropagation & right) const
{
  G4Exception("G4RKPropagation::operator!= meant not to be accessible");
  return 1;
}

//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
void G4RKPropagation::Init(G4V3DNucleus * nucleus)
//----------------------------------------------------------------------------
{

// free theFieldMap memory
  if(theFieldMap) delete_FieldsAndMap(theFieldMap);

// free theEquationMap memory
  if(theEquationMap) delete_EquationsAndMap(theEquationMap);

  if (theField) delete theField;

// Initialize the nuclear field map.
  theNucleus = nucleus;
  theOuterRadius = theNucleus->GetOuterRadius();

  theFieldMap = new G4std::map <G4int, G4VNuclearField*, G4std::less<G4int> >;

  (*theFieldMap)[G4Proton::Proton()->GetPDGEncoding()] = new G4ProtonField(theNucleus);
  (*theFieldMap)[G4Neutron::Neutron()->GetPDGEncoding()] = new G4NeutronField(theNucleus);
  (*theFieldMap)[G4AntiProton::AntiProton()->GetPDGEncoding()] = new G4AntiProtonField(theNucleus);
  (*theFieldMap)[G4KaonPlus::KaonPlus()->GetPDGEncoding()] = new G4KaonPlusField(theNucleus);
  (*theFieldMap)[G4KaonMinus::KaonMinus()->GetPDGEncoding()] = new G4KaonMinusField(theNucleus);
  (*theFieldMap)[G4KaonZero::KaonZero()->GetPDGEncoding()] = new G4KaonZeroField(theNucleus);
  (*theFieldMap)[G4PionPlus::PionPlus()->GetPDGEncoding()] = new G4PionPlusField(theNucleus);
  (*theFieldMap)[G4PionMinus::PionMinus()->GetPDGEncoding()] = new G4PionMinusField(theNucleus);
  (*theFieldMap)[G4PionZero::PionZero()->GetPDGEncoding()] = new G4PionZeroField(theNucleus);
  (*theFieldMap)[G4SigmaPlus::SigmaPlus()->GetPDGEncoding()] = new G4SigmaPlusField(theNucleus);
  (*theFieldMap)[G4SigmaMinus::SigmaMinus()->GetPDGEncoding()] = new G4SigmaMinusField(theNucleus);
  (*theFieldMap)[G4SigmaZero::SigmaZero()->GetPDGEncoding()] = new G4SigmaZeroField(theNucleus);

  theEquationMap = new G4std::map <G4int, G4Mag_EqRhs*, G4std::less<G4int> >;

// theField needed by the design of G4Mag_eqRhs
  theField = new G4KM_DummyField;		//Field not needed for integration
  G4KM_OpticalEqRhs * opticalEq;
  G4KM_NucleonEqRhs * nucleonEq;
  G4double mass;
  G4double opticalCoeff;

  nucleonEq = new G4KM_NucleonEqRhs(theField, theNucleus);
  mass = G4Proton::Proton()->GetPDGMass();
  nucleonEq->SetMass(mass);
  (*theEquationMap)[G4Proton::Proton()->GetPDGEncoding()] = nucleonEq;

  nucleonEq = new G4KM_NucleonEqRhs(theField, theNucleus);
  mass = G4Neutron::Neutron()->GetPDGMass();
  nucleonEq->SetMass(mass);
  (*theEquationMap)[G4Neutron::Neutron()->GetPDGEncoding()] = nucleonEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4AntiProton::AntiProton()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4AntiProton::AntiProton()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4AntiProton::AntiProton()->GetPDGEncoding()] = opticalEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4KaonPlus::KaonPlus()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4KaonPlus::KaonPlus()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4KaonPlus::KaonPlus()->GetPDGEncoding()] = opticalEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4KaonMinus::KaonMinus()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4KaonMinus::KaonMinus()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4KaonMinus::KaonMinus()->GetPDGEncoding()] = opticalEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4KaonZero::KaonZero()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4KaonZero::KaonZero()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4KaonZero::KaonZero()->GetPDGEncoding()] = opticalEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4PionPlus::PionPlus()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4PionPlus::PionPlus()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4PionPlus::PionPlus()->GetPDGEncoding()] = opticalEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4PionMinus::PionMinus()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4PionMinus::PionMinus()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4PionMinus::PionMinus()->GetPDGEncoding()] = opticalEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4PionZero::PionZero()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4PionZero::PionZero()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4PionZero::PionZero()->GetPDGEncoding()] = opticalEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4SigmaPlus::SigmaPlus()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4SigmaPlus::SigmaPlus()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4SigmaPlus::SigmaPlus()->GetPDGEncoding()] = opticalEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4SigmaMinus::SigmaMinus()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4SigmaMinus::SigmaMinus()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4SigmaMinus::SigmaMinus()->GetPDGEncoding()] = opticalEq;

  opticalEq = new G4KM_OpticalEqRhs(theField, theNucleus);
  mass = G4SigmaZero::SigmaZero()->GetPDGMass();
  opticalCoeff =
    (*theFieldMap)[G4SigmaZero::SigmaZero()->GetPDGEncoding()]->GetCoeff();
  opticalEq->SetFactor(mass,opticalCoeff);
  (*theEquationMap)[G4SigmaZero::SigmaZero()->GetPDGEncoding()] = opticalEq;
}



//----------------------------------------------------------------------------
void G4RKPropagation::Transport(G4KineticTrackVector & active,
//----------------------------------------------------------------------------
				const G4KineticTrackVector & spectators,
				G4double timeStep)
{

// Loop over tracks

  G4std::vector<G4KineticTrack *>::iterator i;
  for(i = active.begin(); i != active.end(); ++i)
  {
    G4double currTimeStep = timeStep;
    G4KineticTrack * kt = *i;
    G4int encoding = kt->GetDefinition()->GetPDGEncoding();
    G4std::map <G4int, G4VNuclearField*, G4std::less<G4int> >::iterator fieldIter= theFieldMap->find(encoding);

    G4VNuclearField* currentField=0;
    if ( fieldIter != theFieldMap->end() ) currentField=fieldIter->second;

// debug
//    if ( timeStep > 1e30 ) {
//	G4cout << " Name :" << kt->GetDefinition()->GetParticleName() << G4endl;
//    }
// @hpw@ debugging: free transport..... @@@@@@@@@@@@@@@
//gf	pos = pos+(currTimeStep*c_light/mom.e())*mom.vect();
//gf	kt->SetPosition(pos);
//gf	continue;
// @hpw@ debugging: free transport.....

// Get the time of intersections with the nucleus surface.
    G4double t_enter, t_leave;
// if the particle does not intersecate with the nucleus go to next particle
    if(!GetSphereIntersectionTimes(kt, t_enter, t_leave))
      continue;

/*
 *    G4cout <<" timeStep, Intersection times tenter, tleave  "
 *    	<< currTimeStep << " / " << t_enter << " / " << t_leave <<G4endl;
 */
// if the particle is already outside nucleus go to next
    if(t_leave < 0)
      continue;

// Apply a straight line propagation for particle types
// not included in the model
    if( ! currentField )
    {
      if(currTimeStep == DBL_MAX)currTimeStep = t_leave;
      FreeTransport(kt, currTimeStep);
//      G4cout << " Particle not in model : " << kt->GetDefinition()->GetParticleName() << G4endl;
      continue;
    }

    if(t_enter > 0)  // the particle is out. Transport free to the surface
    {
      if(t_enter > currTimeStep)  // the particle won't enter the nucleus
      {
        FreeTransport(kt, currTimeStep);
	continue;
      }
      else
      {
        FreeTransport(kt, t_enter);   // go to surface
	currTimeStep -= t_enter;
	t_leave	     -= t_enter;  // time left to leave nucleus
// on the surface the particle loose the barrier energy
// 	G4double newE = mom.e()-(*theFieldMap)[encoding]->GetBarrier();
//     GetField = Barrier + FermiPotential
	G4double newE = kt->Get4Momentum().e()-currentField->GetField(kt->GetPosition());
//      G4cout << " enter nucleus, E out/in: " << kt->Get4Momentum().e() << " / " << newE <<G4endl;
//      G4cout << " the Field "<< currentField->GetField(kt->GetPosition()) << " "<< kt->GetPosition()<<G4endl;
//      G4cout << " the particle "<<kt->GetDefinition()->GetParticleName()<<G4endl;
	if(newE <= kt->GetActualMass())  // the particle cannot enter the nucleus
	{
// FixMe: should be "pushed back?"
//      for the moment take it past teh nucleus, so we'll not worry next time..
	  FreeTransport(kt, 1.1*t_leave);   // take past nucleus
	  continue;
	}
//
	G4double newP = sqrt(newE*newE- sqr(kt->GetActualMass()));
	G4LorentzVector new4Mom(newP*kt->Get4Momentum().vect().unit(), newE);
	kt->Set4Momentum(new4Mom);
//     G4cout <<" Enter Nucleus - E/Field/Sum: " <<kt->Get4Momentum().e() << " / "
//    	   << (*theFieldMap)[encoding]->GetField(kt->GetPosition()) << " / "
//	   << kt->Get4Momentum().e()-currentField->GetField(kt->GetPosition())
//	   << G4endl
//	   << " Barrier / field just inside nucleus (0.9999*kt->GetPosition())"
//	   << (*theFieldMap)[encoding]->GetBarrier() << " / "
//	   << (*theFieldMap)[encoding]->GetField(0.9999*kt->GetPosition())
//	   << G4endl;
       }
    }

// FixMe: should I add a control on theCutOnP here?
// Transport the particle into the nucleus
//       G4cerr << "RKPropagation t_leave, curTimeStep " <<t_leave << " " <<currTimeStep<<G4endl;
    G4bool is_exiting=false;
    if(currTimeStep > t_leave)  // particle will exit from the nucleus
    {
    	currTimeStep = t_leave;
	is_exiting=true;
    }

//        G4cerr << "RKPropagation t_leave, curTimeStep " <<t_leave << " " <<currTimeStep<<G4endl;
#ifdef debug_1_RKPropagation
        G4cout << "RKPropagation Ekin, field, p "
	<< kt->Get4Momentum().e() - kt->Get4Momentum().mag() << " "
	<< kt->GetPosition()<<" "
 	<< currentField->GetField(kt->GetPosition())<< G4endl
 	<< kt->Get4Momentum()
 	<< G4endl;
#endif

    G4LorentzVector momold=kt->Get4Momentum();
    G4ThreeVector posold=kt->GetPosition();

    if (! FieldTransport(kt, currTimeStep)) {
        FreeTransport(kt,currTimeStep);
    }
//        G4cout << "RKPropagation Ekin, field, p "
// 	<< kt->Get4Momentum().e() - kt->Get4Momentum().mag() << " "
// 	<< currentField->GetField(kt->GetPosition())<< G4endl
// 	<< kt->Get4Momentum()
// 	<< G4endl;
/*	<< "delta p " << momold-kt->Get4Momentum() << G4endl
	<< "del pos " << posold-kt->GetPosition()
	<< G4endl;
*/

// complete the transport
// FixMe: in some cases there could be a significant
//        part to do still in the nucleus, or we stepped to far... depending on
//        slope of potential
    if(is_exiting)  // particle is exiting
    {
// transport free to a position that is surely out of the nucleus, to avoid
// a new transportation and a new adding the barrier next loop.
      G4double t_in, t_out;
      if(GetSphereIntersectionTimes(kt, t_in, t_out))
      {
        G4double velocity=kt->Get4Momentum().vect().mag()/kt->Get4Momentum().e()*c_light;
	G4double t_min=0.1*fermi/velocity;
	t_out=G4std::max(abs(t_out),t_min); // avoid transport by 0 step not taking it out..
	if(t_in < 0 && t_out >= 0)   //still inside, transport safely out.
	{
	  G4ThreeVector savePos = kt->GetPosition();
	  FreeTransport(kt, 1.1*t_out);
	  // and evaluate the right the energy
	  G4double newE=kt->Get4Momentum().e();

// 	G4cout << " V pos/savePos << "
// 		<< (*theFieldMap)[encoding]->GetField(kt->GetPosition())<< " / "
// 		<< (*theFieldMap)[encoding]->GetField(savePos)
// 		<< G4endl;

	  if ( abs(currentField->GetField(savePos)) > 0. &&
	       abs(currentField->GetField(kt->GetPosition())) > 0.)
	  { // FixMe GF: savePos/pos may be out of nucleus, where GetField(..)=0
	    //           This wrongly adds or subtracts the Barrier here while
	    //           this is done later.
	     newE += currentField->GetField(savePos)
	            - currentField->GetField(kt->GetPosition());
	   }

//       G4cout << " go border nucleus, E in/border: " << kt->Get4Momentum() << " / " << newE <<G4endl;

	   if(newE < kt->GetActualMass())
	   {
//	     G4cout << "RKPropagation-Transport: problem with particle exiting - ignored" << G4endl;
	     continue; // the particle cannot exit the nucleus
	   }
//	   G4cout << "%%%% before update %%%% "<< kt->Get4Momentum()<<G4endl;
	   kt->Update4Momentum(newE);
//	   G4cout << "%%%% beyond update %%%% "<< kt->Get4Momentum()<<G4endl;
//	   G4cout << "Field values: "<<currentField->GetField(savePos)<<" "
//	          <<currentField->GetField(kt->GetPosition())<<" "<<kt->GetDefinition()->GetParticleName()<<G4endl;
	}
      } else
      {
      	  G4cerr << "KineticModel-G4RKPropagation: Positioning problem(ignored)"<< G4endl;
      }

      // add the potential barrier
      // FixMe the Coulomb field is not parallel to mom, this is simple approximation
      G4double newE = kt->Get4Momentum().e()+currentField->GetField(kt->GetPosition());
//      G4cout << " leave nucleus, E in/out: " << kt->Get4Momentum() << " / " << newE <<G4endl;
      if(newE < kt->GetActualMass())
      {  // the particle cannot exit the nucleus
//        G4cout << "HadronKineticModel:RKPropagation: ignoring problem with particle E on exit of nucleus" << G4endl;
	continue;
      }
      kt->Update4Momentum(newE);
    }


  }

}

//----------------------------------------------------------------------------
G4bool G4RKPropagation::FieldTransport(G4KineticTrack * kt, const G4double timeStep)
//----------------------------------------------------------------------------
{
//    G4cout <<"Stepper input"<<kt->Get4Momentum()<<G4endl;
// create the integrator stepper
    //    G4Mag_EqRhs * equation = mapIter->second;
    G4Mag_EqRhs * equation = (*theEquationMap)[kt->GetDefinition()->GetPDGEncoding()];
    G4MagIntegratorStepper * stepper = new G4ClassicalRK4(equation);

// create the integrator driver
    G4double hMin = 1.0e-25*second;   // arbitrary choice. Means 0.03 fm at c
    G4MagInt_Driver * driver = new G4MagInt_Driver(hMin, stepper);

// Temporary: use driver->AccurateAdvance()
  // create the G4FieldTrack needed by AccurateAdvance
    G4double curveLength = 0;
    G4FieldTrack track(kt->GetPosition(),
                       kt->Get4Momentum().vect().unit(), // momentum direction
                       curveLength, // curvelength
		       kt->Get4Momentum().e()-kt->GetActualMass(), // kinetic energy
		       kt->GetActualMass(), // restmass
		       kt->Get4Momentum().beta()*c_light); // velocity
  // integrate
    G4double eps = 0.01;
//    G4cout << "currTimeStep = " << currTimeStep << G4endl;
    if(!driver->AccurateAdvance(track, timeStep, eps))
    {  // cannot track this particle
      G4std::cerr << "G4RKPropagation::FieldTransport() warning: integration error."
         << G4endl << "position " << kt->GetPosition() << " 4mom " <<kt->Get4Momentum()
	 <<G4endl << " timestep " <<timeStep 
		  << G4endl;
      delete driver;
      delete stepper;
      return false;
    }
/*
 *      G4cout <<" E/Field/Sum be4 : " <<mom.e() << " / "
 *     	   << (*theFieldMap)[encoding]->GetField(pos) << " / "
 * 	   << mom.e()+(*theFieldMap)[encoding]->GetField(pos)
 * 	   << G4endl;
 */

 // update the kt
    kt->SetPosition(track.GetPosition());
    G4LorentzVector mom(track.GetMomentum(),sqrt(track.GetMomentum().mag2() + sqr(kt->GetActualMass())));
    kt->Set4Momentum(mom);
//    G4cout <<"Stepper output"<<kt<<" "<<kt->Get4Momentum()<<" "<<kt->GetPosition()<<G4endl;
/*
 *      G4cout <<" E/Field/Sum aft : " <<mom.e() << " / "
 *            << " / " << (*theFieldMap)[encoding]->GetField(pos)<< " / "
 * 	   << mom.e()+(*theFieldMap)[encoding]->GetField(pos)
 * 	   << G4endl;
 */

    delete driver;
    delete stepper;
    return true;
}

//----------------------------------------------------------------------------
G4bool G4RKPropagation::FreeTransport(G4KineticTrack * kt, const G4double timeStep)
//----------------------------------------------------------------------------
{
	G4ThreeVector newpos = kt->GetPosition() +
		               timeStep*c_light/kt->Get4Momentum().e() * kt->Get4Momentum().vect();
	kt->SetPosition(newpos);
	return true;
}

/*
G4bool G4RKPropagation::WillBeCaptured(const G4KineticTrack * kt)
{
  G4double radius = theOuterRadius;

// evaluate the final energy. Il will be captured if newE or newP < 0
  G4ParticleDefinition * definition = kt->GetDefinition();
  G4double mass = definition->GetPDGMass();
  G4ThreeVector pos = kt->GetPosition();
  G4LorentzVector mom = kt->Get4Momentum();
  G4VNuclearField * field = (*theFieldMap)[definition->GetPDGEncoding()];
  G4ThreeVector newPos(0, 0, radius); // to get the field on the surface

  G4double newE = mom.e()+field->GetField(pos)-field->GetField(newPos);

  return ((newE < mass) ? false : true);
}
*/



//----------------------------------------------------------------------------
G4bool G4RKPropagation::GetSphereIntersectionTimes(const G4double radius,
//----------------------------------------------------------------------------
				  const G4ThreeVector & currentPos,
				  const G4LorentzVector & momentum,
				  G4double & t1, G4double & t2)
{
  G4ThreeVector speed = momentum.vect()/momentum.e(); // boost vector
  G4double scalarProd = currentPos.dot(speed);
  G4double speedMag = speed.mag();
  G4double sqrtArg = scalarProd*scalarProd -
    speedMag*speedMag*(currentPos.mag2()-radius*radius);
  if(sqrtArg <= 0.) // particle will not intersect the sphere
  {
//     G4cout << " GetSphereIntersectionTimes sqrtArg negative: " << sqrtArg << G4endl;
     return false;
  }
  t1 = (-scalarProd - sqrt(sqrtArg))/speedMag/speedMag/c_light;
  t2 = (-scalarProd + sqrt(sqrtArg))/speedMag/speedMag/c_light;
  return true;
}

//----------------------------------------------------------------------------
G4bool G4RKPropagation::GetSphereIntersectionTimes(const G4KineticTrack * kt,
				  G4double & t1, G4double & t2)
{
  G4double radius = theOuterRadius + 3*fermi; // "safety" of 3 fermi
  G4ThreeVector speed = kt->Get4Momentum().vect()/kt->Get4Momentum().e(); // bost vector
  G4double scalarProd = kt->GetPosition().dot(speed);
  G4double speedMag = speed.mag();
  G4double sqrtArg = scalarProd*scalarProd -
    speedMag*speedMag*(kt->GetPosition().mag2()-radius*radius);
  if(sqrtArg <= 0.) // particle will not intersect the sphere
  {
//     G4cout << " GetSphereIntersectionTimes sqrtArg negative:  " << sqrtArg << G4endl;
     return false;
  }
  t1 = (-scalarProd - sqrt(sqrtArg))/speedMag/speedMag/c_light;
  t2 = (-scalarProd + sqrt(sqrtArg))/speedMag/speedMag/c_light;
  return true;
}

// Implementation methods

//----------------------------------------------------------------------------
void G4RKPropagation::delete_FieldsAndMap(
//----------------------------------------------------------------------------
	G4std::map <G4int, G4VNuclearField *, G4std::less<G4int> > * aMap)
{
  if(aMap)
  {
    G4std::map <G4int, G4VNuclearField *, G4std::less<G4int> >::iterator cur;
    for(cur = aMap->begin(); cur != aMap->end(); ++cur)
      delete (*cur).second;

    aMap->clear();
    delete aMap;
  }

}

//----------------------------------------------------------------------------
void G4RKPropagation::delete_EquationsAndMap(
//----------------------------------------------------------------------------
	G4std::map <G4int, G4Mag_EqRhs *, G4std::less<G4int> > * aMap)
{
  if(aMap)
  {
    G4std::map <G4int, G4Mag_EqRhs *, G4std::less<G4int> >::iterator cur;
    for(cur = aMap->begin(); cur != aMap->end(); ++cur)
      delete (*cur).second;

    aMap->clear();
    delete aMap;
  }
}
