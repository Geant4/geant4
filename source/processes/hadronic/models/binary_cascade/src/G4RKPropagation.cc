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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
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

#include "G4LorentzRotation.hh"

// unsigned EncodingHashFun(const G4int& aEncoding);

G4RKPropagation::G4RKPropagation() :
theOuterRadius(0), theNucleus(0),
theFieldMap(0), theEquationMap(0),
theField(0)
{ }


G4RKPropagation::~G4RKPropagation()
{
   // free theFieldMap memory
   if(theFieldMap) delete_FieldsAndMap(theFieldMap);

   // free theEquationMap memory
   if(theEquationMap) delete_EquationsAndMap(theEquationMap);

   if (theField) delete theField;
}

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

   theFieldMap = new std::map <G4int, G4VNuclearField*, std::less<G4int> >;

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

   theEquationMap = new std::map <G4int, G4Mag_EqRhs*, std::less<G4int> >;

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


//#define debug_1_RKPropagation 1
//----------------------------------------------------------------------------
void G4RKPropagation::Transport(G4KineticTrackVector & active,
      //----------------------------------------------------------------------------
      const G4KineticTrackVector &,
      G4double timeStep)
{
   //  reset momentum transfer to field
   theMomentumTranfer=G4ThreeVector(0,0,0);

   // Loop over tracks

   std::vector<G4KineticTrack *>::iterator i;
   for(i = active.begin(); i != active.end(); ++i)
   {
      G4double currTimeStep = timeStep;
      G4KineticTrack * kt = *i;
      G4int encoding = kt->GetDefinition()->GetPDGEncoding();

      std::map <G4int, G4VNuclearField*, std::less<G4int> >::iterator fieldIter= theFieldMap->find(encoding);

      G4VNuclearField* currentField=0;
      if ( fieldIter != theFieldMap->end() ) currentField=fieldIter->second;

      // debug
      //    if ( timeStep > 1e30 ) {
      //	G4cout << " Name :" << kt->GetDefinition()->GetParticleName() << G4endl;
      //    }

      // Get the time of intersections with the nucleus surface.
      G4double t_enter, t_leave;
      // if the particle does not intersecate with the nucleus go to next particle
      if(!GetSphereIntersectionTimes(kt, t_enter, t_leave))
      {
         kt->SetState(G4KineticTrack::miss_nucleus);
         continue;
      }


#ifdef debug_1_RKPropagation
      G4cout <<" kt,timeStep, Intersection times tenter, tleave  "
            <<kt<< " / state= " <<kt->GetState() <<" / " <<" "<< currTimeStep << " / " << t_enter << " / " << t_leave <<G4endl;
#endif

      // if the particle is already outside nucleus go to next  @@GF should never happen? check!
      //  does happen for particles added as late....
      //     if(t_leave < 0 )
      //     {
      //        throw G4HadronicException(__FILE__, __LINE__, "G4RKPropagation:: Attempt to track particle past a  nucleus");
      //        continue;
      //     }

      // Apply a straight line propagation for particle types
      // not included in the model
      if( ! currentField )
      {
         if(currTimeStep == DBL_MAX)currTimeStep = t_leave*1.05;
         FreeTransport(kt, currTimeStep);
         if ( currTimeStep >= t_leave )
         {
            if ( kt->GetState() == G4KineticTrack::inside )
            { kt->SetState(G4KineticTrack::gone_out); }
            else
            { kt->SetState(G4KineticTrack::miss_nucleus);}
         } else if (kt->GetState() == G4KineticTrack::outside && currTimeStep >= t_enter ){
            kt->SetState(G4KineticTrack::inside);
         }

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
            G4double newE = kt->GetTrackingMomentum().e()-currentField->GetField(kt->GetPosition());

            if(newE <= kt->GetActualMass())  // the particle cannot enter the nucleus
            {
               // FixMe: should be "pushed back?"
               //      for the moment take it past the nucleus, so we'll not worry next time..
               FreeTransport(kt, 1.1*t_leave);   // take past nucleus
               kt->SetState(G4KineticTrack::miss_nucleus);
               //	   G4cout << "G4RKPropagation: Warning particle cannot enter Nucleus :" << G4endl;
               //	   G4cout << " enter nucleus, E out/in: " << kt->GetTrackingMomentum().e() << " / " << newE <<G4endl;
               //	   G4cout << " the Field "<< currentField->GetField(kt->GetPosition()) << " "<< kt->GetPosition()<<G4endl;
               // 	   G4cout << " the particle "<<kt->GetDefinition()->GetParticleName()<<G4endl;
               continue;
            }
            //
            G4double newP = std::sqrt(newE*newE- sqr(kt->GetActualMass()));
            G4LorentzVector new4Mom(newP*kt->GetTrackingMomentum().vect().unit(), newE);
            G4ThreeVector transfer(kt->GetTrackingMomentum().vect()-new4Mom.vect());
            G4ThreeVector boost= transfer / std::sqrt(transfer.mag2() + sqr(theNucleus->GetMass()));
            new4Mom*=G4LorentzRotation(boost);
            kt->SetTrackingMomentum(new4Mom);
            kt->SetState(G4KineticTrack::inside);

            /*
     G4cout <<" Enter Nucleus - E/Field/Sum: " <<kt->GetTrackingMomentum().e() << " / "
    	   << (*theFieldMap)[encoding]->GetField(kt->GetPosition()) << " / "
	   << kt->GetTrackingMomentum().e()-currentField->GetField(kt->GetPosition())
	   << G4endl
	   << " Barrier / field just inside nucleus (0.9999*kt->GetPosition())"
	   << (*theFieldMap)[encoding]->GetBarrier() << " / "
	   << (*theFieldMap)[encoding]->GetField(0.9999*kt->GetPosition())
	   << G4endl;
             */
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

#ifdef debug_1_RKPropagation
      G4cerr << "RKPropagation is_exiting?, t_leave, curTimeStep " <<is_exiting<<" "<<t_leave << " " <<currTimeStep<<G4endl;
      G4cout << "RKPropagation Ekin, field, projectile potential, p "
            << kt->GetTrackingMomentum().e() - kt->GetTrackingMomentum().mag() << " "
            << kt->GetPosition()<<" "
            << G4endl << currentField->GetField(kt->GetPosition()) << " "
            <<  kt->GetProjectilePotential()<< G4endl
            << kt->GetTrackingMomentum()
            << G4endl;
#endif

      G4LorentzVector momold=kt->GetTrackingMomentum();
      G4ThreeVector posold=kt->GetPosition();

      //    if (currentField->GetField(kt->GetPosition()) > kt->GetProjectilePotential() ||
      if (currTimeStep > 0 &&
            ! FieldTransport(kt, currTimeStep)) {
         FreeTransport(kt,currTimeStep);
      }

#ifdef debug_1_RKPropagation
      G4cout << "RKPropagation Ekin, field, p "
            << kt->GetTrackingMomentum().e() - kt->GetTrackingMomentum().mag() << " "
            << G4endl << currentField->GetField(kt->GetPosition())<< G4endl
            << kt->GetTrackingMomentum()
            << G4endl
            << "delta p " << momold-kt->GetTrackingMomentum() << G4endl
            << "del pos " << posold-kt->GetPosition()
            << G4endl;
#endif

      // complete the transport
      // FixMe: in some cases there could be a significant
      //        part to do still in the nucleus, or we stepped to far... depending on
      //        slope of potential
      G4double t_in=-1, t_out=0;  // set onto boundary.

      // should go out, or are already out by a too long step..
      if(is_exiting ||
            (GetSphereIntersectionTimes(kt, t_in, t_out) &&t_in<0 && t_out<=0 ))  // particle is exiting
      {
         if(t_in < 0 && t_out >= 0)   //still inside, transport safely out.
         {
            // transport free to a position that is surely out of the nucleus, to avoid
            // a new transportation and a new adding the barrier next loop.
            G4ThreeVector savePos = kt->GetPosition();
            FreeTransport(kt, t_out);
            // and evaluate the right the energy
            G4double newE=kt->GetTrackingMomentum().e();

            // 	G4cout << " V pos/savePos << "
            // 		<< (*theFieldMap)[encoding]->GetField(kt->GetPosition())<< " / "
            // 		<< (*theFieldMap)[encoding]->GetField(savePos)
            // 		<< G4endl;

            if ( std::abs(currentField->GetField(savePos)) > 0. &&
                  std::abs(currentField->GetField(kt->GetPosition())) > 0.)
            { // FixMe GF: savePos/pos may be out of nucleus, where GetField(..)=0
               //           This wrongly adds or subtracts the Barrier here while
               //           this is done later.
               newE += currentField->GetField(savePos)
	                              - currentField->GetField(kt->GetPosition());
            }

            //       G4cout << " go border nucleus, E in/border: " << kt->GetTrackingMomentum() << " / " << newE <<G4endl;

            if(newE < kt->GetActualMass())
            {
#ifdef debug_1_RKPropagation
               G4cout << "RKPropagation-Transport: problem with particle exiting - ignored" << G4endl;
               G4cout << " cannot leave nucleus, E in/out: " << kt->GetTrackingMomentum() << " / " << newE <<G4endl;
#endif
               if (kt->GetDefinition() == G4Proton::Proton() ||
                     kt->GetDefinition() == G4Neutron::Neutron() ) {
                  kt->SetState(G4KineticTrack::captured);
               } else {
                  kt->SetState(G4KineticTrack::gone_out);  //@@GF tofix
               }
               continue; // the particle cannot exit the nucleus
            }
            G4double newP = std::sqrt(newE*newE- sqr(kt->GetActualMass()));
            G4LorentzVector new4Mom(newP*kt->GetTrackingMomentum().vect().unit(), newE);
            G4ThreeVector transfer(kt->GetTrackingMomentum().vect()-new4Mom.vect());
            G4ThreeVector boost= transfer / std::sqrt(transfer.mag2() + sqr(theNucleus->GetMass()));
            new4Mom*=G4LorentzRotation(boost);
            kt->SetTrackingMomentum(new4Mom);
         }
         // add the potential barrier
         // FixMe the Coulomb field is not parallel to mom, this is simple approximation
         G4double newE = kt->GetTrackingMomentum().e()+currentField->GetField(kt->GetPosition());
         if(newE < kt->GetActualMass())
         {  // the particle cannot exit the nucleus  @@@ GF check.
#ifdef debug_1_RKPropagation
            G4cout << " cannot leave nucleus, E in/out: " << kt->GetTrackingMomentum() << " / " << newE <<G4endl;
#endif
            if (kt->GetDefinition() == G4Proton::Proton() ||
                  kt->GetDefinition() == G4Neutron::Neutron() ) {
               kt->SetState(G4KineticTrack::captured);
            } else {
               kt->SetState(G4KineticTrack::gone_out);  //@@GF tofix
            }
            continue;
         }
         G4double newP = std::sqrt(newE*newE- sqr(kt->GetActualMass()));
         G4LorentzVector new4Mom(newP*kt->GetTrackingMomentum().vect().unit(), newE);
         G4ThreeVector transfer(kt->GetTrackingMomentum().vect()-new4Mom.vect());
         G4ThreeVector boost= transfer / std::sqrt(transfer.mag2() + sqr(theNucleus->GetMass()));
         new4Mom*=G4LorentzRotation(boost);
         kt->SetTrackingMomentum(new4Mom);
         kt->SetState(G4KineticTrack::gone_out);
      }

   }

}


//----------------------------------------------------------------------------
G4ThreeVector G4RKPropagation::GetMomentumTransfer() const
//----------------------------------------------------------------------------
{
   return theMomentumTranfer;
}


//----------------------------------------------------------------------------
G4bool G4RKPropagation::FieldTransport(G4KineticTrack * kt, const G4double timeStep)
//----------------------------------------------------------------------------
{
   //    G4cout <<"Stepper input"<<kt->GetTrackingMomentum()<<G4endl;
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
         kt->GetTrackingMomentum().vect().unit(), // momentum direction
         curveLength, // curvelength
         kt->GetTrackingMomentum().e()-kt->GetActualMass(), // kinetic energy
         kt->GetActualMass(), // restmass
         kt->GetTrackingMomentum().beta()*c_light); // velocity
   // integrate
   G4double eps = 0.01;
   //    G4cout << "currTimeStep = " << currTimeStep << G4endl;
   if(!driver->AccurateAdvance(track, timeStep, eps))
   {  // cannot track this particle
#ifdef debug_1_RKPropagation
      std::cerr << "G4RKPropagation::FieldTransport() warning: integration error."
            << G4endl << "position " << kt->GetPosition() << " 4mom " <<kt->GetTrackingMomentum()
            <<G4endl << " timestep " <<timeStep
            << G4endl;
#endif
      delete driver;
      delete stepper;
      return false;
   }
   /*
       G4cout <<" E/Field/Sum be4 : " <<mom.e() << " / "
      	   << (*theFieldMap)[encoding]->GetField(pos) << " / "
  	   << mom.e()+(*theFieldMap)[encoding]->GetField(pos)
  	   << G4endl;
    */

   // Correct for momentum ( thus energy) transfered to nucleus, boost particle into moving nucleus frame.
   G4ThreeVector MomentumTranfer = kt->GetTrackingMomentum().vect() - track.GetMomentum();
   G4ThreeVector boost= MomentumTranfer / std::sqrt (MomentumTranfer.mag2() +sqr(theNucleus->GetMass()));

   // update the kt
   kt->SetPosition(track.GetPosition());
   G4LorentzVector mom(track.GetMomentum(),std::sqrt(track.GetMomentum().mag2() + sqr(kt->GetActualMass())));
   mom *= G4LorentzRotation( boost );
   theMomentumTranfer += ( kt->GetTrackingMomentum() - mom ).vect();
   kt->SetTrackingMomentum(mom);

   //    G4cout <<"Stepper output"<<kt<<" "<<kt->GetTrackingMomentum()<<" "<<kt->GetPosition()<<G4endl;
   /*
    *   G4ThreeVector MomentumTranfer2=kt->GetTrackingMomentum().vect() - mom.vect();
    * G4cout << " MomentumTransfer/corrected" <<    MomentumTranfer << " " <<  MomentumTranfer.mag()
    *  	    <<  " " <<    MomentumTranfer2 << " " <<  MomentumTranfer2.mag() << " "
    *	    << MomentumTranfer-MomentumTranfer2 << " "<<
    *	    MomentumTranfer-MomentumTranfer2.mag() << " " << G4endl;
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
         timeStep*c_light/kt->GetTrackingMomentum().e() * kt->GetTrackingMomentum().vect();
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
  G4LorentzVector mom = kt->GetTrackingMomentum();
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
   G4double speedMag2 = speed.mag2();
   G4double sqrtArg = scalarProd*scalarProd -
         speedMag2*(currentPos.mag2()-radius*radius);
   if(sqrtArg <= 0.) // particle will not intersect the sphere
   {
      //     G4cout << " GetSphereIntersectionTimes sqrtArg negative: " << sqrtArg << G4endl;
      return false;
   }
   t1 = (-scalarProd - std::sqrt(sqrtArg))/speedMag2/c_light;
   t2 = (-scalarProd + std::sqrt(sqrtArg))/speedMag2/c_light;
   return true;
}

//----------------------------------------------------------------------------
G4bool G4RKPropagation::GetSphereIntersectionTimes(const G4KineticTrack * kt,
      G4double & t1, G4double & t2)
{
   G4double radius = theOuterRadius + 3*fermi; // "safety" of 3 fermi
   G4ThreeVector speed = kt->GetTrackingMomentum().vect()/kt->GetTrackingMomentum().e(); // bost vector
   G4double scalarProd = kt->GetPosition().dot(speed);
   G4double speedMag2 = speed.mag2();
   G4double sqrtArg = scalarProd*scalarProd -
         speedMag2*(kt->GetPosition().mag2()-radius*radius);
   if(sqrtArg <= 0.) // particle will not intersect the sphere
   {
      return false;
   }
   t1 = (-scalarProd - std::sqrt(sqrtArg))/speedMag2/c_light;
   t2 = (-scalarProd + std::sqrt(sqrtArg))/speedMag2/c_light;
   return true;
}

// Implementation methods

//----------------------------------------------------------------------------
void G4RKPropagation::delete_FieldsAndMap(
      //----------------------------------------------------------------------------
      std::map <G4int, G4VNuclearField *, std::less<G4int> > * aMap)
{
   if(aMap)
   {
      std::map <G4int, G4VNuclearField *, std::less<G4int> >::iterator cur;
      for(cur = aMap->begin(); cur != aMap->end(); ++cur)
         delete (*cur).second;

      aMap->clear();
      delete aMap;
   }

}

//----------------------------------------------------------------------------
void G4RKPropagation::delete_EquationsAndMap(
      //----------------------------------------------------------------------------
      std::map <G4int, G4Mag_EqRhs *, std::less<G4int> > * aMap)
{
   if(aMap)
   {
      std::map <G4int, G4Mag_EqRhs *, std::less<G4int> >::iterator cur;
      for(cur = aMap->begin(); cur != aMap->end(); ++cur)
         delete (*cur).second;

      aMap->clear();
      delete aMap;
   }
}
