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
// $Id: G4LowEnergyPolarizedRayleigh.hh,v 1.3 2005-05-12 09:22:05 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//
// File name:     G4LowEnergyPolarizedRayleigh.hh
//
// Author:        Capra Riccardo
//
// Creation date: May 2005
//
// History:
// -----------
// 02 May 2005  R. Capra         1st implementation
//
//----------------------------------------------------------------

//! \file    G4LowEnergyPolarizedRayleigh.hh
//! \brief   Header file of G4LowEnergyPolarizedRayleigh class
//! \author  Capra Riccardo
//! \date    May 2005
//! \par     History:
//! <TABLE>
//!  <TR><TD> 02 May 2005 </TD><TD> R. Capra	</TD><TD> 1<SUP>st</SUP> implementation </TD></TR>
//! </TABLE>
//! \sa      G4LowEnergyPolarizedRayleigh.cc G4VLowEnergyDiscretePhotonProcess.hh         

#ifndef   G4LowEnergyPolarizedRayleigh_hh
 #define  G4LowEnergyPolarizedRayleigh_hh
 
 #include "G4VLowEnergyDiscretePhotonProcess.hh"

 //! \class   G4LowEnergyPolarizedRayleigh
 //! \brief   The class implements the polarized Rayleigh process
 //!
 //!          Specialization of G4VLowEnergyDiscretePhotonProcess class.
 //!          This is a preliminary implementation of the Rayleigh polarized process.
 //!          Reference articles are:
 //! <UL>
 //!  <LI> Rayleigh process:                     &quot;The Quantum Theory of Radiation,&quot;
 //!                                                   W. Heitler,       Oxford at the Clarendon Press, Oxford (1954)                                                 </LI>
 //!  <LI> Scattering function:                  &quot;A simple model of photon transport,&quot;
 //!                                                   D.E. Cullen,      Nucl. Instr. Meth. in Phys. Res. B 101 (1995) 499-510                                        </LI>
 //!  <LI> Polarization of the outcoming photon: &quot;Beam test of a prototype detector array for the PoGO astronomical hard X-ray/soft gamma-ray polarimeter,&quot; 
 //!                                                   T. Mizuno et al., Nucl. Instr. Meth. in Phys. Res. A 540 (2005) 158-168                                        </LI>
 //! </UL>
 class G4LowEnergyPolarizedRayleigh : public G4VLowEnergyDiscretePhotonProcess
 {
  private:
   //! \brief Hides copy constructor
                                                G4LowEnergyPolarizedRayleigh(const G4LowEnergyPolarizedRayleigh &);

   //! \brief Hides assignment operator
   G4LowEnergyPolarizedRayleigh &               operator=(const G4LowEnergyPolarizedRayleigh &);



  public:
   //! \brief Class constructor
   //! \param processName The name of the process
                                                G4LowEnergyPolarizedRayleigh(const G4String &processName = "polarLowEnRayleigh");

   //! \brief Class destructor
   inline virtual                              ~G4LowEnergyPolarizedRayleigh(void) {}



   //! \brief Calculates the new photon direction and polarization
   //!
   //!        In the code the following definitions are used <IMG src="../img/definitions.gif">
   //! \param aTrack The track
   //! \param aStep  The step
   //! \return The proposed changes for the photon
   virtual G4VParticleChange *                  PostStepDoIt(const G4Track & aTrack, const G4Step & aStep);
   
  private:
   //! \brief Generates \f$cos \left ( \theta\right )\f$ of the scattered photon
   //! \param incomingPhotonEnergy The energy of the incoming photon
   //! \param zAtom Atomic number
   //! \return \f$cos \left ( \theta\right )\f$
   G4double                                     GenerateCosTheta(G4double incomingPhotonEnergy, G4int zAtom) const;
   
   //! \brief Generates \f$\phi\f$ of the scattered photon
   //! \param cosTheta \f$cos \left ( \theta\right )\f$ of the scattered photon
   //! \return  \f$\phi\f$
   G4double                                     GeneratePhi(G4double cosTheta) const;

   //! \brief Generates the polarization direction \f$\beta\f$ in the plane x, y relative to the x direction
   //! \return \f$\beta\f$
   G4double                                     GeneratePolarizationAngle(void) const;
   
   
   
   //! \brief Lower energy limit due to the available data in the physics tables
   G4double                                     intrinsicLowEnergyLimit;

   //! \brief Higher energy limit due to the available data in the physics tables
   G4double                                     intrinsicHighEnergyLimit;
 };

#endif /* G4LowEnergyPolarizedRayleigh_hh */
