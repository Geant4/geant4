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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLNuclearDensity.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  NuclearDensity::NuclearDensity(G4int A, G4int Z, IFunction1D *densityF)
    :theA(A), theZ(Z)
  {
    densityFunction = densityF;
    theRadiusParameter = densityFunction->getRadiusParameter();
    theMaximumRadius = densityFunction->getMaximumRadius();
    theDiffusenessParameter = densityFunction->getDiffusenessParameter();
    computeCentralRadius();

    initializeDensity();
    initializeFirstDerivative();
    initializeTransmissionRadii();
  }

  NuclearDensity::NuclearDensity(G4int A, G4int Z, IFunction1D *densityF,
				     G4double radius,
				     G4double maximumRadius,
				     G4double diffuseness)
    :theA(A), theZ(Z),
     densityFunction(densityF),
     theRadiusParameter(radius), theMaximumRadius(maximumRadius),
     theDiffusenessParameter(diffuseness)
  {
    computeCentralRadius();
    initializeDensity();
    initializeFirstDerivative();
    initializeTransmissionRadii();
  }

  NuclearDensity::~NuclearDensity() {
    delete densityFunction;
  }

  void NuclearDensity::initMaterial(G4int iamat, G4int izmat)
  {
    G4double res_dws = 0.0;
    G4double fnor = 0.0;

    G4double rcour = 0.0;
    G4int nbr = 0;

    G4double f_r = 0.0;

    G4double drws = 0.0;

    // parametres moyens de densite de la cible (fermi 2 parametres)
    densityFunction->setRadiusParameter(ParticleTable::getNuclearRadius(iamat,izmat));
    densityFunction->setDiffusenessParameter(ParticleTable::getSurfaceDiffuseness(iamat,izmat));
    densityFunction->setMaximumRadius(ParticleTable::getMaximumNuclearRadius(iamat,izmat));

    drws = densityFunction->getMaximumRadius()/29.0;

    // preparation de la distribution w.s.:
    G4double step = 0.2;
    if (iamat >= 19) {
      step = 0.2;
      res_dws = G4integrate(0.0, 13.5, step);
    }
    else { 
      // preparation de la distribution m.h.o.:
      if(iamat >= 6) {
	step=0.1;
	res_dws = G4integrate(0.0, 10.0, step);
      }
      else {
	// preparation de la distribution gaussienne:
	//	 G4double cte = std::pow(ws->adif,3)*std::sqrt(2.*3.141592654);
	res_dws = 3.0*(std::pow(densityFunction->getDiffusenessParameter(), 3)
		       *std::sqrt(Math::twoPi))/2.0;
      }
    }
    fnor = res_dws;

    // calcul de q/pf=f(r)      
    nbr = G4int(std::floor((densityFunction->getMaximumRadius())/drws + 1.5));
    rcour = -1*(drws);

    G4int j = 0;
    for(G4int i = 0; i < nbr; i++) { // do i=1,nbr
      rcour = rcour + drws;
      if(i == 0) { // 1->0
	j++;
	f_r = 0.0;
	x.push_back(f_r);
	y.push_back(0.0);
	r_t.push_back(0.0);
	tmin.push_back(f_r);
	res_dws = 0.0;
      } else {
	step = rcour/20.;
	if(step >= 0.05) {
	  step = 0.05;
	}
	res_dws = G4integrate(0.0, rcour, step);
	f_r = res_dws/fnor;

	if(f_r >= 0.0)  { // Safeguard against negative f_r
	  f_r = std::pow(f_r,(1./3.));
	  j++;
	  x.push_back(f_r);
	  y.push_back(rcour);
	  r_t.push_back(rcour);
	  tmin.push_back(f_r);
	} else {
	  //	  x.push_back(0.0);
	  //	  y.push_back(rcour);
	  if(std::abs(f_r) > 0.01) {
	    ERROR("i = " << i <<  " f_r " << f_r
		      << " poG4int has been skipped." << std::endl);
	  }
	}
      }
    }
    //	      numberOfPoG4ints = j;
    x[x.size() - 1] = 1.0; // Set the last value to 1.0 (y = rmax)
  }

  G4double NuclearDensity::G4integrate(G4double ami, G4double ama, G4double step) const {
    G4double res = 0.0;
    G4double x1[5];
    for(G4int init_i = 0; init_i < 5; init_i++) {
      x1[init_i] = 0.0;
    }
    G4double ri = ami;
    G4double ra = ama;
    G4int nb = 0;
    G4double acont = 1.0;
    G4double dr = step;

    if(ama <= ami) {
      acont = -1.0;
      ri = ama;
      ra = ami;
    }
  
    x1[0] = 95.0/288.0;
    x1[1] = 317.0/240.0;
    x1[2] = 23.0/30.0;
    x1[3] = 793.0/720.0;
    x1[4] = 157.0/160.0;
    nb = G4int(std::floor(((ra - ri)/dr + 1.0000000001))); // 1.0000000001 -> 0.0
    dr = (ra - ri)/(G4double(nb - 1)); 
    res = 0.0;

    if(nb < 10) {
      ERROR("Not enough G4integration poG4ints" << std::endl);
      return 0.0;
    }
  
    for(G4int i = 0; i < 5; i++) {
      res = res + (densityFunction->getValue(ri)
		   + densityFunction->getValue(ra)) * x1[i];
      ri = ri + dr;
      ra = ra - dr;
    }

    nb = nb - 10;

    if(nb == 0) {
      return (res*dr*acont);
    }

    for(G4int i = 0; i < nb; i++) {
      res = res + densityFunction->getValue(ri);
      ri = ri + dr;
    }

    return(res*dr*acont);
  }

  void NuclearDensity::initializeDensity() {
    initMaterial(theA, theZ);
  }

  void NuclearDensity::initializeFirstDerivative() {
    if(x.empty()) {
      ERROR("INCL::NuclearDensity: No datapoG4ints in the nuclear density"
        << std::endl);
      return;
    }

    for(unsigned int i = 0; i != (x.size() - 1); ++i) { // For nuclear density r(p)
      if((x[i+1] - x[i]) == 0.0) { // Safeguard against division by zero
	s.push_back(0.0);
	continue;
      }
      s.push_back((y[i+1] - y[i])/(x[i+1] - x[i]));
    }
    s.push_back(s[x.size() - 2]);

    for(unsigned int i = 0; i != (x.size() - 1); ++i) { // For local energy
      if((r_t[i+1] - r_t[i]) == 0.0) { // Safeguard against division by zero
	s_loce.push_back(0.0);
	continue;
      }
      s_loce.push_back((tmin[i+1] - tmin[i])/(r_t[i+1] - r_t[i]));
    }
    s_loce.push_back(s_loce[r_t.size() - 2]);
  }

  void NuclearDensity::initializeTransmissionRadii() {
    const G4double r0 = 1.12;
    const G4double theNucleonTransmissionRadius = r0*Math::pow13((G4double)theA) + 0.88;

    transmissionRadius[Proton] = theNucleonTransmissionRadius;
    transmissionRadius[Neutron] = theNucleonTransmissionRadius;
    transmissionRadius[PiPlus] = theCentralRadius;
    transmissionRadius[PiZero] = theCentralRadius;
    transmissionRadius[PiMinus] = theCentralRadius;
    transmissionRadius[DeltaPlusPlus] = theNucleonTransmissionRadius;
    transmissionRadius[DeltaPlus] = theNucleonTransmissionRadius;
    transmissionRadius[DeltaZero] = theNucleonTransmissionRadius;
    transmissionRadius[DeltaMinus] = theNucleonTransmissionRadius;
    transmissionRadius[Composite] = theCentralRadius;
  }

  G4double NuclearDensity::getMaxRFromPLegacy(G4double p) const {
    G4double flin = 0.0;
    G4double xv = p;
    G4double dgx = 0.0;
    G4int j = 0;
    G4double tz= xv - x[0];
    if(tz < 0.0) {
      goto flin1;
    } else if(tz == 0.0) {
      goto flin2;
    } else {
      goto flin3;
    }
  flin1:
    flin = y[0]+s[0]*tz;
    return flin;
  flin2:
    flin = y[0];
    return flin;
  flin3:

    for(unsigned int i=1; i < x.size(); ++i) {
      j=i;
      tz=xv-x[i];
      if(tz < 0.0) {
	goto flin8;
      }	else if(tz == 0.0) {
	goto flin9;
      } else {
	   //            goto 10
	continue;
      }
    }

    goto flin8;

  flin9:
    flin=y[j];
    return flin;

  flin8:
    j=j-1;
    dgx=xv-x[j];
    flin=y[j]+s[j]*dgx;
    return flin;
  }

  G4double NuclearDensity::getMaxRFromPNew(G4double p) const {
    G4double tz = p - x[0];
    G4int j = 0;

    if(tz < 0) {
      return (y[0] + s[0]*tz);
    } else if(tz == 0) {
      return y[0];
    } else { // tz > 0
      for(unsigned int i = 1; i < x.size(); ++i) {
	j = i;
	tz = p - x[j];
	if(tz <= 0) {
	  break;
	}
      }
      if(tz >= 0) {
	return y[j];
      } else if(tz < 0.0) {
	j = j - 1;
	G4double dgx = p - x[j];
	return(y[j] + s[j]*dgx);
      }
    }
    return 0.0;
  }

  G4double NuclearDensity::getMaxRFromP(G4double p) const {
    //return getMaxRFromPLegacy(p);
    return getMaxRFromPNew(p);
  }

  G4double NuclearDensity::getMaxTFromR(G4double r) const {
    G4double tz = r - r_t[0];
    G4int j = 0;

    if(tz < 0) {
      return (tmin[0] + s_loce[0]*tz);
    } else if(tz == 0) {
      return tmin[0];
    } else { // tz > 0
      for(unsigned int i = 1; i < r_t.size(); ++i) {
	j = i;
	tz = r - r_t[j];
	if(tz <= 0) {
	  break;
	}
      }
      if(tz >= 0) {
	return tmin[j];
      } else if(tz < 0.0) {
	j = j - 1;
	G4double dgx = r - r_t[j];
	return(tmin[j] + s_loce[j]*dgx);
      }
    }
    return 0.0;
  }

}
