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
// By JPW, working, but to be cleaned up. @@@

#include "G4NeutronInelasticCrossSection.hh"
#include "globals.hh"

   G4double G4NeutronInelasticCrossSection::
   GetCrossSection(const G4DynamicParticle* aPart, 
                   const G4Element* anEle, G4double aTemperature)
   {
      G4double atomicNumber = anEle->GetN();
      G4double nOfProtons = anEle->GetZ();
      return GetCrossSection(aPart->GetKineticEnergy(), atomicNumber, nOfProtons);
   }
   
   G4double G4NeutronInelasticCrossSection::
   GetCrossSection(G4double anEnergy, G4double atomicNumber, G4double nOfProtons)
   {
      G4double kineticEnergy = log10(anEnergy/MeV);
      G4double nOfNeutrons = atomicNumber-nOfProtons;
      const G4double p1=1.3773;
      const G4double p2=1.+10./atomicNumber-0.0006*atomicNumber;
      const G4double p3=0.6+13./atomicNumber-0.0005*atomicNumber;
      const G4double p4=7.2449-0.018242*atomicNumber;
      const G4double p5=1.64-1.8/atomicNumber-0.0005*atomicNumber;
      const G4double p6=1.+200./atomicNumber+0.02*atomicNumber;
      const G4double p7=(atomicNumber-70.)*(atomicNumber-200.)/11000.;
      
      double part1 = pi*(p1*p1)*log(nOfNeutrons);
      double part2 = 1.+ pow(atomicNumber, 1./3.) - p2*(1.-1./pow(atomicNumber, 1./3.));

      double firstexp = -p4*(kineticEnergy-p5);
      double first=1.+exp(firstexp);
      double corr = 1.+p3*(1.-1./first); 

      double secondexp = -p6*(kineticEnergy-p7);
      double second=1.+exp(secondexp);
      double corr2 =1./second;

      double xsec = corr*corr2*part1*part2*10.*millibarn;
      return xsec;
    }
