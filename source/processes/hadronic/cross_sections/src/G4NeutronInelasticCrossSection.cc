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
// By JPW, working, but to be cleaned up. @@@

#include "G4NeutronInelasticCrossSection.hh"
#include "globals.hh"

   G4double G4NeutronInelasticCrossSection::
   GetCrossSection(const G4DynamicParticle* aPart, 
                   const G4Element* anEle, G4double )
   {
      G4double atomicNumber = anEle->GetN();
      G4double nOfProtons = anEle->GetZ();
      return GetCrossSection(aPart->GetKineticEnergy(), atomicNumber, nOfProtons);
   }
   
   G4double G4NeutronInelasticCrossSection::
   GetCrossSection(G4double anEnergy, G4double atomicNumber, G4double nOfProtons)
   {
      G4double kineticEnergy = std::log10(anEnergy/MeV);
      G4double nOfNeutrons = atomicNumber-nOfProtons;
      const G4double p1=1.3773;
      const G4double p2=1.+10./atomicNumber-0.0006*atomicNumber;
      const G4double p3=0.6+13./atomicNumber-0.0005*atomicNumber;
      const G4double p4=7.2449-0.018242*atomicNumber;
      const G4double p5=1.64-1.8/atomicNumber-0.0005*atomicNumber;
      const G4double p6=1.+200./atomicNumber+0.02*atomicNumber;
      const G4double p7=(atomicNumber-70.)*(atomicNumber-200.)/11000.;
      
      double part1 = pi*(p1*p1)*std::log(nOfNeutrons);
      double part2 = 1.+ std::pow(atomicNumber, 1./3.) - p2*(1.-1./std::pow(atomicNumber, 1./3.));

      double firstexp = -p4*(kineticEnergy-p5);
      double first=1.+std::exp(firstexp);
      double corr = 1.+p3*(1.-1./first); 

      double secondexp = -p6*(kineticEnergy-p7);
      double second=1.+std::exp(secondexp);
      double corr2 =1./second;

      double xsec = corr*corr2*part1*part2*10.*millibarn;
      return xsec;
    }
