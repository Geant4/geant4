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
/*
 *  \file electromagnetic/TestEm7/src/G4LindhardPartition.cc
 *  \brief Implementation of the G4LindhardPartition class
 *
 *  Created by Marcus Mendenhall on 1/14/08.
 *  2008 Vanderbilt University, Nashville, TN, USA.
 *
 */

//

#include "G4LindhardPartition.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

/*
for a first cut, we will compute NIEL from a Lindhard-Robinson partition
based on the most abundant element in the material.

this is from IEEE Trans. Nucl Science Vol. 48 No.1 February 2001 page 162++
Insoo Jun, "Effects of Secondary Particles on the Total Dose..."
and, by reference,
Lindhard, Nielsen, Scharff & Thompson, 
"Integral Equations Governing Radiation Efects...", 
Mat. Fys. Medd. Dan. Vid. Selsk. vol 33 #10, pp1-42, 1963
and
Robinson, "The dependence of radiation effects on primary recoil energy",
in Proc. Int. Conf. Radiation-Induced Voids in Metal,  
Albany, NY 1972 pp. 397-439
def lindhard_robinson(z1, a1, z2, a2, ke):
el=30.724*z1*z2*math.sqrt(z1**0.6667+z2**0.6667)*(a1+a2)/a2
fl=0.0793*z1**0.6667*math.sqrt(z2)*(a1+a2)**1.5/
((z1**0.6667+z2**0.6667)**0.75*a1**1.5*math.sqrt(a2))
eps=ke*(1.0/el)
return 1.0/(1+fl*(3.4008*eps**0.16667+0.40244*eps**0.75+eps))
*/

G4LindhardRobinsonPartition::G4LindhardRobinsonPartition() 
{
  max_z = 120;
  for(size_t i=1; i<max_z; i++) {z23[i]=std::pow((G4double)i, 2./3.);}
}

G4double G4LindhardRobinsonPartition::PartitionNIEL(
   G4int z1, G4double a1, const G4Material *material, G4double energy) const
{
  size_t nMatElements = material->GetNumberOfElements();
        
  const G4double *atomDensities=material->GetVecNbOfAtomsPerVolume();
  G4double maxdens=0.0;
  size_t maxindex=0;
  for (size_t k=0 ; k < nMatElements ; k++ )
    {
      if(atomDensities[k] > maxdens) {
        maxdens=atomDensities[k];
        maxindex=k;
      }
    }
  const G4Element *element=material->GetElement(maxindex);

  G4int z2=G4int(element->GetZ());
        
  G4double a2=element->GetA()/(Avogadro*amu);
        
  G4double zpow=z23[z1]+z23[z2];
  G4double asum=a1+a2;
        
  G4double el=30.724*z1*z2*std::sqrt(zpow)*asum/a2;
  G4double fl=0.0793*z23[z1]*std::sqrt(z2*asum*asum*asum/(a1*a1*a1*a2))
              /std::pow(zpow, 0.75);
  G4double eps=(energy/eV)*(1.0/el);

  return 
    1.0/(1+fl*(3.4008*std::pow(eps, 0.16667)+0.40244*std::pow(eps, 0.75)+eps));
}

