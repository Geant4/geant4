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
/// \file ClassifiedDamage.cc
/// \brief Implementation of the ClassifiedDamage class

#include "ClassifiedDamage.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ClassifiedDamage::ClassifiedDamage()
{
	fType = fNA;
	fDamage.clear();
	fBp_begin = 0;
	fBp_end = 0;
	fComplexity = -1;
	fIncludeBase = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::AddDamage(Damage pDmg)
{
	fDamage.push_back(pDmg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::ComputeBp()
{
	fBp_begin = fDamage[0].GetCopyNb();
	fBp_end = fDamage[fDamage.size()-1].GetCopyNb();

	fBp_barycenter=0;
	for(auto it=fDamage.begin();it!=fDamage.end();it++)
	{
		fBp_barycenter+=it->GetCopyNb();
	}
	fBp_barycenter = fBp_barycenter/fDamage.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::ComputePosition()
{
 	// TODO
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::ComputeComplexity()
{
	if(fDamage.size()==0)
	{
		fComplexity = -1;
	}
	else
	{
		fComplexity = -1;
		for(auto it=fDamage.begin();it!=fDamage.end();it++)
		{
			if((it->GetDamageType()==Damage::DamageType::fBackbone))
			{
				fComplexity++;
			}
			else
			{
				if(fIncludeBase)
				{
					fComplexity++;
				}
			}
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::ComputeType()
{
	bool firstStrandTouched = false;
	bool secondStrandTouched = false;

	if(fDamage.size()==0)
	{
		fType = fNA;
	}
	else
	{
		for(auto it=fDamage.begin();it!=fDamage.end();it++)
		{
			if((it->GetDamageType()==Damage::DamageType::fBackbone))
			{
				int strand = it->GetStrand();
				if(strand == 1)
				{
					firstStrandTouched = true;
				}
				if(strand == 2)
				{
					secondStrandTouched = true;
				}
			}

			if (it->GetCause() == Damage::DamageCause::fDirect) {
				fIsThereDirectContribution = true;
			}

			if (it->GetCause() == Damage::DamageCause::fIndirect) {
				fIsThereIndirectContribution = true;
			}
		}

		if(firstStrandTouched && secondStrandTouched)
		{
			fType = fDSB;
		}
		else
		{
			fType = fSSB;
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::Reset()
{
	fType = fNA;
	fDamage.clear();
	fBp_begin = 0;
	fBp_end = 0;
	fComplexity = -1;
}
