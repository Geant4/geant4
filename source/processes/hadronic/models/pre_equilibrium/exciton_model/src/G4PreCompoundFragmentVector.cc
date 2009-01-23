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
// $Id: G4PreCompoundFragmentVector.cc,v 1.10 2009-01-23 10:02:11 antoni Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#include "G4PreCompoundFragmentVector.hh"
#include "G4HadronicException.hh"

const G4PreCompoundFragmentVector & 
G4PreCompoundFragmentVector::
operator=(const G4PreCompoundFragmentVector &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4PreCompoundFragmentVector::operator= meant to not be accessable");
    return *this;
}


G4bool G4PreCompoundFragmentVector::
operator==(const G4PreCompoundFragmentVector &) const
{
    return false;
}

G4bool G4PreCompoundFragmentVector::
operator!=(const G4PreCompoundFragmentVector &) const
{
    return true;
}

G4double G4PreCompoundFragmentVector::
CalculateProbabilitiesOriginal(const G4Fragment & aFragment)
{
  TotalEmissionProbabilityOriginal = 0.0;
  pcfvector::iterator aChannel;
  for (aChannel=theChannels->begin(); aChannel != theChannels->end(); 
       aChannel++) 
    {
      // Calculate emission probailities
      // Compute total (integrated over kinetic energy) emission 
      // probability of a fragment and
      // Summing channel emission probabilities
      TotalEmissionProbabilityOriginal += (*aChannel)->CalcEmissionProbability(aFragment);
//JMQ 21/01/09 for printout
//G4cout<<" Landa_c("<<(*aChannel)->GetName()<<" )="<<(*aChannel)->CalcEmissionProbability(aFragment)<<G4endl;
    }
  return TotalEmissionProbabilityOriginal;
}


	//JMQ 15/01/09 new method
G4double G4PreCompoundFragmentVector::
CalculateProbabilities(const G4Fragment & aFragment)
{

//JMQ 15/01/09
  G4double ProtonReferenceProbability=0.;
  G4double Factord=0.8;
  G4double Factort=0.6;
  G4double Factorhe3=0.6;
  G4double Factora=0.8;
//
  TotalEmissionProbability = 0.0;
  pcfvector::iterator aChannel; 
  for (aChannel=theChannels->begin(); aChannel != theChannels->end(); 
       aChannel++) 
    {
      // Calculate emission probailities
      // Compute total (integrated over kinetic energy) emission 
      // probability of a fragment and
      // Summing channel emission probabilities
  //JMQ el 15/01/09
	if ((*aChannel)->GetA()==1 && (*aChannel)->GetZ()==1)
		{
		ProtonReferenceProbability=(*aChannel)->CalcEmissionProbability(aFragment);
		TotalEmissionProbability += ProtonReferenceProbability;
		}
	if ((*aChannel)->GetA()==2 && (*aChannel)->GetZ()==1) 
		{
		  G4double probJMQ = (*aChannel)->GetMaximalKineticEnergy();
                  G4double kk = (*aChannel)->CalcEmissionProbability(aFragment);
		  if (kk> 0.0) TotalEmissionProbability += Factord*ProtonReferenceProbability;
		}
	 if ((*aChannel)->GetA()==3 && (*aChannel)->GetZ()==1) 
		{
		  G4double probJMQ = (*aChannel)->GetMaximalKineticEnergy();
                  G4double kk = (*aChannel)->CalcEmissionProbability(aFragment);
		  if (kk > 0.0) TotalEmissionProbability += Factort*ProtonReferenceProbability;
		}
	 if ((*aChannel)->GetA()==3 && (*aChannel)->GetZ()==2) 
		{
		  G4double probJMQ = (*aChannel)->GetMaximalKineticEnergy();
                  G4double kk = (*aChannel)->CalcEmissionProbability(aFragment);
		  if (kk > 0.0) TotalEmissionProbability += Factorhe3*ProtonReferenceProbability;
		}
         else 	 if ((*aChannel)->GetA()==4 && (*aChannel)->GetZ()==2) 
		{
		  G4double probJMQ = (*aChannel)->GetMaximalKineticEnergy();
                  G4double kk = (*aChannel)->CalcEmissionProbability(aFragment);
		  if (kk > 0.0) TotalEmissionProbability += Factora*ProtonReferenceProbability;
		}
	else
      TotalEmissionProbability += (*aChannel)->CalcEmissionProbability(aFragment);
 //JMQ 13/11/01
// G4cout<<"En el calculo de las probab. d emision"<<G4endl;
// G4double ZZ=(*aChannel)->GetZ();
// G4double AA=(*aChannel)->GetA();
// G4double proba=(*aChannel)->CalcEmissionProbability(aFragment);
// G4cout<<"Z = "<<ZZ<<"   AA = "<<AA<<"  Probabilidad = "<<proba<<G4endl;
// G4cout<<"-----------------------------------------------------------------"<<G4endl;
//

    }
  return TotalEmissionProbability;
}


G4VPreCompoundFragment * G4PreCompoundFragmentVector::
ChooseFragmentOriginal(void)
{
  const G4int NumOfFrags = theChannels->size();
  std::vector<G4double> running;
  running.reserve(NumOfFrags);
  
  pcfvector::iterator i;
  G4double accumulation = 0.0;
  for (i = theChannels->begin(); i != theChannels->end(); ++i) {
    accumulation += (*i)->GetEmissionProbability();

    running.push_back(accumulation);
  }
	
  // Choose an emission channel
  G4double aChannel = G4UniformRand()*TotalEmissionProbabilityOriginal;
  G4int ChosenChannel = -1;
  std::vector<G4double>::iterator ich;
  for (ich = running.begin(); ich != running.end(); ++ich) 
    {
      if (aChannel <= *ich) 
	{
#ifdef G4NO_ISO_VECDIST
          std::vector<G4double>::difference_type n = 0;
          std::distance(running.begin(),ich,n);
          ChosenChannel = n;
#else
	  ChosenChannel = std::distance(running.begin(),ich);
#endif
	  break;
	}
    }
  running.clear();
  if (ChosenChannel < 0) 
    {
 //     G4cerr
	G4cout
	<< "G4PreCompoundFragmentVector::ChooseFragment: I can't determine a channel\n"
	<< "Probabilities: ORIGINAL ";
      for (i = theChannels->begin(); i != theChannels->end(); ++i) 
	{
	  G4cout << (*i)->GetEmissionProbability() << "  ";
	}
      G4cout << '\n';
G4cout << "TotalEmissionProbabilityOriginal = " << TotalEmissionProbabilityOriginal << G4endl;
      return 0;
    }
  else
    {
      for (i = theChannels->begin(); i != theChannels->end(); ++i) 
	{
	  (*i)->IncrementStage();
	}
    }

  return theChannels->operator[](ChosenChannel);
}

	//JMQ 15/01/09 new method
G4VPreCompoundFragment * G4PreCompoundFragmentVector::
ChooseFragment(void)
{
//JMQ 15/01/09
  G4double ProtonReferenceProbability=0.;
  G4double Factord=0.8;
  G4double Factort=0.6;
  G4double Factorhe3=0.6;
  G4double Factora=0.8;
//


  const G4int NumOfFrags = theChannels->size();
  std::vector<G4double> running;
  running.reserve(NumOfFrags);
  
  pcfvector::iterator i;
  G4double accumulation = 0.0;
  for (i = theChannels->begin(); i != theChannels->end(); ++i) {

//JMQ 15/01/09
	if ((*i)->GetA()==1 && (*i)->GetZ()==1)
		{
 		ProtonReferenceProbability=(*i)->GetEmissionProbability();
		accumulation += ProtonReferenceProbability;
		}
	if ((*i)->GetA()==2 && (*i)->GetZ()==1) 
		{
//		  G4double probJMQ = (*i)->GetMaximalKineticEnergy();
//JMQ 22/01/09 if the true total emission probability is cero (i.e. below the barrier => no emission)
		  G4double probJMQ = (*i)->GetEmissionProbability();
		  if (probJMQ > 0.0) accumulation +=Factord*ProtonReferenceProbability;
		}
	 if ((*i)->GetA()==3 && (*i)->GetZ()==1) 
		{
//		  G4double probJMQ = (*i)->GetMaximalKineticEnergy();
		  G4double probJMQ = (*i)->GetEmissionProbability();
		  if (probJMQ > 0.0) accumulation +=Factort*ProtonReferenceProbability;
		}
	 if ((*i)->GetA()==3 && (*i)->GetZ()==2) 
		{
//		  G4double probJMQ = (*i)->GetMaximalKineticEnergy();
		  G4double probJMQ = (*i)->GetEmissionProbability();
		  if (probJMQ > 0.0) accumulation +=Factorhe3*ProtonReferenceProbability;
		}
	else if ((*i)->GetA()==4 && (*i)->GetZ()==2) 
		{
//		  G4double probJMQ = (*i)->GetMaximalKineticEnergy();
		  G4double probJMQ = (*i)->GetEmissionProbability();
//	if ((*i)->GetMaximalKineticEnergy()>0 && (*i)->GetEmissionProbability() ==0)
//{G4cout<<"GetMaximalKineticEnergy="<<(*i)->GetMaximalKineticEnergy()<<"  y GetEmissionProbability="<<(*i)->GetEmissionProbability()<<G4endl;}
		  if (probJMQ > 0.0) accumulation +=Factora*ProtonReferenceProbability;
		}


	else
    accumulation += (*i)->GetEmissionProbability();

 //JMQ 12/01/09
 // G4cout<<"En el sorteo del trozo a emitir"<<G4endl;
 //G4double ZZ=(*i)->GetZ();
 //G4double AA=(*i)->GetA();
 //G4double proba=(*i)->GetEmissionProbability();
 //G4cout<<"Z = "<<ZZ<<"   A = "<<AA<<"  Probabilidad = "<<proba<<G4endl;
 //G4cout<<"-----------------------------------------------------------------"<<G4endl;
//

    running.push_back(accumulation);
  }
	
  // Choose an emission channel
  G4double aChannel = G4UniformRand()*TotalEmissionProbability;
  G4int ChosenChannel = -1;
  std::vector<G4double>::iterator ich;
  for (ich = running.begin(); ich != running.end(); ++ich) 
    {
      if (aChannel <= *ich) 
	{
#ifdef G4NO_ISO_VECDIST
          std::vector<G4double>::difference_type n = 0;
          std::distance(running.begin(),ich,n);
          ChosenChannel = n;
#else
	  ChosenChannel = std::distance(running.begin(),ich);
#endif
	  break;
	}
    }
//MAC  running.clear();
  if (ChosenChannel < 0) 
    {
 //     G4cerr
      G4cout << "G4PreCompoundFragmentVector::ChooseFragment: I can't determine a channel\n"
	<< "Probabilities: ";
//      for (i = theChannels->begin(); i != theChannels->end(); ++i) 
        for (std::vector<G4double>::iterator it = running.begin(); it != running.end(); it++)
	{
//	  G4cout << (*i)->GetEmissionProbability() << "  ";
          G4cout << (*it) << "  ";
	}
//      G4cout << '\n';
      G4cout<<" TotalEmissionProbability ="<<TotalEmissionProbability<<G4endl;
      return 0;
    }
  else
    {
      for (i = theChannels->begin(); i != theChannels->end(); ++i) 
	{
	  (*i)->IncrementStage();
	}
    }
running.clear();
//  return theChannels->operator[](ChosenChannel);
G4VPreCompoundFragment * elegido= theChannels->operator[](ChosenChannel);
//G4cout<<"Canal elegido ="<<elegido->GetName()<<G4endl;
return elegido;
}
