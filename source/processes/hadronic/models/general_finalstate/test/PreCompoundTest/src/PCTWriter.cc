#include "PCTWriter.hh"

#include "PCTCompositeNucleus.hh"
#include "G4Fragment.hh"


void PCTWriter::WriteHeader(const PCTCompositeNucleus& aCN, const G4double aXS)
{
  if ( this->IsOpen() )
    {
      const PCTProjectile * theProjectile = aCN.GetProjectile();
      PCTTarget * theTarget = const_cast<PCTTarget*>(aCN.GetTarget());
      G4int TargetA = 0;
      if (theTarget->IsNucleus()) 
	{
	  TargetA = theTarget->GetA();
	}
      theFile.setf(std::ios::left, std::ios::adjustfield);
      theFile.fill(' ');
      theFile
	<< "# Projectile, Target and Initial Fragment data \n"
	<< "#----------------------------------------------\n"
	<< "% XS = " << std::setprecision(15) << aXS/barn << '\n'
	<< "% PA = " << theProjectile->GetA() << '\n'
	<< "% PZ = " << theProjectile->GetZ() << '\n'
	<< "% PP = (" 
	<< std::setprecision(15) << theProjectile->GetMomentum().x()/MeV << ","
	<< std::setprecision(15) << theProjectile->GetMomentum().y()/MeV << ","
	<< std::setprecision(15) << theProjectile->GetMomentum().z()/MeV << ";"
	<< std::setprecision(15) << theProjectile->GetMomentum().e()/MeV << ")" << '\n'
	<< "% PM = " << std::setprecision(15) << theProjectile->GetMass()/MeV << '\n'
	<< "% TA = " << TargetA << '\n'
	<< "% TZ = " << theTarget->GetZ() << '\n'
	<< "% TM = " <<  std::setprecision(15) << theTarget->GetMass()/MeV << '\n'
	<< "#=====================================================================================\n";
    }
  else 
    {
      std::cout << "PCTWriter::WriteHeader(): Sorry but I can't write anything because there is not open file\n";
    }
  return;
}


void PCTWriter::WriteReaction(const G4int iteration, const G4ReactionProductVector * theDynamicPartVectorPtr,
			      const G4Fragment * theExcitedNucleus)
{
  if ( this->IsOpen() )
    {
      theFile.setf(std::ios::left, std::ios::adjustfield);
      theFile << "* Start Reaction n.: " << iteration << '\n';

      theFile 
	<< std::setw(5) << theExcitedNucleus->GetA()
	<< std::setw(25) << std::setprecision(15) << theExcitedNucleus->GetMomentum().x()/MeV
	<< std::setw(25) << std::setprecision(15) << theExcitedNucleus->GetMomentum().y()/MeV
	<< std::setw(25) << std::setprecision(15) << theExcitedNucleus->GetMomentum().z()/MeV
	<< std::setw(25) << std::setprecision(15) << theExcitedNucleus->GetMomentum().e()/MeV
	<< std::setw(25) << std::setprecision(15) << theExcitedNucleus->GetExcitationEnergy()/MeV
	<< std::setw(25) << std::setprecision(15) << theExcitedNucleus->GetGroundStateMass()/MeV
	<< '\n';

      for (G4ReactionProductVector::const_iterator p = theDynamicPartVectorPtr->begin(); 
	   p != theDynamicPartVectorPtr->end(); ++p) 
	{
	  theFile.setf(std::ios::left, std::ios::adjustfield); 
	  // Particle name
	  theFile << std::setw(20) << (*p)->GetDefinition()->GetParticleName().c_str();
	  theFile.setf(std::ios::right, std::ios::adjustfield); 
	    
	  theFile 
	    // A
	    << std::setw(5) << (*p)->GetDefinition()->GetBaryonNumber() 
	    // Z
	    << std::setw(5) << G4int((*p)->GetDefinition()->GetPDGCharge())
	    // Mass  
	    << std::setw(25) << std::setprecision(15) << (*p)->GetMass()/MeV
	    // Px
	    << std::setw(25) << std::setprecision(15) << (*p)->GetMomentum().x()/MeV
	    // Py
	    << std::setw(25) << std::setprecision(15) << (*p)->GetMomentum().y()/MeV
	    // Pz
	    << std::setw(25) << std::setprecision(15) << (*p)->GetMomentum().z()/MeV
	    // Total energy
	    << std::setw(25) << std::setprecision(15) << (*p)->GetTotalEnergy()/MeV;
	  if ((*p)->GetCreatorModel() != "")
	    // Creator model
	    theFile
	      << std::setw(30) << (*p)->GetCreatorModel().c_str();
	  else 
	    theFile 
	      << std::setw(30) << "Unknown";
	  theFile
	      << '\n';
	}
	theFile << "* End of reaction\n";
    }
  else 
    {
      std::cout << "PCTWriter::WriteReaction(): Sorry but I can't write anything because there is not open file\n";
    }
  return;
}

