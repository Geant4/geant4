
#include "ExamplePrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "G4GeneralParticleSource.hh"

ExamplePrimaryGeneratorAction::ExamplePrimaryGeneratorAction()
{
  //DEBUG
  DebugXmin = -999999.;

  particleGun = new G4GeneralParticleSource ();
}

ExamplePrimaryGeneratorAction::~ExamplePrimaryGeneratorAction()
{
  //  if(verbosityLevel == 2)
  //{
      // Always give the user debug info.
  G4cout << "Output of DEBUG stuff" << G4endl;
  G4cout << "positional stuff" << G4endl;
  G4cout << "Scale, X, Scale, Y, Scale, Z" << G4endl;
  G4double scalex, scaley, scalez;
  for(int i=0; i<100; i++)
    {
      scalex = DebugXmin + (i+1)*DebugXStep;
      scaley = DebugYmin + (i+1)*DebugYStep;
      scalez = DebugZmin + (i+1)*DebugZStep;
      G4cout << scalex << "   " << debugx[i] << "     " << scaley << "   " << debugy[i] << "     " << scalez << "   " << debugz[i] << G4endl;
    }
  
  G4cout << "Scale, Number Px, Py, Pz, Scale, Number Theta, Scale, Number Phi" << G4endl;
  G4double scalep, scalet, scalephi;
  for(i=0; i<100; i++)
    {
      scalep = -1 + (i+1)*0.02;
      scalet = (i+1) * (pi/100.);
      scalephi = (i+1) * (twopi/100.);
      G4cout << scalep << "   " << debugpx[i] << "   " << debugpy[i] << "   " << debugpz[i] << "     " << scalet << "   " << debugtheta[i] << "     " << scalephi << "   " << debugphi[i] << G4endl;
    }
  
  G4cout << "Initial Energy,  No. Of Events" << G4endl;
  G4double ene_out = 0.;
  for(i=0; i<100; i++)
    {
      if(EneDisType == "Mono")
	ene_out = emin/2. + (i+1)*debug_energy_step;
      //      else if(EneDisType == "Arb")
      //{
      //  if(IntType == "Spline")
      //    ene_out = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(0)) + (i+1)*debug_energy_step;
      //  else
      //    ene_out = ArbEnergyH.GetLowEdgeEnergy(size_t(0)) + (i+1)*debug_energy_step;
      //}
      //else if(EnergyDisType == "Epn")
      //{
      //  ene_out = IPDFEnergyH.GetLowEdgeEnergy(size_t(0)) + (i+1)*debug_energy_step;
      //}
      else
	ene_out = emin + (i+1)*debug_energy_step;
      G4cout << ene_out << "       " << debugenergy[i] << G4endl;
    }
  //}
  //  G4cout << "About to delete particleGun " << G4endl;
  delete particleGun;
  //G4cout << "Deleted particleGun " << G4endl;
}

void ExamplePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //G4cout << "About to generate primary vertex" << G4endl;
  particleGun->GeneratePrimaryVertex(anEvent);
  //G4cout << "Back " << DebugXmin << G4endl;

  if(DebugXmin == -999999.)
    {
    SourceType = particleGun->GetPosDisType();
    //G4cout << "SourceType " << SourceType << G4endl;
    SourceShape = particleGun->GetPosDisShape();
    //G4cout << "SourceShape " << SourceShape << G4endl;
    radius = particleGun->GetRadius();
    //G4cout << "radius " << radius << G4endl;
    halfx = particleGun->GetHalfX();
    //G4cout << "halfx " << halfx << G4endl;
    halfy = particleGun->GetHalfY();
    //G4cout << "halfy " << halfy << G4endl;
    halfz = particleGun->GetHalfZ();
    //G4cout << "halfz " << halfz << G4endl;
    centre = particleGun->GetCentreCoords();
    //G4cout << "centre " << centre << G4endl;

    // for use with energy
    EneDisType = particleGun->GetEnergyDisType();
    //G4cout << "EneDisType " << EneDisType << G4endl;
    InterpolationType = particleGun->GetIntType();
    //G4cout << "InterpolationType " << InterpolationType << G4endl;
    if(EneDisType == "Arb")
      {
	emin = particleGun->GetArbEmin();
	//G4cout << "emin " << emin << G4endl;
	emax = particleGun->GetArbEmax();
	//G4cout << "emax " << emax << G4endl;
      }
    else
      {
	emin = particleGun->GetEmin();
	//G4cout << "emin " << emin << G4endl;
	emax = particleGun->GetEmax();
	//G4cout << "emax " << emax << G4endl;
      }
    }
  
  G4ThreeVector ParticlePos = particleGun->GetParticlePosition();
  //G4cout << "ParticlePos " << ParticlePos << G4endl;
  G4ThreeVector ParticleMomDir = particleGun->GetParticleMomentumDirection();
  //G4cout << "ParticleMomDir " << ParticleMomDir << G4endl;

  //G4cout << "Starting DEBUG stuff " << DebugXmin << G4endl;
  // DEBUG SECTION
  //  if(verbosityLevel == 2)
  // G4cout << "Collecting DEBUG info" <<G4endl;
  G4int idebug = 0;
  // position
  if(DebugXmin == -999999.)
    {
      //G4cout << "Here 11 " << SourceType << " " << centre <<G4endl;
      if(SourceType == "Point")
	{
	  // DEBUG - make Xmin etc 0.5 * point and Xmax etc 1.5 * point
	  if(centre.x() == 0.0)
	    {
	      DebugXmin = -2.;
	      DebugXmax = 2.;
	    }
	  else
	    {
	      DebugXmin = centre.x() * 0.5;
	      DebugXmax = centre.x() * 1.5;
	    }

	  if(centre.y() == 0.0)
	    {
	      DebugYmin = -2.;
	      DebugYmax = 2.;
	    }
	  else
	    {
	      DebugYmin = centre.y() * 0.5;
	      DebugYmax = centre.y() * 1.5;
	    }

	  if(centre.z() == 0.0)
	    {
	      DebugZmin = -2.;
	      DebugZmax = 2.;
	    }
	  else
	    {
	      DebugZmin = centre.z() * 0.5;
	      DebugZmax = centre.z() * 1.5;
	    }
	}
      else 
	{
	  //G4cout << "Here 11a " << SourceShape << G4endl;
	  if((SourceShape == "Circle") || (SourceShape == "Annulus") || (SourceShape == "Sphere"))
	    {
	      DebugZmax = radius;
	    }
	  else if((SourceShape == "Ellipse") || (SourceShape == "Ellipsoid"))
	    {
	      DebugZmax = halfx;
	      if(halfy > DebugZmax)
		DebugZmax = halfy;
	      if(halfz > DebugZmax)
		DebugZmax = halfz;
	    }
	  else if(SourceShape == "Square")
	    {
	      DebugZmax = halfx;
	    }
	  else if(SourceShape == "Rectangle")
	    {
	      DebugZmax = halfx;
	      if(DebugZmax < halfy)
		DebugZmax = halfy;
	    }
	  else if(SourceShape == "Cylinder")
	    {
	      if(radius >= halfz)
		DebugZmax = radius;
	      else
		DebugZmax = halfz;
	    }
	  else if(SourceShape == "Para")
	    {
	      DebugZmax = halfx;
	      if(DebugZmax < halfy)
		DebugZmax = halfy;
	      if(DebugZmax < halfz)
		DebugZmax = halfz;
	    }
	  DebugZmax = 3 * DebugZmax;
	  DebugXmin = centre.x() - DebugZmax;
	  DebugYmin = centre.y() - DebugZmax;
	  DebugZmin = centre.z() - DebugZmax;
	  DebugXmax = centre.x() + DebugZmax;
	  DebugYmax = centre.y() + DebugZmax;
	  DebugZmax = centre.z() + DebugZmax;
	}
      DebugXStep = (DebugXmax - DebugXmin)/100.;
      DebugYStep = (DebugYmax - DebugYmin)/100.;
      DebugZStep = (DebugZmax - DebugZmin)/100.;
    }

  //G4cout << "out of first bit" << G4endl;
  idebug = 0;
  G4double X_edge = DebugXmin;
  //G4cout << "entering while loop 1 " << ParticlePos.x() << G4endl;
  while (ParticlePos.x() > X_edge)
    {
      //G4cout << idebug << " " << ParticlePos.x() << " " << X_edge << G4endl;
      X_edge = DebugXmin + (idebug+1)*DebugXStep;
      idebug++;
    }
  debugx[idebug-1] = debugx[idebug-1] + 1;
  idebug = 0;
  G4double Y_edge = DebugYmin;
  //G4cout << "entering while loop 2" << G4endl;
  while (ParticlePos.y() > Y_edge)
    {
      Y_edge = DebugYmin + (idebug+1)*DebugYStep;
      idebug++;
    }
  debugy[idebug-1] = debugy[idebug-1] + 1;
  idebug = 0;
  G4double Z_edge = DebugZmin;
  //G4cout << "entering while loop 3" << G4endl;
  while (ParticlePos.z() > Z_edge)
    {
      Z_edge = DebugZmin + (idebug+1)*DebugZStep;
      idebug++;
    }
  debugz[idebug-1] = debugz[idebug-1] + 1;

  //G4cout << "Here 12" << G4endl;
  // trajectory
  // px, py, pz are unit vectors so arrays run -1 to 1
  idebug = 0;
  G4double Px_edge = -1.;
  //G4cout << "Loop 1" << G4endl;
  while (ParticleMomDir.x() > Px_edge)
    {
      Px_edge = -1 + (idebug+1)*0.02;
      idebug++;
    }
  debugpx[idebug-1] = debugpx[idebug-1] + 1;
  idebug = 0;
  G4double Py_edge = -1.;
  //G4cout << "Loop 2" << G4endl;
  while (ParticleMomDir.y() > Py_edge)
    {
      Py_edge = -1 + (idebug+1)*0.02;
      idebug++;
    }
  debugpy[idebug-1] = debugpy[idebug-1] + 1;
  idebug = 0;
  G4double Pz_edge = -1.;
  //G4cout << "Loop 3" << G4endl;
  while (ParticleMomDir.z() > Pz_edge)
    {
      Pz_edge = -1 + (idebug+1)*0.02;
      idebug++;
    }
  debugpz[idebug-1] = debugpz[idebug-1] + 1;

  G4double theta = particleGun->GetTheta();
  G4double phi = particleGun->GetPhi();
  
  //G4cout << "Here 13" << G4endl;
  // Theta ranges 0-Pi, and Phi goes 0-two pi.
  idebug = 0;
  G4double Theta_edge = 0;
  while (theta > Theta_edge)
    {
      Theta_edge =  (idebug+1) * (pi/100.);
      idebug++;
    }
  debugtheta[idebug-1] = debugtheta[idebug-1] + 1;
  idebug = 0;
  G4double Phi_edge = 0.;
  while (phi > Phi_edge)
    {
      Phi_edge = (idebug+1) * (twopi/100.);
      idebug++;
    }
  debugphi[idebug-1] = debugphi[idebug-1] + 1;

  //G4cout << "Here 14" << G4endl;
  // Energy
  if(EneDisType == "Mono")
    debug_energy_step = emin/100.;
  //  else if(EneDisType == "Arb")
  // {
  //   if(InterpolationType == "Spline")
  //{
  //  G4int len = G4int(IPDFArbEnergyH.GetVectorLength());
  //  debug_energy_step = (IPDFArbEnergyH.GetLowEdgeEnergy(size_t(len-1)) - IPDFArbEnergyH.GetLowEdgeEnergy(size_t(0)))/100.;
  //}
  //  else
  //{
  //  G4int len = G4int(ArbEnergyH.GetVectorLength());
  //  debug_energy_step = (ArbEnergyH.GetLowEdgeEnergy(size_t(len-1)) - ArbEnergyH.GetLowEdgeEnergy(size_t(0)))/100.;
  //}
  //}
  //else if(EneDisType == "Epn")
  // {
  //  G4int len = G4int(IPDFEnergyH.GetVectorLength());
  //  debug_energy_step = (IPDFEnergyH.GetLowEdgeEnergy(size_t(len-1)) - IPDFEnergyH.GetLowEdgeEnergy(size_t(0)))/100.;
  //}
  else
    debug_energy_step = (emax - emin)/100.;

  //G4cout << "Here 15" << G4endl;
  G4double PartEnergy = particleGun->GetParticleEnergy();
  //G4cout << "Energy is " << PartEnergy << " " << emin << " " << emax << G4endl;

  G4double Ebin_edge = 0.;
  idebug = 0;
  while (PartEnergy > Ebin_edge)
    {
      if(EneDisType == "Mono")
	Ebin_edge = emin/2. + (idebug+1)*debug_energy_step;
      //  else if(EneDisType == "Arb")
      //{
      //  if(InterpolationType == "Spline")
      //    Ebin_edge = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(0)) + (idebug+1)*debug_energy_step;
      //  else
      //    Ebin_edge = ArbEnergyH.GetLowEdgeEnergy(size_t(0)) + (idebug+1)*debug_energy_step;
      //}
      // else if(EneDisType == "Epn")
      //{
      //  Ebin_edge = IPDFEnergyH.GetLowEdgeEnergy(size_t(0)) + (idebug+1)*debug_energy_step;
      //}
      else
	Ebin_edge = emin + (idebug+1)*debug_energy_step;
      idebug++;
    }
  debugenergy[idebug-1] = debugenergy[idebug-1] + 1;

  //  G4cout << "debug-energy_step " << debug_energy_step << " " << Ebin_edge << G4endl;
  //G4cout << "Ending thingy" << G4endl;

}



