#include "G4InuclSpecialFunctions.hh"

#include <iostream.h>

G4double G4InuclSpecialFunctions::bindingEnergyExact(G4double A, 
						     G4double Z) {

  G4int verboseLevel = 2;
  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::bindingEnergyExact" << G4endl;
  }

  // calculates the nuclei binding energy using experimental data or
  // asymptotic mass formula if it's impossible to use exact

  // data for precise binding energies

  const G4double DH[4] = {2.225, 8.482, 5.58, 5.79};

  const G4double DHE[6] = {7.718, 28.297, 27.41, 29.268, 28.826, 31.36};

  const G4double DLI[8] = {4.81, 26.33, 31.995, 39.246, 41.279,
			   45.331, 43.03, 43.13};
   
  const G4double DBE[8] = {26.926, 37.602, 56.502, 58.167, 64.979,
			   65.482, 68.78, 66.08};  

  const G4double DB[9] = {24.65, 34.739, 56.317, 64.752, 76.208,
			  79.578, 84.458, 84.86, 87.76};

  const G4double DC[9] = {39.038, 60.319, 73.444, 92.166,
			  97.112, 105.289, 106.507, 110.759, 114.96};

  const G4double DN[9] = {57.86, 74.039, 94.109, 104.663,
			  115.496, 117.986, 123.87, 126.539, 131.54};
    
  const G4double DOX[9] = {75.567, 98.735, 111.954, 127.624,
			   131.766, 139.813, 143.77, 151.374, 152.57};
      
  const G4double DF[8] = {96.37, 111.412, 128.225, 137.375, 147.806, 
			  154.407, 162.51, 167.707};
      
  const G4double DNE[8] = {112.92, 132.146, 143.785, 160.65, 167.412, 
			   177.778, 182.974, 191.844};
      
  const G4double DNA[9] = {131.78, 145.98, 163.081, 174.152, 186.571, 
			   193.53, 202.541, 208.77,  215.91};

  const G4double DMG[9] = {134.53, 149.204, 168.571, 181.732, 198.262, 205.594,
			   216.688, 223.131, 231.635};

  const G4double DAL[9] = {149.36, 168.71, 183.597, 200.532,
			   211.901, 224.959, 232.684, 242.12, 247.87}; 

  const G4double DSI[9] = {172.0, 187.02, 206.056, 219.366,
			   236.544, 245.018, 255.627, 262.216, 271.431}; 

  const G4double DP[8] = {205.98, 221.424, 239.292,
			  250.617, 262.925, 270.861, 280.966, 287.53}; 

  const G4double DS[10] = {224.71, 243.696, 256.701, 
			   271.789, 280.432, 291.847, 298.835, 308.727, 
			   313.04, 321.068};
     
  const G4double DCL[10] = {244.12, 258.255, 274.066, 285.574, 298.22, 
			    306.801, 317.112, 323.222, 331.296, 337.1}; 

  const G4double DAR[10] = {261.68, 278.748, 291.474, 306.727, 315.515, 
			    327.354, 333.952, 343.822, 349.921, 359.35}; 

  const G4double DK[13] = {278.89, 293.031, 308.584, 320.649, 333.734,
			   341.534, 351.631, 359.165, 368.798, 376.092, 
			   384.97, 391.857, 400.207}; 

  //  const G4double DCA[3] = {296.23, 313.098, 326.43};

  const G4double DSC = 326.957;                               

  G4double DM;

  G4double AN = A - Z;

  G4int IZ = int(Z + 0.1); 

  G4int IN = int(AN + 0.1); 
	 
  switch (IZ) {
  case 1: // H
    if(IN <= 4) {
      DM = DH[IN - 1];
    } else {
      DM = bindingEnergyAsymptotic(A, Z); //::: clean repetitive structure
    };

    break;    
  case 2: // He
    if(IN <= 6) {
      DM = DHE[IN - 1];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 3: // Li
    if(IN <= 8) {
      DM = DLI[IN - 1];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 4: // Be
    if(IN <= 9 && IN >= 2) {
      DM = DBE[IN - 2];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 5: // B
    if(IN <= 10 && IN >= 2) {
      DM = DB[IN - 2];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 6: // C
    if(IN <= 11 && IN >= 3) {
      DM = DC[IN - 3];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 7: // N
    if(IN <= 12 && IN >= 4) {
      DM = DN[IN - 4];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 8: // O
    if(IN <= 13 && IN >= 5) {
      DM = DOX[IN - 5];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 9: // F
    if(IN <= 13 && IN >= 6) {
      DM = DF[IN - 6];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 10: // Ne
    if(IN <= 14 && IN >= 7) {
      DM = DNE[IN - 7];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 11: // Na
    if(IN <= 16 && IN >= 8) {
      DM = DNA[IN - 8];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 12: // Mg
    if(IN <= 16 && IN >= 8) {
      DM = DMG[IN - 8];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 13: // Al
    if(IN <= 17 && IN >= 9) {
      DM = DAL[IN - 9];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 14: // Si
    if(IN <= 18 && IN >= 10) {
      DM = DSI[IN - 10];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 15: // P
    if(IN <= 19 && IN >= 12) {
      DM = DP[IN - 12];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 16: // S
    if(IN <= 22 && IN >= 13) {
      DM = DS[IN - 13];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 17: // Cl
    if(IN <= 23 && IN >= 14) {
      DM = DCL[IN - 14];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 18: // Ar
    if(IN <= 24 && IN >= 15) {
      DM = DAR[IN - 15];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 19: // K
    if(IN <= 28 && IN >= 16) {
      DM = DK[IN - 16];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 20: // Ca
    if(IN >= 17) {
      DM = DK[IN - 17];
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  case 21: // Sc
    if(IN == 19) {
      DM = DSC;
    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    };

    break;    
  default: 

    G4cout << " DM(A, Z): wrong exact case: IN " << IN 
	   << " IZ " << IZ << G4endl;

    DM = bindingEnergyAsymptotic(A, Z);
  }; 

  return DM;
}
