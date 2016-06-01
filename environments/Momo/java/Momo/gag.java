// GAG (Geant4 Adaptive GUI)
// Requires : JDK 1.1.3 (or later) +  swing-1.0.2
// Version : alpha  1998/Mar/30
// GEANT4 GUI Group  Toshiaki Kodama
// 1998 July 6 GEANT4 Beta-01

import GAG.GAG;

public class gag {
  public static void main(String args[]){
    GAG gag = new GAG();
    gag.mainLoop(args);
    System.exit(0);
  }
}
