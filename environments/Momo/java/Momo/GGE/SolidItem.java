// GGE (Geant4 Geometry Editor)
// Toshiaki Kodama

import java.io.Serializable;

abstract class SolidItem extends Object implements Serializable {
  boolean isEmpty;
  protected int solidType;
  SolidItem(int solidType){
    this.solidType = solidType;
  }
  int getSolidType(){
    return solidType;
  }
  abstract String getCPP(String name);
}
