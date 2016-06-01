// GGE (Geant4 Geometry Editor)
// 7.9
import java.awt.*;
import java.io.Serializable;
import java.util.*;

class PhysicalItem extends GGEItem {
//Object implements Serializable 
  static final int X=0, Y=1, Z=2;
  static final String vectorComboName[] = {"X", "Y", "Z"};
  static final int T=0, F=1;
  static final String manyComboName[] = {"true","false"};
  static final int MM=0, CM=1, M=2;
  static final String widthUnitName[] = {"mm", "cm", "m"};
  double width;
  int narabi = Y, widthUnit = MM, many = F;
  int  pInt;
  String name, logVolume, mother, x0, y0, z0, kosu;
  Vector vt;

  PhysicalItem(){
    vt = new Vector(11);
  } 

  public String toString(){
    return name;
  }
  
}






