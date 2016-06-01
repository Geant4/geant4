// GGE (Geant4 Geometry Editor)
// List Items
// Toshiaki Kodama

import java.io.*;

abstract class MaterialItem extends GGEItem {
  static final int G=0, KG=1;
  static final String densityUnitName[] = {"g/cm3","kg/m3"};
  static final int UNDEFIND=0, SOLID=1, LIQUID=2, GAS=3, VACUUM=4;
  static final String stateName[] = {"Undefind","Solid","Liquid","Gas","Vacuum"};
  static final String stateCode[] = {"kStateUndefind","kStateSolid","kStateLiquid","kStateGas","kVacuum"};
  static final int K=0;
  static final String tempUnitName[] = {"kelvin"};
  static final int PASCAL=0, BAR=1, ATOM=2;
  static final String pressName[] = {"pasc","bar","atm"};
  static final String pressCode[] = {"pascal","bar","atmosphere"};
  int state = SOLID, pressUnit = PASCAL, densityUnit = G, tempeUnit = K;
  float density, tempe, pressure;

}
