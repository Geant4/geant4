// GGE (Geant4 Geometry Editor)
// Material Combination  Items
// Toshiaki Kodama

class MaterialCombiItem extends MaterialItem {
  ElementItem elements[];
  float ratios[];
  boolean byNum;
  MaterialCombiItem(String name, float density, int densityUnit,int state, float tempe,int tempeUnit, float pressure, ElementItem[] elements, boolean byNum, float[] ratios){
    this.name = name;
    this.density = density;
    this.densityUnit = densityUnit;
    this.state = state;
//    this.tempeOmit = tempeOmit;
    this.tempe = tempe;
    this.tempeUnit = tempeUnit;
    this.pressure = pressure;
    this.elements = elements;
    this.byNum = byNum;
    this.ratios = ratios;
  }
  MaterialCombiItem(ElementItem[] elements){
    ratios = new float[elements.length];
    this.elements = elements;
  }
  public String getRatioNum(int pos){
    if (byNum){
      return ""+(int)ratios[pos];
    }
    return ""+ratios[pos];
  }

}
