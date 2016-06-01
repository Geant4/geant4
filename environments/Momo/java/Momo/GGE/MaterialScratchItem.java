// GGE (Geant4 Geometry Editor)
// Material Scratch Items
// Toshiaki Kodama

class MaterialScratchItem extends MaterialItem {
  int atomNum; float massNum;
  MaterialScratchItem(String name, int atomNum, float massNum, float density, int densityUnit, int state, float tempe, int tempUnit, float pressure){
    this.name = name;
    this.atomNum = atomNum;
    this.massNum = massNum;
    this.density = density;
    this.densityUnit = densityUnit;
    this.state = state;
    this.tempe = tempe;
    this.tempeUnit = tempeUnit;
    this.pressure = pressure;
  }
  MaterialScratchItem(){
  }

}
