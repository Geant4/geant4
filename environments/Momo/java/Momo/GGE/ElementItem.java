// GGE (Geant4 Geometry Editor)
// Element Item
// Toshiaki Kodama

class ElementItem extends GGEItem {
//Object {
  String name, symbol;
  int atomNum;  double massNum;
  boolean isMetal;

  ElementItem(String name, String symbol, int atomNum, double massNum, boolean isMetal){
    this.name = name;
    this.symbol = symbol;
    this.atomNum = atomNum;  this.massNum = massNum;
    this.isMetal = isMetal;
  }
  String getCPP(){
    return new String("G4Element* element" +symbol+ " = new G4Element( \"" + name + "\", \"" + symbol + "\", " + atomNum + ". , " + massNum + "*g/mole );\n");
  }
  public String toString(){
    return symbol;
  }
}
