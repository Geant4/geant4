// GAG (Geant4 Adaptive GUI)
// Requires : JDK 1.1.3 (or later) + swing-1.0.2
// Version : alpha  1998/Mar/31
// GEANT4 GUI Group, Toshiaki Kodama
// 1998 July 6 GEANT4 Beta-01

package GAG;

class GAGcommandItem extends Object {
  static final int NAME=0, GUIDE=1, TYPE=2, OMITTABLE=3, DEFAULT=4, RANGE=5, CANDIDATE=6;
  static final int MAX = 7;
  private String name, guide[], range;
  private int paramEnts;
  String params[][];
  GAGcommandItem(String name, String guide[], String range, int ents){
    this.name = name;
    this.guide = guide;
    paramEnts = ents;
    if (ents > 0) {
      params = new String[paramEnts][MAX];
    }
    this.range = range;
    return;
  }
  String getCommandName(){
    return name;
  }
  String[] getCommandGuide(){
    return guide;
  }
  String getTitle(){
    return guide[0];
  }
  String getCommandRange(){
    return range;
  }
  int getParamEntries(){
    return paramEnts;
  }
  String getParamName(int pos){
    return params[pos][NAME];
  }
  String getParamGuide(int pos){
    return params[pos][GUIDE];
  }
  String getParamType(int pos){
    return params[pos][TYPE];
  }
  boolean isOmittable(int ent){
    return params[ent][OMITTABLE].equals("1");
  }
  String getParamDefault(int pos){
    return params[pos][DEFAULT];
  }
  String getParamRange(int pos){
    return params[pos][RANGE];
  }
  String getParamCandidate(int pos){
    return params[pos][CANDIDATE];
  }
}
