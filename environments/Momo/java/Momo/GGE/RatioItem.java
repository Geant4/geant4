// GGE (Geant4 Geometry Editor)
// Toshiaki Kodama

import java.io.Serializable;

class RatioItem extends Object implements Serializable {
  ElementItem elems[];
  float ratio[];
  int isFraction;
  boolean isEmpty;

  RatioItem(ElementItem eis[]){
    elems = eis;
    ratio = new float[elems.length];
    isEmpty = true;
  }
  int getLength(){
    return elems.length;
  }
  String getEleName(int pos){
    return elems[pos].name;
  }
  public String toString(){
    StringBuffer stb = new StringBuffer();
    for (int i=0; i<elems.length; i++){
      stb.append(elems[i].symbol);
      if (!isEmpty) stb.append(":"+getRatioStr(i)+" "); else stb.append("  ");
    }
    return stb.toString();
  }
  String getRatioStr(int pos){
    if(isFraction==1){
      return ""+ratio[pos];
    }else{
      return ""+Math.round(ratio[pos]);
    }
  }
  public String getCPP(String name){
    StringBuffer cpp = new StringBuffer();
    for (int i=0; i<elems.length; i++){
      cpp.append(name + "->AddElement( element" + elems[i].symbol + ", ");
      cpp.append(getRatioStr(i)+" );\n");
    }
    return cpp.toString();
  }
}
