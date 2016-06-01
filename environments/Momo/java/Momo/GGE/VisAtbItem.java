import java.awt.*;
import java.io.*;
import java.io.Serializable;

class VisAtbItem extends GGEItem {
  String name, atbsValue[];
  boolean atbsIsValid[];
  static final String atbMethodName[] = {"Visibility", "LineStyle", "LineWidth", "ForceWireFrame", "ForceSolid" };
  Color color;
  VisAtbItem(){
    atbsValue = new String[atbMethodName.length];
    atbsIsValid = new boolean[atbMethodName.length];
  }
  public String toString(){
    return name;
  }

  String getCPP(){
    StringBuffer cpp = new StringBuffer("G4VisAttributes * ");
    cpp.append(name+"= new G4VisAttributes( G4colour(");
    cpp.append(color.getRed()+","+color.getGreen()+","+color.getBlue()+" ));\n");
    for (int i=0; i<atbsIsValid.length; i++){
      if (atbsIsValid[i]){
	cpp.append(name + "->Set"+atbMethodName[i]+"("+atbsValue[i]+");\n");
      }
    }
    return cpp.toString();
  }

}



