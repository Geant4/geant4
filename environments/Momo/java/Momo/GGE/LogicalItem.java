// GGE (Geant4 Geometry Editor)
// Toshiaki Kodama

import java.awt.*;

class LogicalItem extends GGEItem {
  static final int MM=0, CM=1, M=2, MICROM=3, NANOM=4;
  static final String lengthName[] = {"m","cm","mm","micrometer","nanometer"};
  static final int RAD=0, MRAD=1, DEG=2;
  static final String angleName[] = {"rad","mrad","deg"};

  static final int BOX=0, TRD=1, TUBS=2, CONS=3, PCONE=4, PGON=5;
  static final String solidName[] = {"Box","Trd","Tubs","Cons", 
                                                          "PCone", "PGon"};
  static final String paraName[][] = {{"pX", "pY", "pZ"},
        {"XHalfLength","YHalfLength","ZHalfLength"},
        {"XHalfLengthAt-dz", "XHalfLengthAt+dz", "YHalfLengthAt-dz",
                                        "YHalfLengthAt+dz", "ZHalfLength"},
        {"InnerRadius", "OuterRadius", "ZHalfLength", "StartPhiAngle",
                                                          "DeltaPhiAngle"},
        {"InnerRadiusAt-dz", "InnerRadiusAt+dz", "OuterRadiusAt-dz",
         "OuterRadiusAt+dz","ZHalfLength","StartPhiAngle","DeltaPhiAngle"},
        {"InnerRadius","OuterRadius","StartPhiAngle","DeltaPhiAngle",
                                      "StartThetaAngle","DeltaThetaAngle"},
        {"XHalfLength","YHalfLength","ZHalfLength","Alpha", "Theta","phi"},
        {"InnerRadius","OuterRadius","SweptRadius","StartPhiAngle",
                                                          "DeltaPhiAngle"},
        {"InnerRadius","OuterRadius","InnerStereo","OuterStereo",
                                                            "ZHalfLength"},


        {"StartPhiAngle", "DeltaPhiAngle","Z of a Section", "OuterRadius",
                                                              "InnerRadius"},
        {"StartPhiAngle", "DeltaPhiAngle", "NumberOfSides","Z of a Section",
                                                "OuterRadius","InnerRadius"},
  };
  static final boolean isAngl[][] = {{false, false, false},
                  {false, false, false, false, false},
                  {false, false, false, true, true},
       	          {false, false, false, false, false, true, true},
                  {false, false, true, true, true, true},
                  {false, false, false, true, true,true},
                  {false, false, false, true, true},
                  {false, false, true, true, false},
                  {true, true, false, false, false },
                  {true, true, false, false, false, false}
  };

  MaterialItem mat;
  Color visColour;
  int solidNum;
  float params[];
  int units[];

  LogicalItem(String name, MaterialItem mat, int sol, Color visColour){
    this(sol);
    this.name = name;
    this.mat = mat;
    this.visColour = visColour;
  }
  LogicalItem(int sol){
    this.solidNum = sol;
    params = new float[paraName[solidNum].length];
    units = new int[isAngl[solidNum].length];
    mat = new MaterialScratchItem();
    mat.name = "";
    visColour = Color.red; //Karidayo
  }
  String[] getChoiceName(int pos){
    return isAngl[solidNum][pos] ? angleName : lengthName ;
  }
  int getParaLength(){
    return paraName[solidNum].length;
  }
  String getSolidName(){
    return solidName[solidNum];
  }
  String getParaName(int pos){
    return paraName[solidNum][pos];
  }
}
