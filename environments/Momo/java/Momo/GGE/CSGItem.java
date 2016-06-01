// GGE (Geant4 Geometry Editor)
// Tetsuya Yamada

import java.lang.Math;

class CSGItem extends SolidItem {
  static final int BOX=0, TRD=1, TUBS=2, CONS=3 , SPHERE=4, PARALLEL=5,
                   TORUS=6, HYPE=7;
  static final String solidName[] = {"Box","SymTrapezoid",
                      "TubeSegment","ConeSegment","SphereSegment",
                                  "Parallelepiped","TorusSegment","HYPE"};

  static final String csgName[] = {"Box","Trd","Tubs",
                           "Cons","Sphere","Para","Torus","Hype"};

  static final String csgDAWNName[] = {"Box","Trd","Tubs",
                           "Cons","SphereSeg","Parallelepiped","Torus","Hype"};

  static final String solidConstr[] = {"Box","Trapezoid","TubeSegment",
                           "ConeSegment","SphereSegment","Parallelepiped",
                                                    "TorusSegment","HYPE"};
  static final String paraName[][] = {
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
    };
  static final boolean isAngl[][] = {{false, false, false},
                  {false, false, false, false, false},
                  {false, false, false, true, true},
       	          {false, false, false, false, false, true, true},
                  {false, false, true, true, true, true},
                  {false, false, false, true, true,true},
                  {false, false, false, true, true},
                  {false, false, true, true, false},
       };

  private double values[];
  private String units[];
  private double mmValue;
  private double bound;


  CSGItem(int solidType){
    super(solidType);
    values = new double[paraName[solidType].length];
    units = new String[isAngl[solidType].length];
    isEmpty = true;
  }
  int getParaCount(){
    return paraName[solidType].length;
  }
  double getValues(int pos){
    return values[pos];
  }
  String getUnits(int pos){
    return units[pos];
  }
  void setValuesUnits(int pos, double value, String unit){
    values[pos] = value;
    units[pos] = unit;
    isEmpty = false;
  }
  String getParaName(int pos){
    return paraName[solidType][pos];
  }
  boolean isAngle(int pos){
    return isAngl[solidType][pos];
  }
  public String toString(){
    return solidName[solidType];
  }
  public String getCPP(String name){
    StringBuffer cpp = new StringBuffer("G4"+csgName[solidType]+
                     " *solid" + name + "= new G4" + 
                     csgName[solidType] +"(\"solid"+name+"\""); 
    for (int i=0; i<values.length; i++){
      cpp.append(", "+ values[i]+"*"+units[i]);
    }
    cpp.append(" );\n");
    return cpp.toString();
  }

  private double changMM(int count){
       
       if(getUnits(count).equals("km")){
          mmValue = Math.pow(values[count]*1000000,2.0);
       }else if(getUnits(count).equals("m")){
          mmValue = Math.pow(values[count]*1000,2.0);
       }else if(getUnits(count).equals("cm")){
          mmValue = Math.pow(values[count]*10,2.0);
       }else if(getUnits(count).equals("micrometer")){
          mmValue = Math.pow(values[count]/1000,2.0);
       }else if(getUnits(count).equals("nanoometer")){
          mmValue = Math.pow(values[count]/1000000,2.0);
       }else if(getUnits(count).equals("fermi")){
          mmValue = Math.pow(values[count]/1000000000,2.0);
       }else if(getUnits(count).equals("mm")){
          mmValue = Math.pow(values[count],2.0);
       }else{
          mmValue = 0;
       }
       return mmValue;
  }

  private double Bound(){ 
    bound = 0.0;
    for (int j=0; j<values.length; j++){
      bound += changMM(j);
    }
    return Math.sqrt(bound)*2;
  } 

  public String getPrim(){
   StringBuffer prim = new StringBuffer("#################################\n");
    prim.append("###         GGESolid.prim      ###\n");
    prim.append("#################################\n#\n");

    prim.append("/BoundingBox ");
    prim.append(" -"+Bound()+"   -"+Bound()+"   -"+Bound()+"   ");
    prim.append(Bound()+"   "+Bound()+"   "+Bound()+"\n");

    prim.append("!SetCamera\n!OpenDevice\n!BeginModeling\n\n");
    prim.append("/"+csgDAWNName[solidType]);
    for (int j=0; j<values.length; j++){
       if(getUnits(j).equals("km")){
          prim.append("    "+values[j]*1000000);
       }else if(getUnits(j).equals("m")){
          prim.append("    "+values[j]*1000);
       }else if(getUnits(j).equals("cm")){
          prim.append("    "+values[j]*10);
       }else if(getUnits(j).equals("micrometer")){
          prim.append("    "+values[j]/1000);
       }else if(getUnits(j).equals("nanoometer")){
          prim.append("    "+values[j]/1000000);
       }else if(getUnits(j).equals("fermi")){
          prim.append("    "+values[j]/1000000000);
       }else if(getUnits(j).equals("mm")){
          prim.append("    "+values[j]);
       }else if(getUnits(j).equals("mrad")){
          prim.append(" "+values[j]/1000);
       }else if(getUnits(j).equals("deg")){
          prim.append(" "+values[j]*2*3.14159265358979323/360);
       }else{
          prim.append(" "+values[j]);
       }
    }
    prim.append("\n\n!EndModeling\n!DrawAll\n!CloseDevice");
    return prim.toString();
  }
}

// pX,pY,pZ - The box's half-widths

// pDx1    Half-length along x at the surface positioned at -dz
// pDx2    Half-length along x at the surface positioned at +dz
// pDy1    Half-length along y at the surface positioned at -dz
// pDy2    Half-length along y at the surface positioned at +dz
// pDz     Half-length along z axis

// pRMin   Inside radius
// pRMax   Outside radius
// pDz     half length in z
// pSPhi   The starting phi angle in radians,
//         adjusted such the fSPhi+fDPhi<=2PI,
//         fSPhi>-2PI
// pDPhi   Delta angle of the segment in radians

// pRMin1  inside radius at  -Dz
// pRMin2  inside radius at  +Dz
// pRMax1  outside radius at -Dz
// pRMax2  outside radius at +Dz
// pDz     half length in z
// pSPhi   starting angle of the segment in radians
// pDPhi   delta angle of the segment in radian









