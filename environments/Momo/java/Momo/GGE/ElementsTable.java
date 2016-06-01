// PERIODIC TABLE OF THE ELEMENTS
// PDG table 4.1

import java.awt.*;
import com.sun.java.swing.*;

class ElementsTable extends AbstractGGEFrame {
  private ElementCell ecel[];
  ElementsTable(){
    super("Elements");
    GridBagLayout gbl = new GridBagLayout();
    GridBagConstraints gbc = new GridBagConstraints();
    getContentPane().setLayout(gbl);
    getContentPane().setBackground(new Color(230, 250, 230));
    int elemPosX[] = {1, 0, 17,
		      0, 1, 12, 13, 14, 15, 16, 17,
		      0, 1, 12, 13, 14, 15, 16, 17,
		      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
		      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
		      0, 1,
		      3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
		      3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
		      0, 1,
		      3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 };
    int elemPosY[] = {0, 0, 0,
		      1, 1, 1, 1, 1, 1, 1, 1,
		      2, 2, 2, 2, 2, 2, 2, 2,
		      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
		      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		      5, 5, 
		      7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
		      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		      6, 6,
		      8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };

    ElementItem el[] = makeElements();
    gbc.insets = new Insets(1,1,1,1);
    gbc.fill = GridBagConstraints.BOTH;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    ecel = new ElementCell[el.length];
    for (int i = 0; i<el.length; i++){
      gbc.gridx = elemPosX[i]; gbc.gridy = elemPosY[i];
      ecel[i] = new ElementCell(el[i]);
      gbl.setConstraints(ecel[i], gbc); getContentPane().add(ecel[i]);
    }
    Color lightRed = new Color(255, 200, 200);
    Color lightYellow = new Color(255, 255, 96);
    JLabel lan = new JLabel("L", JLabel.CENTER);
    lan.setOpaque(true);
    lan.setBackground(lightRed);
    gbc.gridx = 2; gbc.gridy = 5;
    gbl.setConstraints(lan, gbc); getContentPane().add(lan);
    JLabel act = new JLabel("A", JLabel.CENTER);
    act.setOpaque(true);
    act.setBackground(lightRed);
    gbc.gridx = 2; gbc.gridy = 6;
    gbl.setConstraints(act, gbc); getContentPane().add(act);
    JLabel lt = new JLabel("PERIODIC TABLE OF THE ELEMENTS", JLabel.CENTER);
    lt.setOpaque(true);
    lt.setBackground(Color.white);
    gbc.gridx = 2; gbc.gridy = 0;
    gbc.gridwidth = 10; gbc.gridheight = 2;
    gbl.setConstraints(lt, gbc); getContentPane().add(lt);
    JLabel lla = new JLabel("Lanthanides", JLabel.CENTER);
    lla.setOpaque(true);
    lla.setBackground(lightRed);
    gbc.gridx = 0; gbc.gridy = 7;
    gbc.gridwidth = 3; gbc.gridheight = 1;
    gbl.setConstraints(lla, gbc); getContentPane().add(lla);
    JLabel lac = new JLabel("Actinides", JLabel.CENTER);
    lac.setOpaque(true);
    lac.setBackground(lightRed);
    gbc.gridx = 0; gbc.gridy = 8;
    gbl.setConstraints(lac, gbc); getContentPane().add(lac);

    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.anchor = GridBagConstraints.SOUTH;
    gbc.gridwidth = 1; gbc.gridheight = 1;
    JLabel num[] = new JLabel[15];
    for (int i=0; i<15; i++){
      num[i]  = new JLabel(i+3+"", JLabel.CENTER);
      num[i].setBackground(lightYellow);
      num[i].setOpaque(true);
      gbc.gridx = 2+i; 
      gbc.gridy = i < 10 ? 2:0 ;
      gbl.setConstraints(num[i], gbc); getContentPane().add(num[i]);
    }
    pack();
    setResizable(false);
//kashika zokusei
//    setVisible(true);
  }
  ElementItem[] getSelectedElements(){
    int k = 0;
    if (ElementCell.selectedElems == 0) return null;
    ElementItem ei[] = new ElementItem[ElementCell.selectedElems];
    for (int i=0; i<ecel.length; i++){
      if(ecel[i].isSelected()){
	ei[k] = ecel[i].ei;
	ecel[i].clear();
	k++;
      }
    }
    return ei;
  }
  ElementItem[] makeElements(){
    ElementItem el[] = new ElementItem[104];
    el[0] = new ElementItem("Deuterium","D",1,2.014, false);
    el[1] = new ElementItem("Hydrogen","H",1,1.00794, false);
    el[2] = new ElementItem("Helium","He",2,4.002602, false);
    el[3] = new ElementItem("Lithium","Li",3,6.941, true);
    el[4] = new ElementItem("Beryllium","Be",4,9.012182, true);
    el[5] = new ElementItem("Boron","B",5,10.811, false);
    el[6] = new ElementItem("Carbon","C",6,12.011, false);
    el[7] = new ElementItem("Nitrogen","N",7,14.00674, false);
    el[8] = new ElementItem("Oxygen","O",8,15.9994, false);
    el[9] = new ElementItem("Fluorine","F",9,18.9984032, false);
    el[10] = new ElementItem("Neon","Ne",10,20.1797, false);
    el[11] = new ElementItem("Sodium","Na",11,22.989768, true);
    el[12] = new ElementItem("Magnesium","Mg",12,24.305, true);
    el[13] = new ElementItem("Aluminum","Al",13,26.981539, true);
    el[14] = new ElementItem("Silicon","Si",14,28.0855, false);
    el[15] = new ElementItem("Phosphorus","P",15,30.973762, false);
    el[16] = new ElementItem("Sulfur","S",16,32.066, false);
    el[17] = new ElementItem("Chlorine","Cl",17,35.4527, false);
    el[18] = new ElementItem("Argon","Ar",18,39.948, false);
    el[19] = new ElementItem("Potassium","K",19,39.0983, true);
    el[20] = new ElementItem("Calcium","Ca",20,40.078, true);
    el[21] = new ElementItem("Scandium","Sc",21,44.95591, true);
    el[22] = new ElementItem("Titanium","Ti",22,47.867, true);
    el[23] = new ElementItem("Vanadium","V",23,50.9415, true);
    el[24] = new ElementItem("Chromium","Cr",24,51.9961, true);
    el[25] = new ElementItem("Manganese","Mn",25,54.93805, true);
    el[26] = new ElementItem("Iron","Fe",26,55.845, true);
    el[27] = new ElementItem("Cobalt","Co",27,58.9332, true);
    el[28] = new ElementItem("Nickel","Ni",28,58.69, true);
    el[29] = new ElementItem("Copper","Cu",29,63.546, true);
    el[30] = new ElementItem("Zinc","Zn",30,65.39, true);
    el[31] = new ElementItem("Gallium","Ga",31,69.723, true);
    el[32] = new ElementItem("Germanium","Ge",32,72.61, true);
    el[33] = new ElementItem("Arsenic","As",33,74.92159, false);
    el[34] = new ElementItem("Selenium","Se",34,78.96, false);
    el[35] = new ElementItem("Bromine","Br",35,79.904, false);
    el[36] = new ElementItem("Krypton","Kr",36,83.8, false);
    el[37] = new ElementItem("Rubidium","Rb",37,85.4678, true);
    el[38] = new ElementItem("Strontium","Sr",38,87.62, true);
    el[39] = new ElementItem("Yttrium","Y",39,88.90585, true);
    el[40] = new ElementItem("Zirconium","Zr",40,91.224, true);
    el[41] = new ElementItem("Niobium","Nb",41,92.90638, true);
    el[42] = new ElementItem("Molybdenum","Mo",42,95.94, true);
    el[43] = new ElementItem("Technetium","Tc",43,98.91, true);
    el[44] = new ElementItem("Ruthenium","Ru",44,101.07, true);
    el[45] = new ElementItem("Rhodium","Rh",45,102.9055, true);
    el[46] = new ElementItem("Palladium","Pd",46,106.42, true);
    el[47] = new ElementItem("Silver","Ag",47,107.8682, true);
    el[48] = new ElementItem("Cadmium","Cd",48,112.411, true);
    el[49] = new ElementItem("Indium","In",49,114.82, true);
    el[50] = new ElementItem("Tin","Sn",50,118.71, true);
    el[51] = new ElementItem("Antimony","Sb",51,121.75, true);
    el[52] = new ElementItem("Tellurium","Te",52,127.6, false);
    el[53] = new ElementItem("Iodine","I",53,126.90447, false);
    el[54] = new ElementItem("Xenon","Xe",54,131.29, false);
    el[55] = new ElementItem("Cesium","Cs",55,132.90543, true);
    el[56] = new ElementItem("Barium","Ba",56,137.327, true);
    el[57] = new ElementItem("Lanthanum","La",57,138.9055, true);
    el[58] = new ElementItem("Cerium","Ce",58,140.115, true);
    el[59] = new ElementItem("Praseodymium","Pr",59,140.90765, true);
    el[60] = new ElementItem("Neodymium","Nd",60,144.24, true);
    el[61] = new ElementItem("Promethium","Pm",61,145, true);
    el[62] = new ElementItem("Samarium","Sm",62,150.36, true);
    el[63] = new ElementItem("Europium","Eu",63,151.965, true);
    el[64] = new ElementItem("Gadolinium","Gd",64,157.25, true);
    el[65] = new ElementItem("Terbium","Tb",65,158.92534, true);
    el[66] = new ElementItem("Dysprosium","Dy",66,162.5, true);
    el[67] = new ElementItem("Holmium","Ho",67,164.93032, true);
    el[68] = new ElementItem("Erbium","Er",68,167.26, true);
    el[69] = new ElementItem("Thulium","Tm",69,168.93421, true);
    el[70] = new ElementItem("Ytterbium","Yb",70,173.04, true);
    el[71] = new ElementItem("Lutetium","Lu",71,174.967, true);
    el[72] = new ElementItem("Hafnium","Hf",72,178.49, true);
    el[73] = new ElementItem("Tantalum","Ta",73,180.9479, true);
    el[74] = new ElementItem("Tungsten","W",74,183.85, true);
    el[75] = new ElementItem("Rhenium","Re",75,186.207, true);
    el[76] = new ElementItem("Osmium","Os",76,190.2, true);
    el[77] = new ElementItem("Iridium","Ir",77,192.22, true);
    el[78] = new ElementItem("Platinum","Pt",78,195.08, true);
    el[79] = new ElementItem("Gold","Au",79,196.96654, true);
    el[80] = new ElementItem("Mercury","Hg",80,200.59, true);
    el[81] = new ElementItem("Thallium","Tl",81,204.3833, true);
    el[82] = new ElementItem("Lead","Pb",82,207.2, true);
    el[83] = new ElementItem("Bismuth","Bi",83,208.98037, true);
    el[84] = new ElementItem("Polonium","Po",84,209, true);
    el[85] = new ElementItem("Astatine","At",85,210, false);
    el[86] = new ElementItem("Radon","Rn",86,222, false);
    el[87] = new ElementItem("Francium","Fr",87,223, true);
    el[88] = new ElementItem("Radium","Ra",88,226.025, true);
    el[89] = new ElementItem("Actinium","Ac",89,227.028, true);
    el[90] = new ElementItem("Thorium","Th",90,232.0381, true);
    el[91] = new ElementItem("Protactinium","Pa",91,231.03588, true);
    el[92] = new ElementItem("Uranium","U",92,238.0289, true);
    el[93] = new ElementItem("Neptunium","Np",93,237.048, true);
    el[94] = new ElementItem("Plutonium","Pu",94,244, true);
    el[95] = new ElementItem("Americium","Am",95,243, true);
    el[96] = new ElementItem("Curium","Cm",96,247, true);
    el[97] = new ElementItem("Berkelium","Bk",97,247, true);
    el[98] = new ElementItem("Californium","Cf",98,251, true);
    el[99] = new ElementItem("Einsteinium","Es",99,254, true);
    el[100] = new ElementItem("Fermium","Fm",100,257, true);
    el[101] = new ElementItem("Mendelevium","Md",101,258, true);
    el[102] = new ElementItem("Nobelium","No",102,259, true);
    el[103] = new ElementItem("Lawrencium","Lr",103,260, true);
    return el;
  }
}







