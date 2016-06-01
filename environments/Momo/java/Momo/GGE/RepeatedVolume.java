/*
 *
 */

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.table.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

class RepeatedVolume extends JFrame implements ActionListener {
  PhysicalItem phyItem;
  JTable pTable;
  private JFrame frame1, frame2;
  private JMenuItem close;
  private JButton create, ins, dell;
  DefaultTableModel pDataModel;
    //  private CreateRepeatDialog createRepeatDialog;
  private DellRepeatDialog dellRepeat;
  private InsertRepeatDialog insertRepeatDialog;
  private String name, log, type, many, x0, y0, z0, dir, step, elem;
  private String pValue;
  private String tableIdent[] = { "pName","pLogical","MomType","pMother",
              "pMany","X0","Y0","Z0","lengthUnit","Dir Of Tlate","StepSize",
              "Unit","nElem" };
  RepeatedVolume(){
    super("Repeated Placement Translation");
    //    createRepeatDialog =new CreateRepeatDialog(this);	
    insertRepeatDialog =new InsertRepeatDialog(this);
    dellRepeat =new DellRepeatDialog(this);
    
    JPanel panel = new JPanel();
//    phyPanel = new PhysicalPanel(parent);
    panel.setLayout(new BorderLayout(2, 1));
    pTable = new JTable(pDataModel = new DefaultTableModel(tableIdent,0));
    pTable.setAutoCreateColumnsFromModel(false);
    pTable.setColumnSelectionAllowed(false);
    pTable.setRowSelectionAllowed(true);
    pTable.getTableHeader().setReorderingAllowed(false);
    pTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    JComboBox manyCombo = new JComboBox();
       manyCombo.addItem("true");
       manyCombo.addItem("false");
    JComboBox vectorCombo = new JComboBox();
       vectorCombo.addItem("X");
       vectorCombo.addItem("Y");
       vectorCombo.addItem("Z");
    JComboBox motherCombo = new JComboBox();
       motherCombo.addItem("logical");
       motherCombo.addItem("physical");
    LenUnitCombo widthCombo = new LenUnitCombo();
    TableColumn columns[] = new TableColumn[tableIdent.length];
    for (int i=0; i<tableIdent.length; i++){
       columns[i] = pTable.getColumn(tableIdent[i]);
    }
    DefaultTableCellRenderer nameCellRenderer = new DefaultTableCellRenderer();
      nameCellRenderer.setBackground(Color.pink);

    columns[0].setCellEditor(new NameCellEditor());
    columns[0].setCellRenderer(nameCellRenderer);

    columns[0].setMinWidth(60);
    columns[0].setResizable(true);

    columns[1].setCellEditor(new NameCellEditor());
    columns[1].setMinWidth(60);

    columns[2].setCellEditor(new DefaultCellEditor(motherCombo));
    columns[2].setMinWidth(80);

    columns[3].setCellEditor(new NameCellEditor());
    columns[3].setMinWidth(50);

    columns[4].setCellEditor(new DefaultCellEditor(manyCombo));
    columns[4].setMinWidth(60);

    columns[5].setCellEditor(new DoubleCellEditor());
    columns[5].setMinWidth(40);


    columns[6].setCellEditor(new DoubleCellEditor());
    columns[6].sizeWidthToFit();
    columns[6].setMinWidth(40);


    columns[7].setCellEditor(new DoubleCellEditor());
    columns[7].sizeWidthToFit();
    columns[7].setMinWidth(40);

    columns[8].setCellEditor(new DefaultCellEditor(widthCombo));
    columns[8].setMinWidth(75);

    columns[9].setCellEditor(new DefaultCellEditor(vectorCombo));
    columns[9].setMinWidth(80);


    columns[10].setCellEditor(new DoubleCellEditor());
    columns[10].setMinWidth(55);

    columns[11].setCellEditor(new DefaultCellEditor(widthCombo));
    columns[11].setMinWidth(65);

    columns[12].setCellEditor(new NameCellEditor());
    columns[12].setMinWidth(50);

    JScrollPane scrollpane = new JScrollPane(pTable);
    scrollpane.setPreferredSize(new Dimension(815, 180));
    scrollpane.setMinimumSize(new Dimension(50, 25));
    JPanel panel2 = new JPanel();
      panel2.setLayout(new GridLayout(1,3));
      panel2.add(create = new JButton(" Append "));
      panel2.add(ins = new JButton(" Insert "));
      panel2.add(dell = new JButton(" Delete "));

    panel.add("Center",scrollpane);
    panel.add("North",panel2);
    getContentPane().add(panel, BorderLayout.CENTER);
    create.addActionListener(this);
    dell.addActionListener(this);
    ins.addActionListener(this);

    PhysicalItem phyItem = new PhysicalItem();
    setJMenuBar(createMenubar());
    pack();
    setResizable(false);
  }
  private JMenuBar createMenubar(){
    JMenu file;
    JMenuBar mb = new JMenuBar();
    file = new JMenu("File");
//    file.add(save = new JMenuItem("Save"));
    file.add(close = new JMenuItem("Close"));
    mb.add(file);
    close.addActionListener(this);
    return mb;
  }
  
  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();
    if(o==dell){
       dellRepeat.setVisible(true);
       dell.setForeground(Color.red);
    }else if(o==ins){
       insertRepeatDialog.setVisible(true);
       ins.setForeground(Color.red);
    }else if(o==create){
      appendRepeat();
//      create.setForeground(Color.red);
    }else if(o==close){
      setVisible(false);
    }
  }
  public void dellRepeatCloseAct(){dell.setForeground(Color.black);repaint();}
  public void createRepeatCloseAct(){create.setForeground(Color.black);repaint();}
  public void insertRepeatCloseAct(){ins.setForeground(Color.black);repaint();}

  public void repeatedClear(){
     pTable.setModel(pDataModel = new DefaultTableModel(tableIdent,0));
     repaint();
  }

  public void appendRepeat(){
     Object tmp[] = new Object[15];
      tmp[0] = "";
      tmp[1] = "";
      tmp[2] ="physical";
      tmp[3] = "";
      tmp[4] = "false";
      tmp[5] = "0.0";
      tmp[6] = "0.0";
      tmp[7] = "0.0";
      tmp[8] = "mm";
      tmp[9] = "X";
      tmp[10] = "";
      tmp[11] = "mm";
      tmp[12] = "";
     pDataModel.addRow(tmp);
  }

  public void insertRepeat(){
     Object tmp[] = new Object[15];
      tmp[0] = "";
      tmp[1] = "";
      tmp[2] ="physical";
      tmp[3] = "";
      tmp[4] = "false";
      tmp[5] = "0.0";
      tmp[6] = "0.0";
      tmp[7] = "0.0";
      tmp[8] = "mm";
      tmp[9] = "X";
      tmp[10] = "";
      tmp[11] = "mm";
      tmp[12] = "";


     if(pTable.getSelectedRow()==-1){
      frame1 = new JFrame();
      JOptionPane.showMessageDialog(frame1, "Choose the pName","Warning Dialog"
, JOptionPane.WARNING_MESSAGE);
     }else{
      int pRowCount = pTable.getSelectedRow()+1;
      pDataModel.insertRow(pRowCount,tmp);
     }
  }

  public void dellRepeat(){
     if(pTable.getSelectedRow()==-1){
      frame2 = new JFrame();
      JOptionPane.showMessageDialog(frame2, "Choose the pName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      pDataModel.removeRow(pTable.getSelectedRow());
      dellRepeat.setVisible(false);
      dell.setForeground(Color.black);
      repaint();
     }
  }

  //C++ source code
  String getCPP(){

   StringBuffer cpp = new StringBuffer("\n// Repeated Placement Translation \n\n");

   if(pDataModel.getRowCount()==0){
      return cpp.toString();
   }else{

    for(int i=0;i<pDataModel.getRowCount();i++){
       String pMother = (String)pDataModel.getValueAt(i,2)+
                             (String)pDataModel.getValueAt(i,3);      
// pName
      name = pDataModel.getValueAt(i,0).toString();
// pLogical
      log = pDataModel.getValueAt(i,1).toString();
// Mother Type
      type = pDataModel.getValueAt(i,2).toString();
// pMany
      many = pDataModel.getValueAt(i,4).toString();
// X0
      x0 = pDataModel.getValueAt(i,5).toString()+"*"+
                             pDataModel.getValueAt(i,8).toString();
// Y0
      y0 = pDataModel.getValueAt(i,6).toString()+"*"+
                             pDataModel.getValueAt(i,8).toString();
// Z0
      z0 = pDataModel.getValueAt(i,7).toString()+"*"+
                             pDataModel.getValueAt(i,8).toString();
// Dir Of Tlate
      dir = pDataModel.getValueAt(i,9).toString();
// Step Size
      step = pDataModel.getValueAt(i,10).toString()+"*"+
                             pDataModel.getValueAt(i,11).toString();
// nElement
      elem = pDataModel.getValueAt(i,12).toString();

     if(dir.equals("X")){
       cpp.append("G4int copy"+name+";\ncopy"+name+"=0;\n");

       cpp.append("for (G4int "+name+"=1; "+name+"<="+elem+"; "+name+"++){\n");
       cpp.append("  G4double trans"+name+" ="+x0+"+"+step+"*("+name+"-1);\n");

         if(type.equals("physical")){

              cpp.append("  G4VPhysicalVolume * "+"physical"+name+
                        " = new G4PVPlacement(0,      //no rotation \n"+
                        "\t\t G4ThreeVector(trans"+name+", "+y0+", "+z0+"),\n"+
                        "\t\t \"physical"+name+"\",   //its name \n"+
                        "\t\t logical"+log+",         //its logical volume\n"+
                        "\t\t "+pMother+",            //its mother volume\n"+
                        "\t\t "+many+",             //no boolean operation\n"+
                        "\t\t copy"+name+"++);      //copy number \n"+
                        "}\n\n");

         }else if(type.equals("logical")){//mother logical

              cpp.append("  G4VPhysicalVolume * "+"physical"+name+
                        " = new G4PVPlacement(0,      //no rotation \n"+
                        "\t\t G4ThreeVector(trans"+name+", "+y0+", "+z0+"),\n"+
                        "\t\t logical"+log+",         //its logical volume\n"+
                        "\t\t \"physical"+name+"\",   //its name\n"+
                        "\t\t "+pMother+",            //its mother volume\n"+
                        "\t\t "+many+",             //no boolean operation\n"+
                        "\t\t copy"+name+"++);      //copy number \n"+
                        "}\n\n");
	 } // mother logical end

     }else if(dir.equals("Y")){

      cpp.append("G4int copy"+name+";\ncopy"+name+"=0;\n");

      cpp.append("for (G4int "+name+"=1; "+name+"<="+elem+"; "+name+"++){\n");

      cpp.append("  G4double trans"+name+" = "+y0+"+"+step+"*("+name+"-1);\n");

	  if(type.equals("physical")){

              cpp.append("  G4VPhysicalVolume * "+"physical"+name+
                        " = new G4PVPlacement(0,       //no rotation \n"+
                        "\t\t G4ThreeVector("+x0+", trans"+name+", "+z0+"),\n"+
                        "\t\t \"physical"+name+"\",    //its name \n"+ 
                        "\t\t logical"+log+",          //its logical volume\n"+
                        "\t\t "+pMother+",             //its mother volume\n"+
                        "\t\t "+many+",             //no boolean operation\n"+
                        "\t\t copy"+name+"++);      //copy number \n"+
                        "}\n\n");

	  }else if(type.equals("logical")){// mother logical

              cpp.append("  G4VPhysicalVolume * "+"physical"+name+
                        " = new G4PVPlacement(0,       //no rotation \n"+
                        "\t\t G4ThreeVector("+x0+", trans"+name+", "+z0+"),\n"+
                        "\t\t logical"+log+",          //its logical volume\n"+
                        "\t\t \"physical"+name+"\",    //its name \n"+
                        "\t\t "+pMother+",             //its mother volume\n"+
                        "\t\t "+many+",             //no boolean operation\n"+
                        "\t\t copy"+name+"++);      //copy number\n"+
                        "}\n\n");

	  }//mother logical end

     }else if(dir.equals("Z")){

      cpp.append("G4int copy"+name+";\ncopy"+name+"=0;\n");

      cpp.append("for (G4int "+name+"=1; "+name+"<="+elem+";"+name+"++){\n");

      cpp.append("  G4double trans"+name+" = "+z0+"+"+step+"*("+name+"-1);\n");

         if(type.equals("physical")){

           cpp.append("  G4VPhysicalVolume * "+"physical"+name+
                    " = new G4PVPlacement(0,           //no rotation \n"+
                    "\t\t G4ThreeVector("+x0+", "+y0+", trans"+name+"),\n"+
                    "\t\t \"physical"+name+"\",        //its name \n"+
                    "\t\t logical"+log+",              //its logical volume\n"+
                    "\t\t "+pMother+",                 //its mother volume\n"+
                    "\t\t "+many+",                 //no boolean operation\n"+
                    "\t\t copy"+name+"++);          //copy number\n"+
                    "}\n");

         }else if(type.equals("logical")){// mother logical
           cpp.append("  G4VPhysicalVolume * "+"physical"+name+
                    " = new G4PVPlacement(0,           //no rotation \n"+
                    "\t\t G4ThreeVector("+x0+", "+y0+", trans"+name+"),\n"+
                    "\t\t logical"+log+",              //its logical volume\n"+
                    "\t\t \"physical"+name+"\",        //its name \n"+
                    "\t\t "+pMother+",                 //its mother volume\n"+
                    "\t\t "+many+",                 //no boolean operation\n"+
                    "\t\t copy"+name+"++);          //copy number\n"+
                    "}\n\n");
         }//mother logical end


     }
    }
    return cpp.toString();
   }
  }
}














