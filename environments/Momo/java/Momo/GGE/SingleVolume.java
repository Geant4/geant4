/*  GGE Single Positioned Volume
 *
 */
// 1998 Sep 29 frame rotation, defaults setting
// 1998 Sep 20 corrected for type 2 and 4 constructors. Yet imcomplete for any inputs (rot axis, mother etc.)
// add Insert Button

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.table.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

class SingleVolume extends JFrame implements ActionListener {
  JTable singleTable;
  private JFrame frame1, frame2;
  private JMenuItem close;
  private JButton create, insert, dell;
  DefaultTableModel singleDataModel;
//  private CreateSingleDialog createSingleDialog;
  private DellSingleDialog dellSingle;
  private InsertSingleDialog insertSingleDialog;
  private String move, name, type, log, many, x0, y0, z0, rotAxis, angle;
  private String pValue;
  private String tableIdent[] = { "Move","pName","pLogic","MomType",
                                  "pMother","pMany","X0","Y0",
                                  "Z0","lengthUnit","RotAxis","Angle","Unit"};
  SingleVolume(){
    super("Single Positioned Placement");
//    createSingleDialog =new CreateSingleDialog(this);
    dellSingle =new DellSingleDialog(this);
    insertSingleDialog = new InsertSingleDialog(this);
    
    JPanel panel = new JPanel();
    panel.setLayout(new BorderLayout(2, 1));
    singleTable = new JTable(singleDataModel = 
                                new DefaultTableModel(tableIdent,0));
    singleTable.setAutoCreateColumnsFromModel(false);
    singleTable.setColumnSelectionAllowed(false);
    singleTable.setRowSelectionAllowed(true);
    singleTable.getTableHeader().setReorderingAllowed(false);
    singleTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

     JComboBox moveCombo = new JComboBox();
       moveCombo.addItem("Frame");
       moveCombo.addItem("Body");
    JComboBox motherCombo = new JComboBox();
       motherCombo.addItem("logical");
       motherCombo.addItem("physical");
       motherCombo.addItem("NULL");
       motherCombo.addActionListener(this);
    JComboBox manyCombo = new JComboBox();
       manyCombo.addItem("true");
       manyCombo.addItem("false");
    JComboBox vectorCombo = new JComboBox();
       vectorCombo.addItem("X");
       vectorCombo.addItem("Y");
       vectorCombo.addItem("Z");
       //vectorCombo.addItem("    "); // this is removed
    LenUnitCombo widthCombo = new LenUnitCombo();
    AngUnitCombo angCombo = new AngUnitCombo();
                 angCombo.addItem("   ");
    TableColumn columns[] = new TableColumn[tableIdent.length];
    for (int i=0; i<tableIdent.length; i++){
       columns[i] = singleTable.getColumn(tableIdent[i]);
    }
    DefaultTableCellRenderer nameCellRenderer = new DefaultTableCellRenderer();
      nameCellRenderer.setBackground(Color.pink);

    columns[0].setCellEditor(new DefaultCellEditor(moveCombo));
    columns[0].setMinWidth(50);

    columns[1].setCellEditor(new NameCellEditor());
    columns[1].setCellRenderer(nameCellRenderer);
    columns[1].sizeWidthToFit();
    columns[1].setMinWidth(60);
    columns[1].setResizable(true);

    columns[2].setCellEditor(new NameCellEditor());
    columns[2].setMinWidth(50);

    columns[3].setCellEditor(new DefaultCellEditor(motherCombo));
    columns[3].setMinWidth(50);

    columns[4].setCellEditor(new NameCellEditor());
    columns[4].setMinWidth(50);

    columns[5].setCellEditor(new DefaultCellEditor(manyCombo));
    columns[5].setMinWidth(50);

    columns[6].setCellEditor(new DoubleCellEditor());
    columns[6].setMinWidth(40);


    columns[7].setCellEditor(new DoubleCellEditor());
    columns[7].setMinWidth(40);


    columns[8].setCellEditor(new DoubleCellEditor());
    columns[8].setMinWidth(40);

    columns[9].setCellEditor(new DefaultCellEditor(widthCombo));
    columns[9].setMinWidth(50);

    columns[10].setCellEditor(new DefaultCellEditor(vectorCombo));
    columns[10].setMinWidth(50);

    columns[11].setCellEditor(new DoubleCellEditor());
    columns[11].setMinWidth(50);

    columns[12].setCellEditor(new DefaultCellEditor(angCombo));
    columns[12].setMinWidth(50);

    JScrollPane scrollpane = new JScrollPane(singleTable);
    scrollpane.setPreferredSize(new Dimension(815, 180));
    scrollpane.setMinimumSize(new Dimension(50, 25));

    JPanel panel2 = new JPanel();
           panel2.setLayout(new GridLayout(1,3));
           panel2.add(create = new JButton(" Append "));
           panel2.add(insert = new JButton(" Insert "));
           panel2.add(dell = new JButton(" Delete "));
    panel.add("Center",scrollpane);
    panel.add("North",panel2);
    getContentPane().add(panel, BorderLayout.CENTER);
    create.addActionListener(this);
    dell.addActionListener(this);
    insert.addActionListener(this);
    setJMenuBar(createMenubar());
    pack();
//    setBackground(Color.green);
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
       dellSingle.setVisible(true);
       dell.setForeground(Color.red);
    }else if(o==create){
       appendSingle();
//       create.setForeground(Color.red);
    }else if(o==insert){
       insertSingleDialog.setVisible(true);
       insert.setForeground(Color.red);
    }else if(o==close){
      setVisible(false);
    }
  }
  public void dellSingleCloseAct(){dell.setForeground(Color.black);repaint();}
  public void createSingleCloseAct(){create.setForeground(Color.black);repaint();}
  public void insertSingleCloseAct(){insert.setForeground(Color.black);repaint();}

  public void singleClear(){
    singleTable.setModel(singleDataModel =
                                new DefaultTableModel(tableIdent,0));
    repaint();
  }

  public void appendSingle(){      // defaults
     Object tmp[] = new Object[15];
      tmp[0] ="Body";
      tmp[1] = "";
      tmp[2] = "";
      tmp[3] = "NULL";
      tmp[4] = "";
      tmp[5] = "false";
      tmp[6] = "0.0";
      tmp[7] = "0.0";
      tmp[8] = "0.0";
      tmp[9] = "mm";
      tmp[10] = "X";
      tmp[11] = "0.0";
      tmp[12] = "deg";
     singleDataModel.addRow(tmp);
     createSingleCloseAct();
     create.setForeground(Color.black);
     repaint();
  }

  public void insertSingle(){       // defaults
     Object tmp[] = new Object[15];
      tmp[0] ="Body";
      tmp[1] = "";
      tmp[2] = "";
      tmp[3] = "physical";
      tmp[4] = "";
      tmp[5] = "false";
      tmp[6] = "0.0";
      tmp[7] = "0.0";
      tmp[8] = "0.0";
      tmp[9] = "mm";
      tmp[10] = "X";
      tmp[11] = "0.0";
      tmp[12] = "deg";

     if(singleTable.getSelectedRow()==-1){
      frame1 = new JFrame();
      JOptionPane.showMessageDialog(frame1, "Choose the pName","Warning Dialog"
, JOptionPane.WARNING_MESSAGE);
     }else{
      int singleRowCount = singleTable.getSelectedRow()+1;
      singleDataModel.insertRow(singleRowCount,tmp);
      insertSingleDialog.setVisible(false);
      insert.setForeground(Color.black);
      repaint();      
     }
  }

  public void dellSingle(){
     if(singleTable.getSelectedRow()==-1){
      frame2 = new JFrame();
      JOptionPane.showMessageDialog(frame2, "Choose the pName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      singleDataModel.removeRow(singleTable.getSelectedRow());
      dellSingle.setVisible(false);
      dell.setForeground(Color.black);
      repaint();
     }
  }

//C++ source code

  String getCPP(){
   StringBuffer cpp = new StringBuffer("\n// Single Positioned Placement \n\n");
   if(singleDataModel.getRowCount()==0){
      return cpp.toString();
   }else{
    for(int i=0;i<singleDataModel.getRowCount();i++){
      String pMother = (String)singleDataModel.getValueAt(i,3); //pMother 
      if(!pMother.equals("NULL")){
       pMother = pMother+(String)singleDataModel.getValueAt(i,4);
      }else if(pMother.equals("NULL")){
       pMother = "NULL";
      } 
// move
     move = singleDataModel.getValueAt(i,0).toString();
// pName      
     name = singleDataModel.getValueAt(i,1).toString();
// pLogical
     log = singleDataModel.getValueAt(i,2).toString(); 
// Mother Type
     type = singleDataModel.getValueAt(i,3).toString(); 
// pMany
     many = singleDataModel.getValueAt(i,5).toString();
// X0
     x0 = singleDataModel.getValueAt(i,6).toString()+"*"+
                                 singleDataModel.getValueAt(i,9).toString();
// Y0
     y0 = singleDataModel.getValueAt(i,7).toString()+"*"+
                                 singleDataModel.getValueAt(i,9).toString();
// Z0 
     z0 = singleDataModel.getValueAt(i,8).toString()+"*"+
                                 singleDataModel.getValueAt(i,9).toString();
// Rotation Axis
     rotAxis = singleDataModel.getValueAt(i,10).toString();

// Angle
     angle = singleDataModel.getValueAt(i,11).toString()+"*"+
                                 singleDataModel.getValueAt(i,12).toString();

// 1st, 3rd Constructors and  2nd, 4th Constructors
     // rotation matrix
  cpp.append( "G4RotationMatrix rotMatrix"+name+";   // unit rotation matrix\n");
  cpp.append( "G4double angle"+name+" = "+angle+";   // rotational angle\n");
  cpp.append( "rotMatrix"+name+".rotate" + rotAxis + "(angle"+name+");  // rot matrix\n\n");


  if(move.equals("Frame")){     // Start 1st,3rd Constructor, i.e., move frame

       cpp.append( "G4VPhysicalVolume *  physical"+name+
                   "= new G4PVPlacement( new G4RotationMatrix(rotMatrix" + name + ") ,        // Frame rotation \n"+
                   "\t\t G4ThreeVector("+x0+", "+y0+", "+z0+"),\n");

	if(type.equals("logical")){ // 3rd Constructor        
               cpp.append("\t\t logical"+log+",\t// 3rd constructor its logical volume \n"+
                          "\t\t \"physical"+name+"\",    //its name \n");
        }else{                      // Start 1st Constructor, mother physical or NULL
               cpp.append("\t\t \"physical"+name+"\",   // 1st constructor its name \n"+
                          "\t\t logical"+log+",         //its logical volume \n");
	}
        cpp.append("\t\t "+pMother+",\t//its mother volume \n"+
                   "\t\t "+many+",\t//no boolean operation \n"+
                   "\t\t 0);\t//copy number \n\n");


  }else if(move.equals("Body")){     // Start 2nd,4th Constructor, move body

        cpp.append( "G4VPhysicalVolume *  physical"+name+
                    "= new G4PVPlacement(G4Transform3D(rotMatrix"+name+",\t//rotation \n"+
                    "\t\t G4ThreeVector("+x0+", "+y0+", "+z0+")),\n");

	if (type.equals("logical")){      // Start 4th Constructor mother = logical
             cpp.append("\t\t logical"+log+",\t//its current logical volume(4th constructor) \n"+
                        "\t\t \"physical"+name+"\",    //its name \n");
         } else {  // mother is  physical or NULL
             cpp.append("\t\t \"physical"+name+"\",   //its name  (2nd constructor)\n"+
                        "\t\t logical"+log+",         //its logical volume \n");
	 }

         cpp.append(
                    "\t\t "+pMother+",              //its mother volume \n"+
                    "\t\t "+many+",                 //no boolean operation \n"+
                    "\t\t 0);                       //copy number \n\n");
      }  // endif
    }  // end for loop
    return cpp.toString();
   }
  }


}














