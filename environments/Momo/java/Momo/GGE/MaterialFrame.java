/*
 *
 */
//7.3


import com.sun.java.swing.*;
import com.sun.java.swing.border.*;
import com.sun.java.swing.event.*;
import com.sun.java.swing.table.*;
import com.sun.java.swing.preview.*;
import com.sun.java.swing.preview.filechooser.*;
import ExampleFileFilter;
import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;


public class MaterialFrame extends AbstractGGEFrame implements ActionListener {
  ElementsTable et;
  DellCombiDialog dellCombi;
  DellScratchDialog dellScratch;
  MaterialScratchTable msTable;
  MaterialCombiTable mcTable;
  DefaultTableModel msDataModel, mcDataModel;

  String msNames[] = {"Use","Name","A","Z","Density","Unit",
				"State","Temp"," Unit","Press","Unit "};
  String mcNames[] = {"Use","Name","Elements","Density","Unit",
				"State","Temp"," Unit","Press","Unit "};

  JLabel label;
  private JMenuItem load, save, append, clear, tOpen, tClose;
  private JFileChooser saveData, loadData;
  private String fileName, createFile, newFile;
  static String fileDir;

  /* For OptionPane */
  private JFrame frame1,frame2,frame3,frame4,frame5,frame6,frame7,frame8;

  public MaterialFrame() {
      super("Materials");
      et = new ElementsTable();

    msTable = new MaterialScratchTable(msDataModel = 
				new DefaultTableModel(msNames,0), this);
    mcTable = new MaterialCombiTable(mcDataModel = 
				new DefaultTableModel(mcNames,0), this);
    getContentPane().setLayout(new BorderLayout());

    JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT,true,
                 msTable.makeMaterialPanel(), mcTable.makeMaterialPanel());

    getContentPane().add(splitPane, "Center");
    setJMenuBar( createMenubar() );
    setSize(700,580);
    pack();
    setVisible(true);

  }
  private JMenuBar createMenubar(){
      JMenu file, element;
      JMenuBar mb = new JMenuBar();

      file = new JMenu("File");
         file.add(load = new JMenuItem("Load a Material File"));
         file.addSeparator();
         file.add(append = new JMenuItem("Append a Material File"));
         file.addSeparator();
         file.add(save = new JMenuItem("Save Material"));
         file.addSeparator();
         file.add(clear = new JMenuItem("Clear Material"));

      load.addActionListener(this);
      append.addActionListener(this);
      save.addActionListener(this);
      clear.addActionListener(this);

      mb.add(file);

      element = new JMenu("Elements");
         element.add(tOpen = new JMenuItem("Open Periodic Table"));
         element.addSeparator();
         element.add(tClose = new JMenuItem("Close Periodic Table"));
      tOpen.addActionListener(this);
      tClose.addActionListener(this);

      mb.add(element);
      return mb;
  }

  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();

    if( o == load ){MaterialfileLoad();}
    if( o == append ){MaterialfileAppend();}
    if( o == save ){MaterialfileSave();}
    if( o == clear ){MaterialfileClear();}
    if( o == tOpen ){
               et.setVisible(true);
    }
    if( o == tClose ){
               et.setVisible(false);
    }
  }

  private synchronized void resume(){
    notify();
  }

//DATA LOAD ,SAVE -------------------------------------------------
//Append
  void MaterialfileAppend(){
      JFileChooser loadData = new JFileChooser(".");
      loadData.cancelSelection();
      loadData.resetChoosableFileFilters();

     ExampleFileFilter mtFile = new ExampleFileFilter("g4mt","MaterialSource");
      loadData.setFileFilter(mtFile);
      loadData.setDialogTitle("Load Material");
      loadData.setMultiSelectionEnabled(true);
      if (loadData.showOpenDialog(this) == -1 ) return;
      File fml = loadData.getSelectedFile();
      fileName = fml.getPath();
      if (fml.isFile()) {
        resume();
        try{
          FileInputStream filein = new FileInputStream(fileName);
          ObjectInputStream objin =new ObjectInputStream(filein);
          Vector mtData = (Vector)objin.readObject();
          Vector msData =  (Vector)mtData.firstElement();
          Vector mcData =  (Vector)mtData.lastElement();

          if(msDataModel.getRowCount()==0 && msData.size()==0){
             //--patern 1--
          }else{
             //--patern 2--
             for(int i=0; i<msData.size(); i++){
               Vector ms = (Vector)msData.elementAt(i);
               msDataModel.addRow(ms);
             }
             repaint();
          }

          if(mcDataModel.getRowCount()==0 && mcData.size()==0){
              //--patern 3--
          }else{
              //--patern 4--
              for(int i=0; i<mcData.size(); i++){
               Vector mc = (Vector)mcData.elementAt(i);
               mcDataModel.addRow(mc);
              }
              repaint();
          }
          objin.close();
          repaint();
        }catch(Exception e){
          System.out.println(e.toString());
        }
      }else{
        System.out.println("error="+fml.toString()+"isNotaFile");
      }
      return;
  }


// Load
  void MaterialfileLoad(){
      JFileChooser loadData = new JFileChooser(".");
      loadData.cancelSelection();
      loadData.resetChoosableFileFilters();

     ExampleFileFilter mtFile = new ExampleFileFilter("g4mt","MaterialSource");
      loadData.setFileFilter(mtFile);
      loadData.setDialogTitle("Load Material");
      loadData.setMultiSelectionEnabled(true);
      if (loadData.showOpenDialog(this) == -1 ) return;
      File fml = loadData.getSelectedFile();
      fileName = fml.getPath();
      if (fml.isFile()) {
        resume();
        try{
          FileInputStream filein = new FileInputStream(fileName);
          ObjectInputStream objin =new ObjectInputStream(filein);
          Vector mtData = (Vector)objin.readObject();
          Vector msData =  (Vector)mtData.firstElement();
          Vector mcData =  (Vector)mtData.lastElement();

          if(msDataModel.getRowCount()==0 && msData.size()==0){
             //--patern 1--
          }else{
             //--patern 2--
             msTable.setModel(msDataModel = new DefaultTableModel(msNames,0));
             for(int i=0; i<msData.size(); i++){
               Vector ms = (Vector)msData.elementAt(i);
               msDataModel.addRow(ms);
             }
             repaint();
          }

          if(mcDataModel.getRowCount()==0 && mcData.size()==0){
              //--patern 3--
          }else{
              //--patern 4--
              mcTable.setModel(mcDataModel = new DefaultTableModel(mcNames,0));
              for(int i=0; i<mcData.size(); i++){
               Vector mc = (Vector)mcData.elementAt(i);
               mcDataModel.addRow(mc);
              }
              repaint();
          }
          objin.close();
          repaint();
        }catch(Exception e){
          System.out.println(e.toString());
        }
      }else{
        System.out.println("error="+fml.toString()+"isNotaFile");
      }
      return;
  }
/* --Save-------------------------------- */
  void MaterialfileSave(){
     JFileChooser saveData = new JFileChooser(".");
     ExampleFileFilter mtFile = new ExampleFileFilter("g4mt","MaterialSource");
     saveData.setFileFilter(mtFile);
     saveData.setDialogTitle("Save Material");
     if ( saveData.showSaveDialog(this) == -1 ) return;
     File fms = saveData.getSelectedFile();
     fileName = fms.getPath();
     if (fileName != null){
       try{
         FileOutputStream fileout = new FileOutputStream(fileName);
         ObjectOutputStream objout = new ObjectOutputStream(fileout);
         Vector mtVector = new Vector();
         Vector msVector = new Vector();
//         Object mtsave[] = new Object[1];
//         mtsave[0] = "";
//         for (int i=0; i<msTable.getRowCount(); i++){
//          msTable.setValueAt(mtsave[0],i,0);
//         }
//         for (int i=0; i<mcTable.getRowCount(); i++){
//          mcTable.setValueAt(mtsave[0],i,0);
//         }
         mtVector.addElement((Object)msDataModel.getDataVector());
         mtVector.addElement((Object)mcDataModel.getDataVector());
         objout.writeObject(mtVector);
         objout.flush();
         fileout.close();
       }catch(IOException e){ System.out.println(e.getMessage()); }
         return;
     }
     if ( fileName == null){
       return;
    }
  }

//Clear
  void MaterialfileClear(){
    msTable.setModel(msDataModel = new DefaultTableModel(msNames,0));
    mcTable.setModel(mcDataModel = new DefaultTableModel(mcNames,0));
  }

/* Material APPEND , INSERT , DELETE(START)------------------------- 8.10 */
//----------------Start the Append 8.10
  void appendMC(){     
   ElementItem[] eic = et.getSelectedElements();
   if (eic == null){
     frame1 = new JFrame();
     JOptionPane.showMessageDialog(frame1, "Choose Elements by Opening the PeriodicTable","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     return;
   }else{
     Object tmp[] = new Object[10];
      tmp[0] = "";
      tmp[1] = "";
      tmp[2] = new RatioItem(eic);
      tmp[3] = "";
      tmp[4] = "g/cm3";
      tmp[5] = "kStateUndefined";
      tmp[6] = "273.15";
      tmp[7] = "kelvin";
      tmp[8] = "1.0";
      tmp[9] = "atmosphere";
     mcDataModel.addRow(tmp);
   }  
  }
  void appendMS(){
   ElementItem[] eis = et.getSelectedElements();
    if (eis == null){
     frame2 = new JFrame();
     JOptionPane.showMessageDialog(frame2, "Choose a Element by Opening the Periodic Table","Warning Dialog", JOptionPane.WARNING_MESSAGE);

       return;
    }else{
     Object tmp[] = new Object[11];
      tmp[0] = "";
      tmp[1] = new String(eis[0].name);
      tmp[2] = new Integer(eis[0].atomNum);
      tmp[3] = new Double(eis[0].massNum);
      tmp[4] = "";
      tmp[5] = "g/cm3";
      tmp[6] = "kStateUndefined";
      tmp[7] = "273.15";
      tmp[8] = "kelvin";
      tmp[9] = "1.0";
      tmp[10] = "atmosphere";
     msDataModel.addRow(tmp);
    }
  }

// ------------ Start the Insert Material 8.10 ---------
  void insertMS(){
   ElementItem[] eis = et.getSelectedElements();
    if (eis == null){
     frame3 = new JFrame();
     JOptionPane.showMessageDialog(frame3, "Choose a Element by Opening the Periodic Table","Warning Dialog", JOptionPane.WARNING_MESSAGE);

       return;
    }else{
     Object tmp[] = new Object[11];
      tmp[0] = "";
      tmp[1] = new String(eis[0].name);
      tmp[2] = new Integer(eis[0].atomNum);
      tmp[3] = new Double(eis[0].massNum);
      tmp[4] = "";
      tmp[5] = "g/cm3";
      tmp[6] = "kStateUndefined";
      tmp[7] = "273.15";
      tmp[8] = "kelvin";
      tmp[9] = "1.0";
      tmp[10] = "atmosphere";
     if(msTable.getSelectedRow()==-1){
      frame4 = new JFrame();
      JOptionPane.showMessageDialog(frame4, "Choose the MaterialName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      int msRowCount = msTable.getSelectedRow()+1;
      msDataModel.insertRow(msRowCount,tmp);
     }
    }
  }
  void insertMC(){
   ElementItem[] eic = et.getSelectedElements();
   if (eic == null){
     frame5 = new JFrame();
     JOptionPane.showMessageDialog(frame5, "Choose Elements by Opening the Periodic Table","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     return;
   }else{
     Object tmp[] = new Object[10];
      tmp[0] = "";
      tmp[1] = "";
      tmp[2] = new RatioItem(eic);
      tmp[3] = "";
      tmp[4] = "g/cm3";
      tmp[5] = "kStateUndefined";
      tmp[6] = "273.15";
      tmp[7] = "kelvin";
      tmp[8] = "1.0";
      tmp[9] = "atmosphere";
     if(mcTable.getSelectedRow()==-1){
      frame6 = new JFrame();
      JOptionPane.showMessageDialog(frame6, "Choose the MaterialName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      int mcRowCount = mcTable.getSelectedRow()+1;
      mcDataModel.insertRow(mcRowCount,tmp);
     }
   }
  }
//----------------------DELETE (Start DELL 8.10)------------------
  void dellMS(){
     if(msTable.getSelectedRow()==-1){
      frame7 = new JFrame();
      JOptionPane.showMessageDialog(frame7, "Choose the MaterialName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      msDataModel.removeRow(msTable.getSelectedRow());
      dellScratch = new DellScratchDialog(this);
      msTable.dellMatCloseAct();      
      repaint();
      dellScratch.setVisible(false);
     }   
  }
  void dellMC(){
     if(mcTable.getSelectedRow()==-1){
      frame8 = new JFrame();
      JOptionPane.showMessageDialog(frame8, "Choose the MaterialName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      mcDataModel.removeRow(mcTable.getSelectedRow());
      dellCombi = new DellCombiDialog(this);
      mcTable.dellMatCloseAct();
      repaint();
      dellCombi.setVisible(false);
     }
  } 
/* Material APPEND , INSERT , DELETE(END)--------------------------- 8.10 */


//Data Load Method
  void setLoadData(Vector materialData){
    GGEItem item;
  }

  void scratchSelected(){
    mcTable.stopEditing();
  }
  void combiSelected(){
    msTable.stopEditing();
  }
  public void windowIconified(WindowEvent event){
    et.setVisible(false);
    msTable.stopEditing();
    mcTable.stopEditing();
  }
  public void windowDeiconified(WindowEvent event){
    et.setVisible(true);
  }

//------- C++ source code ------------------------------------
  String getCPP(){
    StringBuffer smatcpp = new StringBuffer("\n// Materials from Scratch\n\n");
    Vector data = msDataModel.getDataVector();
    Vector rowData;

// Density must be checked before generating C++ code. 8.13 not yet!!

    for (int i=0; i<data.size(); i++){
      rowData = (Vector)data.elementAt(i);
      if(!rowData.elementAt(0).equals("Used")){
//       System.out.println("Next");
//Use G4 Units
      }else if(rowData.elementAt(0).equals("Used")){
       smatcpp.append("G4Material* " + rowData.elementAt(1) 
		+ " = new G4Material(\"" + rowData.elementAt(1)
		+ "\", " +rowData.elementAt(2) 
		+ ", " + rowData.elementAt(3) + "*g/mole, " 
                + rowData.elementAt(4)+"*" +rowData.elementAt(5)+","
                +rowData.elementAt(6));
       if (!rowData.elementAt(7).toString().equals("")){
	    smatcpp.append(", "+ rowData.elementAt(7)
                    + "*" + rowData.elementAt(8));
	 if (!rowData.elementAt(9).toString().equals("")){
	    smatcpp.append(", "+rowData.elementAt(9)+"*"
                    +rowData.elementAt(10));
	 }
      }
      smatcpp.append(" );\n");
     }
    }
    StringBuffer cmatcpp = new StringBuffer("\n// Materials from Combination\n\n");
    data = mcDataModel.getDataVector();
    Vector elements = new Vector();

    for (int i=0; i<data.size(); i++){
      rowData = (Vector)data.elementAt(i);
      if(!rowData.elementAt(0).equals("Used")){
//       System.out.println("Next");
      }else if(rowData.elementAt(0).equals("Used")){
       RatioItem rate = (RatioItem)rowData.elementAt(2);
       String name = (String)rowData.elementAt(1);
       cmatcpp.append("G4Material* " + name + " = new G4Material(\"" + name 
		 + "\",  " + rowData.elementAt(3)
		 +"*"+rowData.elementAt(4)+", "+ rate.getLength() + ", " 
                 + rowData.elementAt(5));
       if (!rowData.elementAt(6).toString().equals("")){
	 cmatcpp.append(", "+rowData.elementAt(6)+"*"
                      +rowData.elementAt(7));
	 if (!rowData.elementAt(8).toString().equals("")){
	  cmatcpp.append(", "+rowData.elementAt(8)+"*"+rowData.elementAt(9));
	 }
       }
       cmatcpp.append(" );\n");
       cmatcpp.append(rate.getCPP(name));
       for (int j=0; j<rate.elems.length; j++){
	if (!elements.contains(rate.elems[j])) elements.addElement(rate.elems[j]);
      }
     }
    }
    StringBuffer elecpp = new StringBuffer("// Elements\n");
    for (int i=0; i<elements.size(); i++){
      elecpp.append(((ElementItem)elements.elementAt(i)).getCPP());
    }
    elecpp.append(cmatcpp.toString());
    elecpp.append(smatcpp.toString());
    return elecpp.toString();
  }
}









