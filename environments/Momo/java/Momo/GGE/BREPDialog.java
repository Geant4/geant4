/*  GGE (Geant4 Geometry Editor)
 *  Tetsuyay Yamada
 *  9/9
 */
// title no hyouji no henkou 1998 10 30



import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;
import java.lang.Runtime.*;

import com.sun.java.swing.*;
import com.sun.java.swing.table.*;
import com.sun.java.swing.event.*;

class BREPDialog extends SolidDialog implements ActionListener {
  private BREPItem solid;
  private JTextField phiField[], sidesField;
  private JComboBox sphiAngUnit, dphiAngUnit;
  private JButton addButton, delButton, view, okButton;
  private JFrame frame;
  private JTable zTable;
  private String tableHead[] = {"Z","RMax","RMin","Unit"};
  private DefaultTableModel dataModel;
  DAWNSourceFrame dawnSFrame;

  BREPDialog(SolidCellEditor parent, BREPItem solid){
    super(parent.parent, solid.toString());
    this.solid = solid;
    dawnSFrame = new DAWNSourceFrame();
    getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));    

    JPanel brepButtonPanel = new JPanel();
           brepButtonPanel.setLayout(new GridLayout(1,2));
           brepButtonPanel.add(okButton = new JButton(" OK "));
           brepButtonPanel.add(view = new JButton(" VIEW "));

    getContentPane().add(brepButtonPanel);

    GridBagLayout gbl = new GridBagLayout();
    GridBagConstraints gbc = new GridBagConstraints();
    JPanel panel = new JPanel(gbl);
    JLabel phiLabel[] = new JLabel[2];
    phiField = new JTextField[2];
    sphiAngUnit = new AngUnitCombo();
    dphiAngUnit = new AngUnitCombo();
    for (int i=0; i<solid.phi.length; i++){
      gbc.gridx = 0; gbc.gridy = i;
      gbc.gridwidth = 1; gbc.gridheight = 1;
      gbc.weightx = 0.0; gbc.weighty = 0.0;
      gbc.fill = GridBagConstraints.NONE;
      phiLabel[i] = new JLabel(solid.getParaName(i));
      gbl.setConstraints(phiLabel[i], gbc); panel.add(phiLabel[i]);
      gbc.fill = GridBagConstraints.HORIZONTAL;
      gbc.weightx = 1.0;
      gbc.gridx = 1;
      phiField[i] = new JTextField(8);
      gbl.setConstraints(phiField[i], gbc); panel.add(phiField[i]);
    }
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridx = 2; gbc.gridy = 0;
    gbc.gridwidth = 1; gbc.gridheight = 1;
    gbc.weightx = 0.0;
    gbl.setConstraints(sphiAngUnit, gbc); panel.add(sphiAngUnit);

    gbc.fill = GridBagConstraints.NONE;
    gbc.gridx = 2; gbc.gridy = 1;
    gbc.gridwidth = 1; gbc.gridheight = 1;
    gbc.weightx = 0.0;
    gbl.setConstraints(dphiAngUnit, gbc); panel.add(dphiAngUnit);

    if (solid.getSolidType() == BREPItem.PGON){
      JLabel sidesLabel = new JLabel("Sides");
      sidesField = new JTextField();
      gbc.gridx = 0; gbc.gridy = 2;
      gbc.gridwidth = 1; gbc.gridheight = 1;
      gbc.weightx = 0.0;
      gbc.fill = GridBagConstraints.NONE;
      gbl.setConstraints(sidesLabel, gbc); panel.add(sidesLabel);
      gbc.gridx = 1;
      gbc.fill = GridBagConstraints.HORIZONTAL;
      gbc.weightx = 1.0;
      gbl.setConstraints(sidesField, gbc); panel.add(sidesField);
    }
    zTable = new JTable(dataModel = new DefaultTableModel(tableHead ,0));
    zTable.setAutoCreateColumnsFromModel(false);
    zTable.setColumnSelectionAllowed(false);
    zTable.setRowSelectionAllowed(true);
    zTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    DefaultCellEditor zCellEditor = new IntCellEditor();
    DefaultCellEditor textCellEditor = new DoubleCellEditor();
    LenUnitCombo lenUnitCombo = new LenUnitCombo();
    textCellEditor.setClickCountToStart(2);
    zCellEditor.setClickCountToStart(2);
    TableColumn columns[] = new TableColumn[4];
    for(int i=0;i<4;i++){    
      columns[i] = zTable.getColumn(tableHead[i]);
    }
    DefaultTableCellRenderer zCellRenderer = new DefaultTableCellRenderer();
    zCellRenderer.setBackground(Color.pink);
    columns[0].setCellEditor(textCellEditor);
    columns[0].setCellRenderer(zCellRenderer);    
    columns[1].setCellEditor(textCellEditor);
    columns[2].setCellEditor(textCellEditor);
    columns[3].setCellEditor(new DefaultCellEditor(lenUnitCombo));
    JScrollPane scrollpane = new JScrollPane(zTable);
    scrollpane.setPreferredSize(new Dimension(180, 80));
    scrollpane.setMinimumSize(new Dimension(40, 20));
    gbc.gridx = 0; gbc.gridy = 3;
    gbc.gridwidth = 3; gbc.gridheight = 1;
    gbc.weightx = 1.0; gbc.weighty = 0.0;
    gbc.fill = GridBagConstraints.BOTH;
    getContentPane().add(panel);

    JPanel panel2 = new JPanel();
     panel2.setLayout(new GridLayout(1,2));
     panel2.add(addButton = new JButton("Add a Section"));
    addButton.addActionListener(this);
     panel2.add(delButton = new JButton("Del a Section"));
    delButton.addActionListener(this);
    gbl.setConstraints(panel2, gbc); panel.add(panel2);
    getContentPane().add(panel);

    gbc.gridx = 0; gbc.gridy = 4;
    gbc.gridwidth = 3; gbc.gridheight = 1;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    gbc.fill = GridBagConstraints.BOTH;
    gbl.setConstraints(scrollpane, gbc); panel.add(scrollpane);
//    getContentPane().add(scrollpane);
    view.addActionListener(this);
    okButton.addActionListener(parent);
    setValues();
    pack();
    setVisible(true);
  }

  private void setValues(){
    if (solid.isEmpty){
      return;
    }else{
      for(int i=0; i<solid.phi.length; i++){
          phiField[i].setText(Double.toString(solid.phi[i]));
      }

      if (solid.getSolidType() == BREPItem.PGON){
          sidesField.setText(Integer.toString(solid.nSides));
      }

      for(int j=0; j<solid.zVector.size(); j++){
         Vector br = (Vector)solid.zVector.elementAt(j);
         dataModel.addRow(br);
      }

    }
    sphiAngUnit.setSelectedItem(solid.sAngUnit);
    dphiAngUnit.setSelectedItem(solid.dAngUnit);
  }

  void editStop(){
    setVisible(false);
  }
  SolidItem getValues(){
    for(int i=0; i<2; i++){
      try{
	solid.phi[i] = new Double(phiField[i].getText()).doubleValue();
      }catch(NumberFormatException e){
	solid.phi[i] = 0.0;
      }
    }
    if (solid.getSolidType() == BREPItem.PGON){
      try{
	solid.nSides = Integer.parseInt(sidesField.getText());
      }catch(NumberFormatException e){
	solid.nSides = 0;
      }
    }
    solid.zVector = dataModel.getDataVector();

    solid.sAngUnit = (String)sphiAngUnit.getSelectedItem();
    solid.dAngUnit = (String)dphiAngUnit.getSelectedItem();
    solid.isEmpty = false;
    dispose();
    return solid;
  }

  private void getDAWN(){
    dawnSFrame.setVisible(true);
    dawnSFrame.dawnEditor.setText(solid.getPrim());
    try{
      FileOutputStream fileout = new FileOutputStream("GGESolid.prim");
      PrintWriter outf = new PrintWriter(new DataOutputStream(fileout));
      outf.println(solid.getPrim());
      outf.close();
    }catch(IOException e){
      System.out.println(e.toString());
    }
    try{
      String cmd="dawn -v GGESolid.prim";
      Runtime r = Runtime.getRuntime();
      Process p = r.exec(cmd);
    }catch(IOException e){
      System.out.println(e.toString());
    }
  }

  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();

    if (o == addButton){
       Object tmp[] = new Object[4];
       tmp[0] = "";
       tmp[1] = "";
       tmp[2] = "";
       tmp[3] = "mm";
       dataModel.addRow(tmp);
    }
    if(o == view){
       getValues();
       getDAWN();

    }
    if (o == delButton) {
       if(zTable.getSelectedRow()==-1){
          frame = new JFrame();
          JOptionPane.showMessageDialog(frame, "Choose the  Z  ","Warning Dialog", JOptionPane.WARNING_MESSAGE);
       }else{
          dataModel.removeRow(zTable.getSelectedRow());
          return;
       }
    }
  }
}






