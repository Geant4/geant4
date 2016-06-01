/*GGE (Geant4 Geometry Editor)
 *Tetsuyay Yamada
 */

import java.io.*;
import java.lang.Runtime.*;
import java.awt.*;
import java.awt.event.*;

import com.sun.java.swing.*;
import com.sun.java.swing.event.*;

class CSGDialog extends SolidDialog implements ActionListener{
  private JTextField paramField[];
  private JComboBox paramUnit[];
  private CSGItem solid;
  private JButton viewButton, okButton;
  DAWNSourceFrame dawnSFrame;
  private String csgTitle;
  private Process p;

  CSGDialog(SolidCellEditor parent, CSGItem solid){
     super(parent.parent, solid.toString());
     this.solid = solid;


     JLabel lbl;
     dawnSFrame = new DAWNSourceFrame();
 
     GridBagConstraints gbc = new GridBagConstraints();
     GridBagLayout gbl = new GridBagLayout();
     getContentPane().setLayout(gbl);

     GridBagLayout csgGBL = new GridBagLayout();
     GridBagConstraints csgGBC = new GridBagConstraints();
                 JPanel csgPanel = new JPanel(csgGBL);

     gbc.gridx = 0; csgGBC.gridy = 0;
     gbc.gridwidth = GridBagConstraints.REMAINDER; gbc.gridheight = 1;
     gbc.weightx = 1; gbc.weighty = 1;
     gbc.fill = GridBagConstraints.HORIZONTAL;
     gbc.anchor = GridBagConstraints.CENTER;
      
       JPanel csgButtonPanel = new JPanel();
              csgButtonPanel.setLayout(new GridLayout(1,2));
              csgButtonPanel.add(okButton = new JButton(" OK "));
              csgButtonPanel.add(viewButton = new JButton(" VIEW "));

     gbl.setConstraints(csgButtonPanel, gbc);
     getContentPane().add(csgButtonPanel);

     int paraCount = solid.getParaCount();
     paramField = new JTextField[paraCount];
     paramUnit = new JComboBox[paraCount];
    for (int i=0; i<paraCount; i++){
      csgGBC.gridx = 0; csgGBC.gridy = i+1;
        csgGBC.gridwidth = 1; csgGBC.gridheight = 1;
        csgGBC.weightx = 1; csgGBC.weighty = 1;
        csgGBC.fill = GridBagConstraints.NONE;
        csgGBC.anchor = GridBagConstraints.EAST;
        lbl = new JLabel(solid.getParaName(i));
        csgGBL.setConstraints(lbl, csgGBC);
        csgPanel.add(lbl);
      csgGBC.gridx = 1;
        csgGBC.weightx = 1; csgGBC.weighty = 1;
        csgGBC.fill = GridBagConstraints.HORIZONTAL;
        paramField[i] = new JTextField(8);
        csgGBL.setConstraints(paramField[i], csgGBC);
        csgPanel.add(paramField[i]);
      csgGBC.gridx = 2;
        csgGBC.weightx = 1; csgGBC.weighty = 1;
        csgGBC.fill = GridBagConstraints.NONE;
        csgGBC.anchor = GridBagConstraints.WEST;
        if (solid.isAngle(i)) paramUnit[i] = new AngUnitCombo();
          else  paramUnit[i] = new LenUnitCombo();
        csgGBL.setConstraints(paramUnit[i], csgGBC);
        csgPanel.add(paramUnit[i]);
    }

    gbc.gridx = 0; gbc.gridy = 1;
         gbc.gridwidth = GridBagConstraints.REMAINDER; gbc.gridheight = 1;
         gbc.weightx = 0; gbc.weighty = 0;
         gbc.fill = GridBagConstraints.HORIZONTAL;
         gbc.anchor = GridBagConstraints.CENTER;
         gbl.setConstraints(csgPanel, gbc);
         getContentPane().add(csgPanel);

    okButton.addActionListener(parent);
    viewButton.addActionListener(this);
    setValues();
    pack();
    setVisible(true);
    setResizable(false);
  }
  void editStop(){
    dispose();
  }

  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();
    if(o == viewButton ){ 
      getValues();
      getDAWN();
    }
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
  private void setValues(){
    if (solid.isEmpty) return;
    for (int i=0; i<solid.getParaCount(); i++){
      paramField[i].setText(Double.toString(solid.getValues(i)));
      paramUnit[i].setSelectedItem(solid.getUnits(i));
    }
  }
  SolidItem getValues(){
    double tmp;
    for (int i=0; i<solid.getParaCount(); i++){
      try{
	tmp = new Double(paramField[i].getText()).doubleValue();
      }catch(NumberFormatException e){
	tmp = 0.0;
      }
      solid.setValuesUnits(i, tmp, (String)paramUnit[i].getSelectedItem());
    }
//    dispose();
    return solid;
  }
}




