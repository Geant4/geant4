/*
 *
 */

import java.awt.*;
import java.awt.event.*;

import com.sun.java.swing.*;
import com.sun.java.swing.event.*;

class RatioDialog extends JDialog {
  private JTextField paramField[];
  private JComboBox ratioType;
  private CSGItem solid;
  private RatioItem ratioItem;
  RatioDialog(RatioCellEditor parent, RatioItem ratioItem){
    super(parent.parent, "Mix Ratio");
    this.ratioItem = ratioItem;
    JLabel lbl;
    GridBagConstraints gbc = new GridBagConstraints();
    GridBagLayout gbl = new GridBagLayout();
    getContentPane().setLayout(gbl);

    gbc.gridx = 0; gbc.gridy = 0;
    gbc.gridwidth = 1; gbc.gridheight = 1;
    gbc.weightx = 0.0; gbc.weighty = 0.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.EAST;
    lbl = new JLabel("Ratio");
    gbl.setConstraints(lbl, gbc);
    getContentPane().add(lbl);
    gbc.gridx = 1;
    gbc.weightx = 1.0; gbc.weighty = 0.0;
    gbc.fill = GridBagConstraints.HORIZONTAL;
    ratioType = new JComboBox();
    ratioType.addItem("Number");
    ratioType.addItem("Fraction");
    gbl.setConstraints(ratioType, gbc);
    getContentPane().add(ratioType);

    int paraCount = ratioItem.getLength();
    paramField = new JTextField[paraCount];
    for (int i=0; i<paraCount; i++){
      gbc.gridx = 0; gbc.gridy = i+1;
      gbc.gridwidth = 1; gbc.gridheight = 1;
      gbc.weightx = 0.0; gbc.weighty = 0.0;
      gbc.fill = GridBagConstraints.NONE;
      gbc.anchor = GridBagConstraints.EAST;
      lbl = new JLabel(ratioItem.getEleName(i));
      gbl.setConstraints(lbl, gbc);
      getContentPane().add(lbl);
      gbc.gridx = 1;
      gbc.weightx = 1.0; gbc.weighty = 0.0;
      gbc.fill = GridBagConstraints.HORIZONTAL;
      paramField[i] = new JTextField(8);
      gbl.setConstraints(paramField[i], gbc);
      getContentPane().add(paramField[i]);
    }
    gbc.gridx = 0; gbc.gridy = paraCount+1;
    gbc.gridwidth = 2; gbc.gridheight = 1;
    gbc.weightx = 1.0; gbc.weighty = 0.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.CENTER;
    JButton okButton = new JButton(" OK ");
    gbl.setConstraints(okButton, gbc);
    getContentPane().add(okButton);
    okButton.addActionListener(parent);
    setValues();
    pack();
    setVisible(true);
  }
  void editStop(){
    setVisible(false);
  }
  private void setValues(){
    if (ratioItem.isEmpty) return;
    ratioType.setSelectedIndex(ratioItem.isFraction);
    for (int i=0; i<ratioItem.getLength(); i++){
      paramField[i].setText( ratioItem.getRatioStr(i) );
    }
  }
  RatioItem getValues(){
    ratioItem.isFraction = ratioType.getSelectedIndex();
    for (int i=0; i<ratioItem.getLength(); i++){
      try{
	ratioItem.ratio[i] = new Float(paramField[i].getText()).floatValue();
      }catch(NumberFormatException e){
	ratioItem.ratio[i] = 0f;
      }
    }
    ratioItem.isEmpty = false;
    dispose();
    return ratioItem;
  }
}
