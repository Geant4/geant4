/*
 *
 */

import java.awt.*;
import java.awt.event.*;

import com.sun.java.swing.*;
import com.sun.java.swing.border.*;
import com.sun.java.swing.event.*;
import com.sun.java.swing.preview.*;
import com.sun.java.swing.preview.JColorChooser;

class VisAtbDialog extends JDialog implements ActionListener {
  boolean isCanceled;
  VisMethodPanel atb[];
  private JButton ok, cancel;
  JTextField nameField;
  JColorChooser colorChooser;
  private VolumesFrame parent;
  VisAtbDialog(VolumesFrame parent) {
    super(parent, "VisAttribute", true);
    this.parent = parent;
    atb = new VisMethodPanel[5];
    getContentPane().setLayout(new BorderLayout(3,0));
    getContentPane().add("East", colorChooser = new JColorChooser());
    JPanel p1,p2,p3,c;
    p1 = new JPanel();
    p1.setLayout(new BorderLayout(1,8));
    p1.add("West",new JLabel("Name"));
    p1.add("Center", nameField = new JTextField());
    JComboBox lineStyle = new JComboBox();
    lineStyle.addItem("unbroken");
    lineStyle.addItem("dashed");
    lineStyle.addItem("dotted");
    c = new JPanel();
    c.setLayout(new BoxLayout(c, BoxLayout.Y_AXIS));
    c.setBorder(new LineBorder(Color.darkGray, 1));
    c.add(atb[0] = new VisMethodPanel(VisAtbItem.atbMethodName[0], new BooleanCombo()));
    c.add(atb[1] = new VisMethodPanel(VisAtbItem.atbMethodName[1], lineStyle));
    c.add(atb[2] = new VisMethodPanel(VisAtbItem.atbMethodName[2], new JTextField() ));
    c.add(atb[3] = new VisMethodPanel(VisAtbItem.atbMethodName[3], new BooleanCombo()));
    c.add(atb[4] = new VisMethodPanel(VisAtbItem.atbMethodName[4], new BooleanCombo()));
    p1.add("South", c);
    p2 = new JPanel(new BorderLayout());
    p2.add("North", p1);
    
    p3 = new JPanel(new BorderLayout(3,0));
    p3.add("Center", ok = new JButton("OK"));
    p3.add("East", cancel = new JButton("Cancel"));
    ok.addActionListener(this);
    cancel.addActionListener(this);
    p2.add("South", p3);
    getContentPane().add("Center", p2);
    pack();
    setResizable(false);
  }
  VisAtbItem getNewAtb(){
    nameField.setText("");
    setVisible(true);
    if (isCanceled) return null;
    return getValues();
  }
  VisAtbItem getChangedAtb(VisAtbItem visItem){
    setValues(visItem);
    setVisible(true);
    if (isCanceled) return null;
    return getValues();
  }
  private void setValues(VisAtbItem visItem){
    colorChooser.setColor(visItem.color);
    nameField.setText(visItem.name);
    for (int i=0; i<5; i++){
      atb[i].setAtbState(visItem.atbsIsValid[i]);
      atb[i].setValues(visItem.atbsValue[i]);
    }
  }
  private VisAtbItem getValues(){
    VisAtbItem newVisItem = new VisAtbItem();
    newVisItem.name = nameField.getText().trim();
    if (newVisItem.name.equals("")) return null;
    newVisItem.color = colorChooser.getColor();
    for (int i=0; i<5; i++){
      newVisItem.atbsIsValid[i] = atb[i].getAtbState();
      newVisItem.atbsValue[i] = atb[i].toString();
    }
    return newVisItem;
  }
  public void actionPerformed(ActionEvent ae){
    setVisible(false);
    parent.visPanel.visCloseAct();
    isCanceled = (ae.getSource() == cancel);
  }
}

class VisMethodPanel extends JPanel implements ItemListener {
  private JCheckBox checkbox;
  private JComponent compo;
  VisMethodPanel(String label, JComponent compo){
    this.compo = compo;
    setLayout(new BorderLayout());
    add("West", checkbox = new JCheckBox(label));
    add("Center", compo);
    compo.setEnabled(false);
    checkbox.addItemListener(this);
  }
  void setValues(String label){
    if (compo instanceof JTextField){
      ((JTextField)compo).setText(label);
      return;
    }
    if (compo instanceof JComboBox){
      ((JComboBox)compo).setSelectedItem(label);
      return;
    }
  }
  void setAtbState(boolean state){
    checkbox.setSelected(state);
  }
  boolean getAtbState(){
    return checkbox.isSelected();
  }
  public String toString(){
    if (!compo.isEnabled()) return null;
    if (compo instanceof JTextField){
      return ((JTextField)compo).getText();
    }
    if (compo instanceof JComboBox){
      return (String)((JComboBox)compo).getSelectedItem();
    }
    return null;
  }
  public void itemStateChanged(ItemEvent i){
    compo.setEnabled(i.getStateChange() == ItemEvent.SELECTED);
  }
}

class BooleanCombo extends JComboBox {
  BooleanCombo(){
    addItem("true");
    addItem("false");
  }
}







