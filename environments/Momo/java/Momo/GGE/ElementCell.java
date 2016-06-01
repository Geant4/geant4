// Table of Elements

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;

class ElementCell extends JPanel implements ActionListener {
  static int selectedElems = 0;
  ElementItem ei;
  private boolean isSelected;
  private JButton e;
  private Color bg;
  ElementCell(ElementItem ei){
    this.ei = ei;
    isSelected = false;
    setLayout(new GridLayout(2,1));
    if (ei.isMetal){
      bg =  new Color(255, 175, 175);
    }else{
      bg  = new Color(175, 175, 255);
    }
    setBackground(bg);
    e = new JButton(ei.symbol);
    e.setMargin(new Insets(1,3,1,3));
    e.setToolTipText(ei.name);
    e.addActionListener(this);
    add(e);
    add(new JLabel(Integer.toString(ei.atomNum), JLabel.CENTER));
  }
  public void actionPerformed(ActionEvent ev){
    Object source = ev.getSource();
    if (isSelected){
      selectedElems--;
      isSelected = false;
      setBackground(bg);
    }else{
      if (selectedElems >= 8) return;
      selectedElems++;
      isSelected = true;
      setBackground(Color.lightGray);
    }
    repaint();
  }
  void clear(){
    selectedElems = 0;
    isSelected = false;
    setBackground(bg);
    repaint();
  }
  boolean isSelected(){
    return isSelected;
  }
}
