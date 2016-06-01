/*
 *
 */

import com.sun.java.swing.*;
import com.sun.java.swing.table.*;

import java.util.*;
import java.awt.*;

class DoubleCellEditor extends DefaultCellEditor {
  JFrame frame;
  DoubleCellEditor(){
    super(new JTextField());
    setClickCountToStart(1);
  }
  public Object getCellEditorValue(){
    Double value;
    String valueStr = (String)super.getCellEditorValue();
    if (valueStr.equals("")) return "";
    try{
      value = new Double(valueStr);
    }catch(NumberFormatException e){
      frame = new JFrame();
      JOptionPane.showMessageDialog(frame,"Not a Number", "Warning Dialog",
                                            JOptionPane.WARNING_MESSAGE,null);
      return "";
//      value = new Double(0.0);
    }
    return value;
  } 
}



