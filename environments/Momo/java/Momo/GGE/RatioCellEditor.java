import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import com.sun.java.swing.table.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.io.Serializable;

public class RatioCellEditor extends Object implements TableCellEditor, ActionListener, java.io.Serializable {
  private Component compo;
  private Vector editorListener;
  private Object value;
  private RatioDialog ratioDialog;
  Frame parent;
  RatioCellEditor(Frame parent){
    this.parent = parent;
    editorListener = new Vector();
    compo = new JLabel();
  }
  public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int rowIndex, int columnIndex){
    this.value = value;
    return compo;
  }
  public Object getCellEditorValue(){
    ratioDialog.getValues();
    return value;
  } 
  public boolean isCellEditable(EventObject ev){
    return true;
  } 
  public boolean shouldSelectCell(EventObject ev){
    ratioDialog = new RatioDialog(this, (RatioItem)value );
    return true;
  } 
  public boolean stopCellEditing(){
    ratioDialog.editStop();
    castStop();
    return true;
  } 
  public void cancelCellEditing(){
    ratioDialog.editStop();
    castCancel();
  }
  private void castStop(){
    for (int i=0; i<editorListener.size(); i++){
      CellEditorListener cel = (CellEditorListener)editorListener.elementAt(i);
      cel.editingStopped(new ChangeEvent(this));
    }
  }
  private void castCancel(){
    for (int i=0; i<editorListener.size(); i++){
      CellEditorListener cel = (CellEditorListener)editorListener.elementAt(i);
      cel.editingCanceled(new ChangeEvent(this));
    }
  }
  public void addCellEditorListener(CellEditorListener ce){
    editorListener.addElement(ce);
  }
  public void removeCellEditorListener(CellEditorListener ce){
    editorListener.removeElement(ce);
  } 
  public void actionPerformed(ActionEvent ae){
    stopCellEditing();
  }
}
