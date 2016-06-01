import com.sun.java.swing.*;

class LenUnitCombo extends JComboBox {
  LenUnitCombo(){
    addItem("mm");
    addItem("cm");
    addItem("m");
    addItem("km");
    addItem("micrometer");
    addItem("nanoometer");
    addItem("fermi");
    setMaximumRowCount(3);
  }
}
