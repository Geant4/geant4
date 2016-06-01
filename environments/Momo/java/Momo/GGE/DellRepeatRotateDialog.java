
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class DellRepeatRotateDialog extends JDialog implements ActionListener{
   private JButton b1, b2;
   private JLabel label1;
   private RepeatedRotateVolume repeatFrame;
   public DellRepeatRotateDialog(RepeatedRotateVolume repeatFrame){
      super(repeatFrame,"Delete Repeated Volume",false);
      this.repeatFrame = repeatFrame;
      JPanel createRepeatPanel = new JPanel();
       createRepeatPanel.setLayout(new BorderLayout());
       b1 = new JButton("OK");
       b2 = new JButton("CANCEL");
       b1.addActionListener(this);
       b2.addActionListener(this);

        label1 = new JLabel("   Click the pName to delete");
        label1.setFont(new Font("Serif",Font.BOLD, 13));

       createRepeatPanel.add(label1, BorderLayout.NORTH);
       createRepeatPanel.add(b1, BorderLayout.CENTER);
       createRepeatPanel.add(b2, BorderLayout.EAST);
      getContentPane().add(createRepeatPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 80);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        repeatFrame.dellRepeat();
//        repeatFrame.dellRepeatCloseAct();
//        this.setVisible(false); 
      }else if(e.getSource() == b2){
        this.setVisible(false);       
        repeatFrame.dellRepeatCloseAct();
      }
   }
 
}








