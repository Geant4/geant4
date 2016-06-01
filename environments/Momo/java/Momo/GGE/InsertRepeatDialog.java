
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class InsertRepeatDialog extends JDialog implements ActionListener{
   private JButton b1, b2;
   private JLabel label1, label2;
   private RepeatedVolume repeatFrame;
   public InsertRepeatDialog(RepeatedVolume repeatFrame){
      super(repeatFrame,"Insert Repeated Volume",false);
      this.repeatFrame = repeatFrame;
      JPanel createRepeatPanel = new JPanel();
       createRepeatPanel.setLayout(new BorderLayout());
       b1 = new JButton("OK");
       b2 = new JButton("CANCEL");
       b1.addActionListener(this);
       b2.addActionListener(this);
       JPanel labelPanel = new JPanel();
        labelPanel.setLayout(new BorderLayout());
        label1 = new JLabel("   Click the pName after which ");
        label1.setFont(new Font("Serif",Font.BOLD, 12));
        label2 = new JLabel("   you are going to insert");
        label2.setFont(new Font("Serif",Font.BOLD, 13));
        labelPanel.add(label1, BorderLayout.NORTH);
        labelPanel.add(label2, BorderLayout.CENTER);
       createRepeatPanel.add(labelPanel, BorderLayout.NORTH);
       createRepeatPanel.add(b1, BorderLayout.CENTER);
       createRepeatPanel.add(b2, BorderLayout.EAST);
      getContentPane().add(createRepeatPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 100);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        repeatFrame.insertRepeat();
//        this.setVisible(false); 
      }else if(e.getSource() == b2){
        this.setVisible(false);       
        repeatFrame.insertRepeatCloseAct();
      }
   }
 
}









