/*
 * Copyright (c) 1997 John Jensen. All rights reserved.
 *
 * This software is FREE FOR COMMERCIAL AND NON-COMMERCIAL USE,
 * provided the following condition is met.
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted,
 * provided that any copy or derivative of this software or documentation
 * retaining the name "John Jensen" also retains this condition and the
 * following disclaimer.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * CopyrightVersion 1.0
 */

package TextEdit;

import java.awt.*;
import java.awt.event.*;

// This extension of Dialog implements modal dialogs with a callback interface
// to work around the "modal dialog problem" on Windows.  The "NoteListener"
// interface is used to return an action code to the parent when "continue" is
// pressed.

class Note extends Dialog implements WindowListener, ActionListener
{
	int noteAction;
	NoteListener noteListener;
	Button accepted;
	Button whatever;

	// Constructor takes a string and an action.  Note that action == 0 implies
	// that there is only an "ok" button, while action != 0 implies a "Continue" and
	// a "Cancel" button.

	public Note(Frame parent, String s, int action)
	{
		super(parent, "Note", true);

		noteListener = (NoteListener)parent;
		noteAction = action;

		Panel center = new Panel();
		center.add(new Label(s));
		add("Center", center);

		Panel south = new Panel();

		if (action != 0)
		{
			accepted = new Button("Continue");
			accepted.addActionListener(this);
			south.add(accepted);
			whatever = new Button("Cancel");
			whatever.addActionListener(this);
			south.add(whatever);
		}
		else
		{
			whatever = new Button("Ok");
			whatever.addActionListener(this);
			south.add(whatever);
		}

		add("South", south);

		setSize(300, 100);

		addWindowListener(this);
	}

	// if we get a "Continue" button, call the parent's noteListener

    public void actionPerformed(ActionEvent event)
	{
		Object source;

		source = event.getSource();

        if (source == accepted)
		{
			setVisible(false);
			setModal(false);
			noteListener.noteListen(noteAction);
		}

        if ((source == accepted) || (source == whatever))
			dispose();
	}

	// add the 1.1 WindowListener stuff

    public void windowDeiconified(WindowEvent event) {}
    public void windowIconified(WindowEvent event) {}
    public void windowActivated(WindowEvent event) {}
    public void windowDeactivated(WindowEvent event) {}
    public void windowOpened(WindowEvent event) {}
    public void windowClosed(WindowEvent event) {}

    public void windowClosing(WindowEvent event) {
			dispose();
    }

}

