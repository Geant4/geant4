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

class TextCursor extends Thread
{
	private		TextCanvas	canvas;
	private		boolean		flash;
	private		boolean		undraw;		// true when the next draw erases
	private		Rectangle	r;

	public TextCursor(TextCanvas c)
	{
		canvas = c;
		r = new Rectangle();
	}
	
	public void run()
	{
		while (true)
		{
			sync_draw();
			try { sleep(400); }
			catch(InterruptedException e) {}
		}
	}

	private synchronized void sync_draw()
	{
		if (flash)
			draw_or_undraw();
	}

	public boolean getUndraw()
	{
		return undraw;
	}

	public synchronized void pause_cursor(boolean draw)
	{
//		if (draw && undraw)
		if (undraw)
			draw_or_undraw();
		flash = false;
	}
	
	public synchronized void release_cursor(boolean draw)
	{
		if (!undraw)
			draw_or_undraw();
		flash = true;
	}

	private boolean draw_or_undraw()
	{
		Graphics g = canvas.getGraphics();

		if (g != null)
		{
			if (!undraw)
				canvas.getCursorPos(g,r);
			g.setXORMode(Color.white);
			g.drawLine(r.x, r.y, r.x, r.y+r.height);
			undraw = !undraw;
			g.dispose();
		}
		return undraw;
	}

//	public static void dumpStack() {new Exception("Stack trace").printStackTrace();}

}
