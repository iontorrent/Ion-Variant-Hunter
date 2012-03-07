/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2flowgram.flowalign;

import java.io.PrintStream;
import org.iontorrent.sam2flowgram.util.SamToFlowgramAlignUtil;

/**
 * Represents an base sequence in flow space.
 *
 * @author nils.homer@lifetech.com
 */
public class FlowSeq {
    /**
     * The gap alignment character.
     */
    public static char GAP = '-';

    /**
     * The flow sequence in integers (100x SFF format).
     */
    public int flow[] = null;

    /**
     * The flow sequence length.
     */
    public int length = 0;

    /**
     * The amount of memory allocated for the flow sequence.
     */
    private int mem = 0;

    /** 
     * The first non-empty flow.
     */
    public int nonEmptyFlowFirst;

    /**
     * The last non-empty flow.
     */
    public int nonEmptyFlowLast;

    /**
     * Creates a new flow sequence.
     * @param seq the sequence in integer format.
     * @param flowOrder the flow order in integer format.
     */
    public FlowSeq(byte seq[], byte flowOrder[])
        throws Exception
    {
        this(seq, flowOrder, 0);
    }
    
    /**
     * Creates a new flow sequence.
     * @param seq the sequence in integer format.
     * @param flowOrder the flow order in integer format.
     * @param startFlowIndex the zero-based index in which to start in the flow order.
     */
    public FlowSeq(byte seq[], byte flowOrder[], int startFlowIndex)
        throws Exception
    {
        this(seq, null, flowOrder, startFlowIndex);
    }
    
    /**
     * Creates a new flow sequence.
     * @param seq the sequence in integer format.
     * @param signals the flow signals for the read (100x), null if not present
     * @param flowOrder the flow order in integer format.
     */
    public FlowSeq(byte seq[], int signals[], byte flowOrder[])
        throws Exception
    {
        this(seq, signals, flowOrder, 0);
    }
    
    /**
     * Creates a new flow sequence.
     * @param seq the sequence in integer format.
     * @param signals the flow signals for the read (100x), null if not present
     * @param flowOrder the flow order in integer format.
     * @param startFlowIndex the zero-based index in which to start in the flow order.
     */
    public FlowSeq(byte seq[], int signals[], byte flowOrder[], int startFlowIndex)
        throws Exception
    {
        int i, k, l, nextI;

        if(null != signals) {
            this.mem = signals.length;
        }
        else {
            // ? is this a good approximation of the flow length ?
            this.mem = seq.length * flowOrder.length; 
            // get rid of Ns
            for(i=0;i<seq.length;i++) {
                if(3 < seq[i]) {
                    seq[i] = 0; // Ns to As
                }
            }
        }
        this.flow = new int[this.mem];

        i = 0;
        k = startFlowIndex;
        while(i < seq.length) {
            int before = i;

            // move beyond the initial gaps
            while(i < seq.length && GAP == seq[i]) {
                i++;
            }
            if(seq.length <= i) break;

            // skip over empty flow
            while(flowOrder[k] != seq[i] && seq[i] <= 3) {
                this.flow[this.length] = 0;
                this.length++;
                k = (k+1) % flowOrder.length;
            }

            // get the number of bases in current flow
            nextI = i+1;
            l = 1;
            while(nextI < seq.length && (flowOrder[k] == seq[nextI] || 3 < seq[nextI])) {
                if(flowOrder[k] == seq[nextI] || 3 < seq[i]) {
                    l++;
                }
                nextI++;
            }
            this.flow[this.length] = l * 100;
            this.length++;
            k = (k+1) % flowOrder.length;
            i = nextI;
            if(i <= before) {
                throw new Exception("i <= before ["+i+"<="+before+"]");
            }
        }

        this.nonEmptyFlowFirst = 0;
        while(this.nonEmptyFlowFirst < this.flow.length && 0 == this.flow[this.nonEmptyFlowFirst]) {
            this.nonEmptyFlowFirst++;
        }
        
        this.nonEmptyFlowLast = this.length-1;
        while(0 < this.nonEmptyFlowLast && 0 == this.flow[this.nonEmptyFlowLast]) {
            this.nonEmptyFlowLast--;
        }

        // copy over flow signals
        if(null != signals) {
            for(i=startFlowIndex, k=0;i<signals.length;i++,k++) {
                // the read sequence has the priority
                if(SamToFlowgramAlignUtil.getBaseCallFromFlowSignal(this.flow[k]) == SamToFlowgramAlignUtil.getBaseCallFromFlowSignal(signals[k])) {
                    this.flow[k] = signals[k];
                }
                else {
                    this.flow[k] += (signals[k] - (SamToFlowgramAlignUtil.getBaseCallFromFlowSignal(signals[k]) * 100));
                    if(this.flow[k] < 0) {
                        this.flow[k] = 0;
                    }
                }
            }
        }
    }

    /**
     * Debuggin print function.
     * @param out the output stream.
     */
    public void print(PrintStream out)
    {
        int i;

        for(i=0;i<this.length;i++) {
            if(0 < i) {
                out.print(",");
            }
            out.print(this.flow[i]);
        }
        out.println("");
    }
} 
