/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2flowgram.flowalign;

import java.io.PrintStream;
import org.iontorrent.sam2flowgram.util.SamToFlowgramAlignUtil;

/**
 * Stores the flow order for a give read or query.
 *
 * @author nils.homer@lifetech.com
 */
public class FlowOrder
{
    /**
     * The flow order in integer values.
     */
    public byte flowOrder[] = null;

    /**
     * The length or period of the flow order.
     */
    public int length = 0;

    /**
     * Forward jump table.  Indicates the number bases after a given base the 
     * same base is observed (cyclically).
     */
    public int jumpFwd[] = null;

    /**
     * Reverse jump table.  Indicates the number bases before a given base the 
     * same base is observed (cyclically).
     */
    public int jumpRev[] = null;

    /**
     * The key sequence, if any.
     */
    public byte keySequence[] = null;

    /**
     * Creates the forward and reverse jump table from the flow order.
     */
    private void createJumpTables()
    {
        int i, j, k;
        this.jumpFwd = new int[this.length];
        this.jumpRev = new int[this.length];

        for(i=0;i<this.length;i++) {
            k = 1;
            j = (i + 1) % this.length;
            while(this.flowOrder[i] != this.flowOrder[j] && i != j) {
                j = (j + 1) % this.length;
                k++;
            }
            this.jumpFwd[i] = this.jumpRev[j] = k;
        }
    }
    
    /**
     * Creates a new flow order, removing homopolymers if desired.
     * @param flowOrder the flow order in integer format.
     * @param compress compress homopolymers to one base if true.
     */
    public FlowOrder(byte flowOrder[], boolean compress)
    {
        int i, l;

        if(compress) {
            // 1st pass, get the compressed length
            for(i=l=1;i<flowOrder.length;i++) {
                if(flowOrder[i-1] != flowOrder[i]) {
                    l++;
                }
            }
            // 2nd pass, get the compressed length
            this.flowOrder = new byte[l];
            this.flowOrder[0] = flowOrder[0];
            for(i=1,l=1;i<flowOrder.length;i++) {
                if(flowOrder[i-1] != flowOrder[i]) {
                    this.flowOrder[l] = flowOrder[i];
                    l++;
                }
            }
            this.length = this.flowOrder.length;
            this.createJumpTables();
        }
        else {
            this.flowOrder = new byte[flowOrder.length];
            this.length = this.flowOrder.length;
            System.arraycopy(flowOrder, 0, this.flowOrder, 0, this.length);
            this.createJumpTables();
        }
    }

    /**
     * Creates a new flow order.
     * @param flowOrder the flow order in integer format.
     */
    public FlowOrder(byte flowOrder[])
    {
        this.flowOrder = new byte[flowOrder.length];
        this.length = this.flowOrder.length;
        System.arraycopy(flowOrder, 0, this.flowOrder, 0, this.length);
        this.createJumpTables();
    }
    
    /**
     * Creates a new flow order.
     * @param flowOrder the flow order.
     * @param keySequence the key sequence.
     */
    public FlowOrder(String flowOrder, String keySequence)
    {
        int i, j;

        // flow order
        this.flowOrder = flowOrder.getBytes();
        this.length = this.flowOrder.length;
        SamToFlowgramAlignUtil.ntToInt(this.flowOrder);
        
        // key sequence
        this.keySequence = keySequence.getBytes();
        SamToFlowgramAlignUtil.ntToInt(this.keySequence);
    }

    /**
     * Creates a new flow order.
     * @param flowOrder the flow order.
     */
    public FlowOrder(String flowOrder)
    {
        this.flowOrder = flowOrder.getBytes();
        this.length = this.flowOrder.length;
        SamToFlowgramAlignUtil.ntToInt(this.flowOrder);
        this.createJumpTables();
    }

    /**
     * Creates a new flow order.
     * @param flowOrder the flow order.
     */
    public FlowOrder(FlowOrder flowOrder)
    {
        this.length = flowOrder.length;
        this.flowOrder = new byte[flowOrder.flowOrder.length];
        System.arraycopy(flowOrder.flowOrder, 0, this.flowOrder, 0, this.length);
        this.keySequence = new byte[flowOrder.keySequence.length];
        System.arraycopy(flowOrder.keySequence, 0, this.keySequence, 0, this.keySequence.length);
        this.createJumpTables();
    }

    /**
     * Cyclically shifts the flow order to match the starting flow index.
     * @param startFlowIndex the zero-based index in the flow order to which to shift.
     * @param newJumpTables true if we are to create new jump tables, false otherwise
     */
    public void rotate(int startFlowIndex, boolean newJumpTables)
    {
        int i, j;
        byte[] tmp = null;

        if(0 != startFlowIndex) {
            tmp = new byte[this.length];
            System.arraycopy(this.flowOrder, 0, tmp, 0, this.length);
            for(i=0;i<this.length;i++) {
                j = (i + startFlowIndex) % this.length;
                this.flowOrder[i] = tmp[j];
            }
            this.createJumpTables();
        }
    }

    /**
     * Debugging print function.
     * @param out the output stream.
     */
    public void print(PrintStream out)
    {
        int i;
        out.println("flow order:");
        for(i=0;i<this.flowOrder.length;i++) {
            if(0 < i) {
                out.print(",");
            }
            out.print(this.flowOrder[i]);
        }
        out.println("");
        
        out.println("jump fwd:");
        for(i=0;i<this.jumpFwd.length;i++) {
            if(0 < i) {
                out.print(",");
            }
            out.print(this.jumpFwd[i]);
        }
        out.println("");

        out.println("jump rev:");
        for(i=0;i<this.jumpRev.length;i++) {
            if(0 < i) {
                out.print(",");
            }
            out.print(this.jumpRev[i]);
        }
        out.println("");
    }
}
