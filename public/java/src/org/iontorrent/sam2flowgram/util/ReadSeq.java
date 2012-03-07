/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2flowgram.util;

import org.iontorrent.sam2flowgram.flowalign.*;
import java.io.PrintStream;
import java.util.*;
import net.sf.samtools.*;

/**
 * Stores information about the read sequence after initial
 * alignment, partition soft-clipping sequence.
 *
 * All fields relative to the forward genomic strand.
 *
 * @author nils.homer@lifetech.com
 */
public class ReadSeq {
    /**
     * The read bases free of soft-clipped bases.
     */
	public String readBases = null;
    
    /**
     * The read base qualities free of soft-clipped base qualities.
     */
	public String readQuals = null;

    /**
     * The soft-clipped bases at the start of the read.
     */
	public String softClipStartBases = null;

    /**
     * The soft-clipped base qualities at the start of the read.
     */
	public String softClipStartQuals = null;

    /**
     * The soft-clipped bases at the end of the read.
     */
	public String softClipEndBases = null;

    /**
     * The soft-clipped base qualities at the end of the read.
     */
	public String softClipEndQuals = null;

    /**
     * The number of hard-clipped bases at the start of the read.
     */
    public int hardClipStartBasesLength = 0;

    /**
     * The number of hard-clipped bases at the end of the read.
     */
    public int hardClipEndBasesLength = 0;

    /**
     * The read bases in integer format free of soft-clipped bases.
     */
	public byte readBytes[] = null;

    /**
     * The index into the flow order for the first non-empty read base.
     */
    private int flowOrderIndexStart = -1;

    /**
     * If the last base of the key sequence is mixed with the first base of the read sequence.
     * NB: getFlowOrderIndexStart must be set.
     */
    private int startSequenceOverlapAmount = -1;

    /**
     * Given an initial SAM record, partitions the read bases and read
     * base qualities according to soft-clipping.
     */
	public ReadSeq(SAMRecord record)
    {
        List<CigarElement> cigarElements = null;
		CigarElement eStart = null, eEnd = null;
		int size, length, startL, endL;

		startL = endL = 0;

		// read string
		this.readBases = record.getReadString();
		this.readQuals = record.getBaseQualityString();

        // ignore unmapped reads
        if(!record.getReadUnmappedFlag()) {
            // hard and soft clipping
            cigarElements = record.getCigar().getCigarElements();
            size = cigarElements.size();

            // start hard and soft clip
            eStart = cigarElements.get(0);
            // first hard clip
            if(1 < size && CigarOperator.H == eStart.getOperator()) {
                this.hardClipStartBasesLength = eStart.getLength();
                // move to the next cigar operator
                eStart = cigarElements.get(1);
            }
            // first soft clip
            if(CigarOperator.S == eStart.getOperator()) {
                length = eStart.getLength();
                softClipStartBases = this.readBases.substring(0, length);
                softClipStartQuals = this.readQuals.substring(0, length);
                this.readBases = this.readBases.substring(length, this.readBases.length());
                this.readQuals = this.readQuals.substring(length, this.readQuals.length());
                startL = length;
            }

            // end hard and soft clip
            eEnd = cigarElements.get(size-1);
            // last hard clip
            if(1 < size && CigarOperator.H == eEnd.getOperator()) {
                this.hardClipEndBasesLength = eEnd.getLength();
                // move to the previous cigar operator
                eEnd = cigarElements.get(size-2);
            }
            // last soft clip
            if(CigarOperator.S == eEnd.getOperator()) {
                length = eEnd.getLength();
                softClipEndBases = this.readBases.substring(this.readBases.length() - length);
                softClipEndQuals = this.readQuals.substring(this.readQuals.length() - length);
                this.readBases = this.readBases.substring(0, this.readBases.length() - length);
                this.readQuals = this.readQuals.substring(0, this.readQuals.length() - length);
                endL = length;
            }
        }

		// read bytes
		this.readBytes = new byte[record.getReadLength() - startL - endL]; // allocate
		System.arraycopy(record.getReadBases(), startL, this.readBytes, 0, record.getReadLength() - startL - endL); // copy
		SamToFlowgramAlignUtil.ntToInt(this.readBytes); // to integer
	}

    /**
     * Given the integer flow order, integer soft-clipped bases, and the first integer read base,
     * computes the flow order index (zero-based) where the first read base starts.  This assumes
     * the soft-clipped bases are in sequencing order.
     * @param flowOrder the integer flow order.
     * @param softClipBases the soft-clipped bases before the first base, null if none exist.
     * @param readBase the first sequenced readBase.
     * @return the flow index.
     */
    private int updateFlowInformationHelper(byte[] flowOrder, byte[] keySequence, int[] flowSignals, 
            int hardClipBasesLength, byte[] softClipBases, byte readBase)
    {
        int i, j, l, k;
        int lastBaseCall = 0;
        
        // Start at the first flow
        // i - array counter
        // j - flow order index 
        // k - tmp variable 
        // l - flow sequence index 
        j = k = l = 0;

        // Go through the key sequence
        if(null != keySequence) {
            i = 0;
            while(i<keySequence.length) {
                while(keySequence[i] != flowOrder[j]) {
                    j = (j + 1) % flowOrder.length;
                    l++;
                }
                k = 1;
                i++;
                while(i < keySequence.length && keySequence[i] == flowOrder[j]) {
                    i++;
                    k++;
                }
                lastBaseCall = k;
            }
        }

        // Go through the hard clipped bases
        if(0 < hardClipBasesLength) {
            i = hardClipBasesLength;
            while(0 < i) {
                k = SamToFlowgramAlignUtil.getBaseCallFromFlowSignal(flowSignals[l]); // TODO: this assumes that the base calls correspond to the implied base calls from the flow signals
                if(0 <= i - k) {
                    j = (j + 1) % flowOrder.length;
                    l++;
                }
                else { // mix between the last hard clip base and the first template base
                    lastBaseCall = k;
                }
                i -= k;
            }
        }
        
        if(null != softClipBases) { // Go through the soft-clipped bases
            for(i=0;i<softClipBases.length;i++) {
                while(softClipBases[i] != flowOrder[j]) {
                    j = (j + 1) % flowOrder.length;
                    l++;
                }
            }
            this.startSequenceOverlapAmount = 0; // default: no overlap
        }
        else if(readBase == flowOrder[j]) {
            this.startSequenceOverlapAmount = SamToFlowgramAlignUtil.getFlowSignalFromBaseCall(lastBaseCall);
        }
        else {
            this.startSequenceOverlapAmount = 0; // default: no overlap
        }

        // Move to the first base in the read
        while(readBase != flowOrder[j]) {
            j = (j + 1) % flowOrder.length;
            l++;
        }

        return l;
    }

    /**
     * Given the flow order and strand of the read, gives the zero-based index
     * into the flow order corresponding to the first non-soft-clipped base in the read
     * (sequencing order).  This will also consider the key seuqence
     * @param strand the original strand of the alignment.
     * @param flowOrder the flow order 
     * @return the flow index.
     */
    public int[] updateFlowInformation(boolean strand, FlowOrder flowOrder, int[] flowSignals)
        throws Exception
    {
        byte[] softClipBases = null;
        int hardClipBasesLength = 0;
        byte readBase;
        int i;
        int[] newFlowSignals = null;

        if(this.flowOrderIndexStart < 0) { // update
            // cycle through the key sequence
            if(strand) {
                // We need the hard clipped bases
                if(0 < this.hardClipEndBasesLength) {
                    hardClipBasesLength = this.hardClipEndBasesLength;
                }
                // We need to reverse the soft end clipped bases
                if(null != this.softClipEndBases) {
                    softClipBases = SamToFlowgramAlignUtil.basesToInt(this.softClipEndBases);
                    SamToFlowgramAlignUtil.reverseCompliment(softClipBases);
                }

                readBase = SamToFlowgramAlignUtil.NTINT2COMP[this.readBytes[this.readBytes.length-1]];
            }
            else {
                // We need the hard clipped bases
                if(0 < this.hardClipStartBasesLength) {
                    hardClipBasesLength = this.hardClipStartBasesLength;
                }
                // We need the soft clipped bases
                if(null != this.softClipStartBases) {
                    softClipBases = SamToFlowgramAlignUtil.basesToInt(this.softClipStartBases);
                }
                readBase = this.readBytes[0];
            }
            this.flowOrderIndexStart = updateFlowInformationHelper(flowOrder.flowOrder, flowOrder.keySequence, flowSignals, hardClipBasesLength, softClipBases, readBase);
            // update the flow signals
            newFlowSignals = new int[flowSignals.length - this.flowOrderIndexStart];
            for(i=0;i<flowSignals.length - this.flowOrderIndexStart;i++) {
                newFlowSignals[i] = flowSignals[i+this.flowOrderIndexStart];
            }
        }
        
        return newFlowSignals;
    }

    /**
     * @return the key sequence overlap
     */
    public int getStartSequenceOverlap()
        throws Exception
    {
        if(this.startSequenceOverlapAmount < 0) {
            throw new Exception("Key sequence overlap not set");
        }
        return this.startSequenceOverlapAmount;
    }

    /**
     * @return the flow index start, skipping over the key sequence and hard clip bases
     */
    public int getFlowOrderIndexStart()
    {
        return this.flowOrderIndexStart;
    }
    
    /**
     * @return the read bases in integer format.
     */
    public byte[] getReadBytes()
    {
        return this.readBytes;
    }

    /**
     * Debugging function.
     * @param out the output stream.
     */
	public void print(PrintStream out)
	{
		out.println(this.readBases);
		out.println(this.readQuals);
		out.println(this.softClipStartBases);
		out.println(this.softClipStartQuals);
		out.println(this.softClipEndBases);
		out.println(this.softClipEndQuals);
	}
}
