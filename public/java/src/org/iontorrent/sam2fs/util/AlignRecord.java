/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
//here be dragons
package org.iontorrent.sam2fs.util;

import java.io.IOException;
import java.lang.reflect.Array;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import org.iontorrent.sam2fs.flowspace.*;

import java.util.*;
import net.sf.samtools.*;
import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.DateParser;
import net.sf.samtools.util.RuntimeIOException;
import net.sf.samtools.util.StringUtil;

/**
 * Stores information for query alignment and re-alignment.
 *
 * @author nils.homer@lifetech.com
 */
public class AlignRecord implements Cloneable {

    /**
     * The SAM record associated with this query.
     */
    public SAMRecord record;
    /**
     * The zero-based input file index from which this query was read.
     */
    public int fileIndex; // to which input/output file does this belong?
    /**
     * The one-based start position of the alignment, where we include any upstream
     * matching bases as to always start at the beginning of a homopolymer.
     */
    public int positionStart;
    /**
     * The one-based end position of the alignment.
     */
    public int positionEnd; // one-based
    /**
     * The zero-based reference sequence index.
     */
    public int referenceIndex; // zero-based
    /**
     * The strand (true - reverse, false - forward).
     */
    public boolean strand;
    /**
     * The flow order for the read, phased to match the first sequenced base
     * in the read.
     */
    public FlowOrder flowOrder;
    /**
     * The flow signals for the read according to the flow order, in sequencing
     * order.
     */
    public FlowSeq readFlows;
    /**
     * The read sequence, with soft-clipping inforation partitioned, relative
     * to the genomic forward strand.
     */
    public ReadSeq readSeq;
    /**
     * The read bases in integer format relative to the forward genomic strand.
     */
    public byte[] readBases;
    /**
     * The flow space representation of the alignment relative to the forward
     * genomic strand.
     */
    public FlowSpaceAlignment alignment;
    /**
     * True if we are to perform a base space start local alignment.
     */
    public int flowOrderIndexStart;

    /**
     * @param record the original SAM record.
     * @param fileIndex the zero-based file index from which the query was read.
     */
    public AlignRecord(SAMRecord record, int fileIndex, FlowOrder flowOrder) {
        this.record = record;
        this.fileIndex = fileIndex;
        this.readSeq = new ReadSeq(record);
        this.readBases = readSeq.getReadBytes();
        this.referenceIndex = record.getReferenceIndex();
        this.positionStart = record.getAlignmentStart();
        this.positionEnd = record.getAlignmentEnd();
        this.strand = record.getReadNegativeStrandFlag();
        this.flowOrder = flowOrder;
        this.alignment = null;
    }

    
    public AlignRecord(AlignRecord toCopy) {
     
    }
    
    /**
     * Note that this does a shallow copy of everything, except for the attribute list, for which a copy of the list
     * is made, but the attributes themselves are copied by reference.  This should be safe because callers should
     * never modify a mutable value returned by any of the get() methods anyway.
     */
    @Override
    public Object clone() throws CloneNotSupportedException {
        final AlignRecord newRecord = (AlignRecord)super.clone();
        newRecord.fileIndex = this.fileIndex;
        newRecord.referenceIndex = this.referenceIndex;
        newRecord.flowOrder = this.flowOrder;
        newRecord.alignment = this.alignment;
        try{
            newRecord.record = (SAMRecord)this.record.clone();
            newRecord.readSeq = new ReadSeq(this.record);
            newRecord.readBases = readSeq.getReadBytes();
            newRecord.positionStart = this.record.getAlignmentStart();
            newRecord.positionEnd = this.record.getAlignmentEnd();
            newRecord.strand = this.record.getReadNegativeStrandFlag();
            
        } catch (Exception e) {
            
        }
        return newRecord;
    }
    
    
    /**
     * @param start true if we are to examine from the start of the alignment, false if we want to examine from the end of the alignment.
     * @return the first non-empty target flow in the alignment.
     */
    private int getTSeqJ(boolean start) {
        int j;
        if(start) {
            j = 0;
            while(FlowSpaceAlignment.ALN_INS == this.alignment.aln[j] || SamToFlowSpaceUtil.getBaseCallFromFlowSignal(this.alignment.tseq[j]) <= 0) {
                j++;
                if(j == this.alignment.length) {
                    break;
                }
            }
        }
        else {
            j = this.alignment.length - 1;
            while(FlowSpaceAlignment.ALN_INS == this.alignment.aln[j] || SamToFlowSpaceUtil.getBaseCallFromFlowSignal(this.alignment.tseq[j]) <= 0) {
                j--;
                if(j == 0) {
                    break;
                }
            }
        }
        return j;
    }

    /**
     * Adjust the start position based on hps.
     * @param sequence the reference sequence.
     * @param adjustTseq adjust the tseq value in the alignment.
     */
    private void adjustPositionStart(ReferenceSequence sequence, boolean adjustTseq)
        throws Exception {
        int i, j;
        byte base, curBase;
        if(1 < this.positionStart) {
            if(!sequence.isIupac(this.referenceIndex, this.positionStart)) {
                base = sequence.getBase(this.referenceIndex, this.positionStart, true, false);
                i = 0;
                while(1 <= this.positionStart - i - 1) { // examine the previous bases
                    curBase = sequence.getBase(this.referenceIndex, this.positionStart - 1 - i, true, false);
                    if(curBase != base) {
                        break;
                    }
                    i++;
                }
                // adjust tseq
                if(0 < i) {
                    this.positionStart -= i;
                    if(adjustTseq) {
                        j = getTSeqJ(true);
                        this.alignment.tseq[j] += i * 100;
                    }
                }
            }
        }
    }

    /**
     * Adjust the end position based on hps.
     * @param sequence the reference sequence.
     * @param adjustTseq adjust the tseq value in the alignment.
     */
    private void adjustPositionEnd(ReferenceSequence sequence, boolean adjustTseq)
        throws Exception {
        int i, j;
        byte base, curBase;
        if(this.positionEnd < sequence.length(this.referenceIndex)) {
            if(!sequence.isIupac(this.referenceIndex, this.positionEnd)) {
                base = sequence.getBase(this.referenceIndex, this.positionEnd, true, false);
                i = 0;
                while(this.positionEnd + i + 1 <= sequence.length(this.referenceIndex)) {
                    curBase = sequence.getBase(this.referenceIndex, this.positionEnd + i + 1, true, false);
                    if(curBase != base) {
                        break;
                    }
                    i++;
                }
                // adjust tseq
                if(0 < i) {
                    this.positionEnd += i;
                    if(adjustTseq) {
                        j = getTSeqJ(false);
                        this.alignment.tseq[j] += i * 100;
                    }
                }
            }
        }
    }
    
    public FlowSpaceAlignment setAlignment(ReferenceSequence sequence, int phasePenalty)
        throws Exception
    {
        return this.setAlignment(sequence, phasePenalty, false);
    }
    
    /**
     * Fills in the alignment and flow information for this read.
     * @param sequence the reference sequence.
     * @param phasePenalty the flow space phase penalty.
     * @param validateAlignments validate alignments.
     * @return the flow alignment.
     */
    public FlowSpaceAlignment setAlignment(ReferenceSequence sequence, int phasePenalty, boolean validateAlignments)
        throws Exception
    {
        return this.setAlignment(sequence, phasePenalty, 0, validateAlignments);
    }

    /**
     * Fills in the alignment and flow information for this read.
     * @param sequence the reference sequence.
     * @param phasePenalty the flow space phase penalty.
     * @param offset the alignment offset in bases
     * @return the flow alignment.
     */
    public FlowSpaceAlignment setAlignment(ReferenceSequence sequence, int phasePenalty, int offset, boolean validateAlignments)
        throws Exception 
    {
        int i, j, k, n;
        byte base, curBase;
        byte[] referenceBases = null;
        byte[] readBases = null;
        int[] flowSignals = null;
        int[] tmpFlowSignals = null;
        int[] hardClipBases = null;
        boolean tmpBoolean;
        boolean startLocalAlignment, endLocalAlignment;
        int score;
        
        this.alignment = null;

        // default to local alignment
        startLocalAlignment = endLocalAlignment = true; // prove otherwise

        // NB: we could extend to the previous/next homopolymers
        this.positionStart -= offset;
        if(this.positionStart < 1) {
            this.positionStart = 1;
        }
        this.positionEnd += offset;
        if(sequence.length(this.referenceIndex) < this.positionEnd) {
            this.positionEnd = sequence.length(this.referenceIndex);
        }

        // extend the start
        base = sequence.getBase(this.referenceIndex, this.positionStart, true, false);
        while(1 < this.positionStart && base < 4) {
            curBase = sequence.getBase(this.referenceIndex, this.positionStart-1, true, false);
            if(4 <= curBase || curBase != base) {
                break;
            }
            this.positionStart--;
        }
        // extend the end
        base = sequence.getBase(this.referenceIndex, this.positionEnd, true, false);
        while(this.positionEnd < sequence.length(this.referenceIndex) && base < 4) {
            curBase = sequence.getBase(this.referenceIndex, this.positionEnd+1, true, false);
            if(4 <= curBase || curBase != base) {
                break;
            }
            this.positionEnd++;
        }

        // read and reference bases for the flow space alignment
        readBases = new byte[this.readBases.length];
        System.arraycopy(this.readBases, 0, readBases, 0, this.readBases.length);
        referenceBases = sequence.getBases(this.referenceIndex, this.positionStart, this.positionEnd);
        
        // Note: readBases and referenceBases are on the forward strand
        // Reverse compliment
        if(this.strand) {
            SamToFlowSpaceUtil.reverseCompliment(readBases);
            SamToFlowSpaceUtil.reverseCompliment(referenceBases);
            // swap
            tmpBoolean = startLocalAlignment;
            startLocalAlignment = endLocalAlignment;
            endLocalAlignment = tmpBoolean;
        }

        /*
        System.err.println("this.positionStart=" + this.positionStart);
        System.err.println("this.positionEnd=" + this.positionEnd);
        for(i=0;i<readBases.length;i++) {
        System.err.print(SamToFlowSpaceUtil.DNA[(int)readBases[i]]);
        }
        System.err.println("");
        for(i=0;i<referenceBases.length;i++) {
        System.err.print(SamToFlowSpaceUtil.DNA[(int)referenceBases[i]]);
        }
        System.err.println("");
        */

        // retrieve the flow signals
        flowSignals = SamToFlowSpaceUtil.getFlowSignals(this.record);
        if(null == flowSignals) {
            throw new Exception("The FZ optional tag (flow signals) was not present");
        }

        // a. if there are hard clip bases, retrieve them from the flow signals
        // b. get the flow order start index, given the key sequence, hard clip bases, and soft clip bases
        // c. adjust the flow signals based on the start index, and last non-read base (last base of the key, or last base of the hard clip)
        // TODO
        flowSignals = this.readSeq.updateFlowInformation(this.strand, this.flowOrder, flowSignals);
        this.flowOrderIndexStart = this.readSeq.getFlowOrderIndexStart();

        // rotate the flow order to match the first flow in readBases for the flow space pre-alignment
        this.flowOrder.rotate(this.flowOrderIndexStart, true);

        // read flows
        //this.readFlows = new FlowSeq(readBases, this.flowOrder.flowOrder);
        this.readFlows = new FlowSeq(readBases, flowSignals, this.flowOrder.flowOrder);

        // represent the alignment in flow space
        this.alignment = new FlowSpaceAlignment(this.readFlows, referenceBases, this.flowOrder,
            startLocalAlignment,
            endLocalAlignment,
            phasePenalty);

        //this.alignment.print(System.err);

        // swap back
        if(this.strand) {
            // swap
            tmpBoolean = startLocalAlignment;
            startLocalAlignment = endLocalAlignment;
            endLocalAlignment = tmpBoolean;
            // reverse the alignment
            this.alignment.reverseCompliment();
        }
          
        // adjust start position and end position based on alignment
        if(0 < this.alignment.tseqStart) {
            this.positionStart += this.alignment.tseqStart;
        }
        if(this.alignment.tseqEnd < this.alignment.tseqLength - 1) {
            this.positionEnd -= (this.alignment.tseqLength - this.alignment.tseqEnd - 1);
        }

        // adjust if we start with an insertion
        if(FlowSpaceAlignment.ALN_INS == this.alignment.aln[0]) {
            this.positionStart--;
        }
        
        if(validateAlignments) {
            // get the score before validating
            score = recalculateAlignmentScore(phasePenalty);
        
            // This is required for the graph ordering to work.
            this.alignment.rightAdjustIndels();
        
            // validate the alignment
            validateAlignment(sequence, score);
        }
        else {
            // This is required for the graph ordering to work.
            this.alignment.rightAdjustIndels();
        }

        return this.alignment;
    }

    private String validateAlignmentTestOne()
    {
        /** 
         * Test 1
         * Check that the target length is equal to the flow space alignment target length
         */
        int i, length;
        int positionStart;

        for(i=length=0;i<this.alignment.length;i++) {
            if(FlowSpaceAlignment.ALN_INS != this.alignment.aln[i]) {
                length += this.alignment.tseq[i];
            }
        }
        // convert from flow signals
        length = SamToFlowSpaceUtil.getBaseCallFromFlowSignal(length);

        positionStart = this.positionStart;
        if(FlowSpaceAlignment.ALN_INS == this.alignment.aln[0]) {
            positionStart++;
        }

        if(this.positionEnd - positionStart + 1 != length) {
            return new String("Test one: target lengths do not match!\n(" + this.positionEnd + " - " + positionStart + " + 1) != " + length + "\n");
        }
        else {
            return "";
        }
    } 
    
    private String validateAlignmentTestTwo(ReferenceSequence sequence)
        throws Exception
    {
        /** 
         * Test 2
         * Check that the target sequence is equal to the target sequence implied by the flow space alignment
         */
        int i, j, length;
        byte[] referenceBases = null;
        byte[] fsReferenceBases = null;
        boolean differ = false;
        int positionStart;

        positionStart = this.positionStart;
        if(FlowSpaceAlignment.ALN_INS == this.alignment.aln[0]) {
            positionStart++;
        }
        
        referenceBases = sequence.getBases(this.referenceIndex, positionStart, this.positionEnd, false, true);

        fsReferenceBases = new byte[this.positionEnd - this.positionStart + 1];
        for(i=j=0;i<this.alignment.length;i++) {
            if(FlowSpaceAlignment.ALN_INS != this.alignment.aln[i]) {
                length = SamToFlowSpaceUtil.getBaseCallFromFlowSignal(this.alignment.tseq[i]);
                while(0 < length) {
                    fsReferenceBases[j] = (byte)SamToFlowSpaceUtil.DNA[this.alignment.flowOrder[i]];
                    if(referenceBases.length <= j || fsReferenceBases[j] != referenceBases[j]) {
                        differ = true;
                    }
                    j++;
                    length--;
                }
            }
        }

        if(differ) {
            return new String("Test two: target sequences do not match!\n" + new String(referenceBases) + "\n" + new String(fsReferenceBases) + "\n");
        }
        else {
            return "";
        }
    } 
    
    private String validateAlignmentTestThree(ReferenceSequence sequence)
        throws Exception
    {
        /** 
         * Test 3
         * Check that the first and last target homopolymers are fully extended.
         */
        byte curBase;
        
        positionStart = this.positionStart;
        if(FlowSpaceAlignment.ALN_INS == this.alignment.aln[0]) {
            positionStart++;
        }
        
        curBase = sequence.getBase(this.referenceIndex, positionStart);
        if(1 < positionStart && !sequence.isIupac(this.referenceIndex, positionStart) && sequence.getBase(this.referenceIndex, positionStart) == sequence.getBase(this.referenceIndex, positionStart-1))  {
            return new String("Test three: target sequence is not fully extended at the start (" + positionStart + ")!\n");
        }
        if(this.positionEnd < sequence.length(this.referenceIndex) && !sequence.isIupac(this.referenceIndex, this.positionEnd) && sequence.getBase(this.referenceIndex, this.positionEnd) == sequence.getBase(this.referenceIndex, this.positionEnd+1))  {
            return new String("Test three: target sequence is not fully extended at the end (" + this.positionEnd + ")!\n");
        }
        return "";
    } 
    
    private String validateAlignmentTestFour()
        throws Exception
    {
        /** 
         * Test 4
         * Check that the query's flow signals match the record flow signals.
         */
        int i, j, length;
        int[] flowSignals = null;
        int[] fsFlowSignals = null;
        boolean differ = false;

        length = this.readFlows.nonEmptyFlowLast - this.readFlows.nonEmptyFlowFirst + 1;
        flowSignals = new int[length];
        System.arraycopy(this.readFlows.flow, this.readFlows.nonEmptyFlowFirst, flowSignals, 0, length);
        if(this.strand) {
            // swap
            for(i=0;i<flowSignals.length / 2;i++) {
                j = flowSignals[i];
                flowSignals[i] = flowSignals[flowSignals.length-i-1];
                flowSignals[flowSignals.length-i-1] = j;
            }
        }

        fsFlowSignals = new int[this.alignment.length];
        for(i=j=0;i<this.alignment.length;i++) {
            if(FlowSpaceAlignment.ALN_DEL != this.alignment.aln[i]) {
                if(flowSignals[j] != this.alignment.qseq[i]) {
                    differ = true;
                }
                fsFlowSignals[j] = this.alignment.qseq[i];
                j++;
            }
        }
        length = j;

        if(differ) {
            String s = new String("Test four: flow signals differ!\n");
            for(i=0;i<flowSignals.length;i++) {
                if(0 < i) {
                    s += ",";
                }
                s += flowSignals[i];
            }
            s += "\n";
            for(i=0;i<length;i++) {
                if(0 < i) {
                    s += ",";
                }
                s += fsFlowSignals[i];
            }
            s += "\n";
            return s;
        }
        else {
            return "";
        }
    } 

    // NB: compute this before right-justifying
    private int recalculateAlignmentScore(int phasePenalty)
        throws Exception
    {
        int i, score;
        byte prevBase = 'N';
        int n;

        int start, end, by;

        score = n = 0;
        if(strand) { // reverse
            start=this.alignment.length-1;
            end=-1;
            by=-1;
        }
        else {
            start=0;
            end=this.alignment.length;
            by=1;
        }
        for(i=start;i!=end;i+=by) {
            if(FlowSpaceAlignment.ALN_INS == this.alignment.aln[i]) {
                score -= this.alignment.qseq[i];
                if(prevBase == this.alignment.flowOrder[i]) {
                    n++;
                    score -= phasePenalty; 
                }
            }
            else if(FlowSpaceAlignment.ALN_DEL == this.alignment.aln[i]) {
                score -= this.alignment.tseq[i];
                prevBase = this.alignment.flowOrder[i];
            }
            else if(FlowSpaceAlignment.ALN_MATCH == this.alignment.aln[i] || FlowSpaceAlignment.ALN_MISMATCH == this.alignment.aln[i]) {
                if(0 == i || this.alignment.length - 1 == i) {
                    int s = SamToFlowSpaceUtil.getFlowSignalFromBaseCall(SamToFlowSpaceUtil.getBaseCallFromFlowSignal(this.alignment.qseq[i]));
                    if(s < this.alignment.qseq[i]) {
                        score -= (this.alignment.qseq[i] - s);
                    }
                    else {
                        score -= (s - this.alignment.qseq[i]);
                    }
                }
                else {
                    if(this.alignment.tseq[i] <  this.alignment.qseq[i]) {
                        score -= (this.alignment.qseq[i] - this.alignment.tseq[i]);
                    }
                    else {
                        score -= (this.alignment.tseq[i] - this.alignment.qseq[i]);
                    }
                }
                prevBase = this.alignment.flowOrder[i];
            }
            else {
                throw new Exception("Unknown operator: " + this.alignment.aln[i]);
            }
            /*
               System.err.println("score=" + score
               + " " + this.alignment.qseq[i]
               + " " + this.alignment.aln[i]
               + " " + this.alignment.tseq[i]
               + " " + this.alignment.flowOrder[i]);
               */
        }
        return score;
    }
    
    private String validateAlignmentTestFive(int score)
        throws Exception
    {
        /** 
         * Test 5
         * Check that the alignment score is correct
         */
        if(score != this.alignment.getScore()) {
            return new String("Test five: alignment score was wrong (" + score + " != " + this.alignment.getScore() + ")!\n");
        }
        else {
            return "";
        }
    } 

    public void validateAlignment(ReferenceSequence sequence, int score)
        throws Exception
    {
        String out = new String("");

        // run tests
        out += validateAlignmentTestOne();
        out += validateAlignmentTestTwo(sequence);
        out += validateAlignmentTestThree(sequence);
        out += validateAlignmentTestFour();
        out += validateAlignmentTestFive(score);

        if(0 < out.length()) {
            out += this.getSAMString(); 
            out += this.alignment.getAlignmentString(150, true);
            throw new Exception(out);
        }
    }

    public void clear() {
        this.record = null;
        this.flowOrder = null;
        this.readFlows = null;
        this.readSeq = null;
        this.readBases = null;
        this.alignment = null;
    }
    
     public String getSAMString() {
        StringBuilder sb = new StringBuilder();
        //sb.append( 
        String FIELD_SEPARATOR = "\t";
        sb.append(this.record.getReadName());
        sb.append(FIELD_SEPARATOR);
        sb.append(Integer.toString(this.record.getFlags()));
        sb.append(FIELD_SEPARATOR);
        sb.append(this.record.getReferenceName());
        sb.append(FIELD_SEPARATOR);
        sb.append(Integer.toString(this.record.getAlignmentStart()));
        sb.append(FIELD_SEPARATOR);
        sb.append(Integer.toString(this.record.getMappingQuality()));
        sb.append(FIELD_SEPARATOR);
        sb.append(this.record.getCigarString());
        sb.append(FIELD_SEPARATOR);

        //  == is OK here because these strings are interned
        if (this.record.getReferenceName() == this.record.getMateReferenceName() &&
                SAMRecord.NO_ALIGNMENT_REFERENCE_NAME != this.record.getReferenceName()) {
            sb.append("=");
        } else {
            sb.append(this.record.getMateReferenceName());
        }

        sb.append(FIELD_SEPARATOR);
        sb.append(Integer.toString(this.record.getMateAlignmentStart()));
        sb.append(FIELD_SEPARATOR);
        sb.append(Integer.toString(this.record.getInferredInsertSize()));
        sb.append(FIELD_SEPARATOR);
        sb.append(this.record.getReadString());
        sb.append(FIELD_SEPARATOR);
        sb.append(this.record.getBaseQualityString());
        //SAMBinaryTagAndValue attribute = this.record.getBinaryAttributes();
        
        //StringBuilder mdTag = new StringBuilder();
        String[] stringAttributes       = {"MD", "PG", "RG" }; 
        String[] integerAttributes      = {"NM", "AS", "XS", "XT" };
        String[] byteArrayAttributes    = {"FZ"};
        String stringAttr = ":Z:";
        String intAttr = ":i:";
        String byteArrayAttr = ":B:S,";
        
        for (String attr : stringAttributes) {
            sb.append(FIELD_SEPARATOR);
            sb.append(attr);
            sb.append(stringAttr);
            //System.out.println( attr+":Z:" + this.record.getStringAttribute( attr ) );
            sb.append( this.record.getStringAttribute(attr) );
        }

        for (String attr : integerAttributes) {
            sb.append( FIELD_SEPARATOR );
            sb.append(attr);
            sb.append(intAttr);
            //System.out.println( attr+":i:" + this.record.getIntegerAttribute( attr ) );

            sb.append( this.record.getIntegerAttribute(attr) );
        }
       
        int[] flowSignals = SamToFlowSpaceUtil.getFlowSignals(this.record);
        sb.append( FIELD_SEPARATOR );
        sb.append("FZ");
        sb.append(byteArrayAttr);
        for (int flow : flowSignals) {
            sb.append(flow);
            sb.append(",");
        }
        
        sb.append("\n");
        
        return sb.toString();
     }
}
