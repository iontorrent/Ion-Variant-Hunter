/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2fs.util;

import java.io.*;
import java.util.*;
import net.sf.samtools.*;
import net.sf.samtools.*;
import net.sf.picard.reference.*;

/**
 * Wrapper for the Picard reference bases, to enable IUPAC supprt.
 *
 * @author nils.homer@lifetech.com
 */
public class ReferenceSequence {

    /**
     * The reference sequence file.
     */
    private net.sf.picard.reference.ReferenceSequenceFile referenceSequenceFile = null;

    /**
     * The current reference sequence itself, one contig only.
     */
    private net.sf.picard.reference.ReferenceSequence referenceSequence = null;

    /**
     * The reference sequences, all contigs.
     */
    private List<net.sf.picard.reference.ReferenceSequence> referenceSequences = null;

    /**
     * The reference dictionary.
     */
    private net.sf.samtools.SAMSequenceDictionary referenceDictionary = null;

    /**
     * True if we are to only store one contig at a time, false otherwise
     */
    private boolean oneAtATime = false;
    
    /**
     * The reference sequence file must have an accompanying index (.fai) 
     * and a dictionary file (.dict).
     * @param input the input reference sequence file.
     * @param oneAtATime true if we are to only store one contig at a time, false otherwise
     */
    public ReferenceSequence(File input, boolean oneAtATime)
        throws Exception
    {
        this.referenceSequenceFile = new IndexedFastaSequenceFile(input);
        if(!this.referenceSequenceFile.isIndexed()) {
            throw new Exception("Reference sequence file was not indexed.");
        }
        this.referenceDictionary = this.referenceSequenceFile.getSequenceDictionary();
        if(null == this.referenceDictionary) {
            throw new Exception("Could not find FASTA dictionary file.");
        }
        this.oneAtATime = oneAtATime;
        if(!this.oneAtATime) {
            this.referenceSequences = new ArrayList<net.sf.picard.reference.ReferenceSequence>();
        }
    }

    /**
     * The reference sequence file must have an accompanying index (.fai) 
     * and a dictionary file (.dict).
     * @param input the input reference sequence file.
     */
    public ReferenceSequence(File input)
        throws Exception
    {
        this(input, true);
    }

    /**
     * Moves the reference sequence to a new contig.
     * @param referenceIndex the zero-based reference index to which to move.
     */
    public void moveTo(int referenceIndex)
        throws Exception
    {
        int i;
        if(this.oneAtATime) {
            if(null == this.referenceSequence || referenceIndex != this.referenceSequence.getContigIndex()) {
                this.referenceSequence = this.referenceSequenceFile.getSequence(this.referenceDictionary.getSequence(referenceIndex).getSequenceName());
                if(null == this.referenceSequence) {
                    throw new Exception("Premature EOF in the reference sequence");
                }
                else if(this.referenceSequence.getContigIndex() != referenceIndex) {
                    throw new Exception("Could not find the reference sequence");
                }
            }
        }
        else {
            // add them up to the current
            for(i=this.referenceSequences.size(); i <= referenceIndex; i++) {
                this.referenceSequences.add(this.referenceSequenceFile.getSequence(this.referenceDictionary.getSequence(i).getSequenceName()));
            }
        }
    }

    /**
     * @return the reference sequence dictionary.
     */
    public net.sf.samtools.SAMSequenceDictionary getDictionary()
    {
        return this.referenceDictionary;
    }

    /**
     * @param referenceIndex the zero-based reference index to which to move.
     * @return true if the contig is contained, false otherwise
     */
    public boolean contains(int referenceIndex)
    {
        if(this.oneAtATime) {
            if(null == this.referenceSequence || this.referenceSequence.getContigIndex() != referenceIndex) {
                return false;
            }
        }
        else {
            if(this.referenceSequences.size() <= referenceIndex) {
                return false;
            }
        }
        return true;
    }
    
    public boolean isIupac(int referenceIndex, int position)
        throws Exception
    {
        return (4 <= getBase(referenceIndex, position, true, false));
    }

    /**
     * The base will be returned in integer form.  If the base is an non-DNA IUPAC base,
     * it will be converted to the lexicographically smallest non-matching base.  In the 
     * case of an N, it will be converted to an A.
     * @param referenceIndex the zero-based reference index to which to move.
     * @param position the one-based position.
     * @return the base at the given position.
     */
    public byte getBase(int referenceIndex, int position)
        throws Exception
    {
        return this.getBase(referenceIndex, position, true, true);
    }

    /**
     * @param referenceIndex the zero-based reference index to which to move.
     * @param position the one-based position.
     * @param asInt true if the base is to be returned in integer form, false otherwise.
     * @param iupac true if the base is to be converted to the lexicographically
     * smallest non-matching base (for iupac bases), or in the case of an N, it will be
     * converted to an A, otherwise false.
     * @return the base at the given position.
     */
    public byte getBase(int referenceIndex, int position, boolean asInt, boolean iupac)
        throws Exception
    {
        byte base;
    
        if(!this.contains(referenceIndex)) {
            throw new Exception("Contig not found in the reference:" + referenceIndex);
        }
        
        if(this.oneAtATime) {
            base = this.referenceSequence.getBases()[position-1];
        }
        else {
            base = this.referenceSequences.get(referenceIndex).getBases()[position-1];
        }

        if(iupac) {
            base = SamToFlowSpaceUtil.IUPAC2NT[base];
        }
        if(asInt) {
            base = SamToFlowSpaceUtil.NT2INT[(int)base];
        }

        return base;
    }

    /**
     * The bases will be returned in integer form.  If any base is an non-DNA IUPAC base,
     * it will be converted to the lexicographically smallest non-matching base.  In the 
     * case of an N, it will be converted to an A.
     * @param referenceIndex the zero-based reference index to which to move.
     * @param start int the one-based start position.
     * @param end int the one-based end position.
     * @return a new array of bases.
     */
    public byte[] getBases(int referenceIndex, int start, int end)
        throws Exception
    {
        return this.getBases(referenceIndex, start, end, true, true);
    }

    /**
     * @param referenceIndex the zero-based reference index to which to move.
     * @param start int the one-based start position.
     * @param end int the one-based end position.
     * @param asInt true if the base is to be returned in integer form, false otherwise.
     * @param iupac true if the base is to be converted to the lexicographically
     * smallest non-matching base (for iupac bases), or in the case of an N, it will be
     * converted to an A, otherwise false.
     * @return a new array of bases.
     */
    public byte[] getBases(int referenceIndex, int start, int end, boolean asInt, boolean iupac)
        throws Exception
    {
        byte[] bases = null;
        bases = new byte[end - start + 1];

        if(!this.contains(referenceIndex)) {
            throw new Exception("Contig not found in the reference:" + referenceIndex);
        }

        if(this.oneAtATime) {
            System.arraycopy(this.referenceSequence.getBases(), start-1, bases, 0, bases.length);
        }
        else {
            System.arraycopy(this.referenceSequences.get(referenceIndex).getBases(), start-1, bases, 0, bases.length);
        }

        if(iupac) {
            SamToFlowSpaceUtil.iupacToNT(bases);
        }
        if(asInt) {
            SamToFlowSpaceUtil.ntToInt(bases);
        }
        return bases;
    }

    /**
     * @param referenceIndex the zero-based reference index to which to move.
     * @return the length of the current reference sequence.
     */
    public int length(int referenceIndex)
        throws Exception
    {
        if(!this.contains(referenceIndex)) {
            throw new Exception("Contig not found in the reference:" + referenceIndex);
        }
        if(this.oneAtATime) {
            return this.referenceSequence.length();
        }
        else {
            return this.referenceSequences.get(referenceIndex).length();
        }
    }

}  
