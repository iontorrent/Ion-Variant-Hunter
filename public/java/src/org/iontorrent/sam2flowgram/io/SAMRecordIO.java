/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2flowgram.io;

import org.iontorrent.sam2flowgram.util.*;
import org.iontorrent.sam2flowgram.flowalign.FlowOrder;

import java.util.*;
import java.io.*;
import net.sf.samtools.*;
import net.sf.samtools.util.*;
import net.sf.picard.sam.*;
import net.sf.picard.io.IoUtil;

/**
 * Handles reading/writing to/from SAM/BAM files.
 *
 * @author nils.homer@lifetech.com
 */
public class SAMRecordIO 
{
    private static final String HEADER_FILE_INDEX_TAG = "XI";

    /**
     * The program record for the output.
     */
    public SAMProgramRecord programRecord;

    /**
     * The list of input SAM/BAM file readers.
     */
    private List<SAMFileReader> readers;
    
    /**
     * The list of SAM headers from the input SAM/BAM file readers.
     */
    private List<SAMFileHeader> readersHeaders;

    /**
     * The flow orders for each input SAM/BAM.
     */
    private List<Map<String,FlowOrder>> flowOrders;

    /**
     * The list of output SAM/BAM file writers.
     */
    private List<SAMFileWriter> writers;

    /**
     * The SAM/BAM output writers for unmapped reads.
     */
    private List<SAMFileWriter> writersUnmapped;

    /**
     * The list of record iterators for each input SAM/BAM file.
     */
    private List<CloseableIterator<SAMRecord>> recordsIters = null;

    /**
     * A record buffer for each input SAM/BAM file.
     */
    private List<FlowAlignRecord> buffer = null;

    /**
     * True if header comments were added to the output file, false otherwise.
     */
    private boolean headerCommentsAdded = false;

    /**
     * Opens the input and output file(s), merges and modifies the header(s), and initializes the buffer(s) for processing.
     * @param inputs the input SAM/BAM file(s).
     * @param outputs the output SAM/BAM file(s).
     * @param outputsUnmapped the output SAM/BAM file(s) for unmapped reads.
     * @param programVeresion the SamToFlowgramAlign program version.
     * @param useRanges true if the SAM/BAM file is to be range queried, false otherwise.
     * @param referenceDictionary the reference sequence dictionary.
     */
    public SAMRecordIO(List<File> inputs, 
            List<File> outputs, 
            List<File> outputsUnmapped,
            String programVersion, 
            boolean useRanges, 
            SAMSequenceDictionary referenceDictionary)
        throws Exception
    {
        ListIterator<File> inputsIter = null;
        ListIterator<File> outputsIter = null;
        ListIterator<File> outputsUnmappedIter = null;
        ListIterator<SAMFileReader> readersIter = null;
        ListIterator<String> headerCommentsIter = null;
        SAMFileHeader mergedHeader = null;
        List<String> headerComments = null;
        int fileIndex = 0;

        this.readers = new ArrayList<SAMFileReader>();
        this.readersHeaders = new ArrayList<SAMFileHeader>();
        this.flowOrders = new ArrayList<Map<String,FlowOrder>>();
        this.writers = new ArrayList<SAMFileWriter>();
        this.writersUnmapped = new ArrayList<SAMFileWriter>();
        
        headerComments = new ArrayList<String>();

        programVersion = new String("trap-" + programVersion); // append "trap-" so we know it was trap

        inputsIter = inputs.listIterator();
        if(null != outputs && 1 < outputs.size()) { // to multiple files 
            outputsIter = outputs.listIterator();
            if(outputs.size() != inputs.size()) {
                throw new Exception("There must be the same # of inputs as outputs");
            }
        }
        if(null != outputsUnmapped && 1 < outputsUnmapped.size()) { // to multiple files 
            outputsUnmappedIter = outputsUnmapped.listIterator();
            if(outputsUnmapped.size() != inputs.size()) {
                throw new Exception("There must be the same # of inputs as unmapped outputs");
            }
        }
        while(inputsIter.hasNext()) {
            File file = inputsIter.next();
            SAMFileReader fileReader = null;
            SAMFileHeader fileHeader = null;
            String headerComment = null;
            List<SAMReadGroupRecord> readGroups = null;
            String flowOrder = null, keySequence = null;
            Map<String,FlowOrder> readGroupToFlowOrder = null;

            IoUtil.assertFileIsReadable(file);

            fileReader = new SAMFileReader(file, true);
            if(useRanges && !fileReader.hasIndex()) {
                throw new Exception("BAM files and BAM indexes are required when using the RANGE or RANGES option"); 
            }
            fileHeader = fileReader.getFileHeader();
            if(!checkHeaderAgainstReferenceDictionary(fileHeader, referenceDictionary)) {
                throw new Exception("FASTA sequence dictionary and SAM/BAM file dictionary are in different orders");
            }

            this.programRecord = fileHeader.getProgramRecord("trap");
            
            // Flow order and key sequence
            readGroups = fileHeader.getReadGroups();
            if(0 == readGroups.size()) {
                throw new Exception("Could not find a read group in the SAM/BAM header");
            }
            readGroupToFlowOrder = new HashMap<String,FlowOrder>();
            for(SAMReadGroupRecord readGroup : readGroups) {
                flowOrder = readGroup.getFlowOrder();
                if(null == flowOrder) {
                    throw new Exception("Could not find the flow order in the SAM/BAM header");
                }
                keySequence = readGroup.getKeySequence();
                if(null == keySequence) {
                    throw new Exception("Could not find the key sequence in the SAM/BAM header");
                }
                readGroupToFlowOrder.put(readGroup.getId(), new FlowOrder(flowOrder, keySequence));
            }
            this.flowOrders.add(readGroupToFlowOrder);

            if(null == programRecord) { // create a new one
                this.programRecord = new SAMProgramRecord("trap");
                this.programRecord.setProgramVersion(programVersion);
                fileHeader.addProgramRecord(this.programRecord);
            }
            else if(0 != programVersion.compareTo(programRecord.getProgramVersion())) { // new version, but trap exists
                this.programRecord = fileHeader.createProgramRecord();
                this.programRecord.setProgramVersion(programVersion);
            }

            // Always set to coordinate sorted
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

            // Create the header comment
            headerComment = new String("@CO\t" + HEADER_FILE_INDEX_TAG + ":" + fileIndex + "\tPATH:" + file.getPath());

            this.readers.add(fileReader);
            this.readersHeaders.add(fileHeader);
            headerComments.add(headerComment);
            if(null != outputs && 1 < outputs.size()) { // to multiple files 
                this.writers.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(fileHeader, true, outputsIter.next())); 
            }
            if(null != outputsUnmapped && 1 < outputsUnmapped.size()) { // to multiple files
                this.writersUnmapped.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(fileHeader, true, outputsUnmappedIter.next()));
            }
            fileIndex++;
        }

        // Merge headers
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, this.readersHeaders, true);
        mergedHeader = headerMerger.getMergedHeader();
        // Always set to coordinate sorted
        mergedHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        // Add the header commments
        if(1 < inputs.size()) {
            headerCommentsIter = headerComments.listIterator();
            while(headerCommentsIter.hasNext()) {
                mergedHeader.addComment(headerCommentsIter.next());
            }
            headerCommentsAdded = true;
        }

        if(null != outputs) {
            if(0 == outputs.size()) { // to STDOUT
                this.writers.add(new SAMFileWriterFactory().makeSAMWriter(mergedHeader, true, System.out));
            }
            else if(1 == outputs.size()) { // one output file
                this.writers.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(mergedHeader, true, outputs.get(0)));
            }
        }

        if(null != outputsUnmapped) {
            if(0 == outputsUnmapped.size()) {
                // suppress
            }
            else if(1 == outputsUnmapped.size()) { // one output file
                this.writersUnmapped.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(mergedHeader, true, outputsUnmapped.get(0)));
            }
        }

        // Default iterators
        this.recordsIters = new ArrayList<CloseableIterator<SAMRecord>>();
        readersIter = this.readers.listIterator();
        while(readersIter.hasNext()) {
            this.recordsIters.add(readersIter.next().iterator());
        }
    }

    /**
     * Validates the sequence dictionary in the header versus the given reference sequence dictionary.
     * @param header the SAM/BAM header
     * @param referenceDictionary the reference sequence dictionary.
     * @return true if they are the same, false otherwise.
     */
    private boolean checkHeaderAgainstReferenceDictionary(SAMFileHeader header,  SAMSequenceDictionary referenceDictionary)
    {
        int i;
        SAMSequenceDictionary headerDict;
        if(null == header || null == referenceDictionary) {
            return true;
        }

        headerDict = header.getSequenceDictionary();
        if(headerDict.size() != referenceDictionary.size()) {
            return false;
        }
        for(i=0;i<headerDict.size();i++) {
            // check name and length
            SAMSequenceRecord headerRec = headerDict.getSequence(i);
            SAMSequenceRecord referenceRec = referenceDictionary.getSequence(i);

            if(!headerRec.getSequenceName().equals(referenceRec.getSequenceName()) 
                    || headerRec.getSequenceLength() != referenceRec.getSequenceLength()) 
            {
                return false;
            }
        }

        return true;
    }

    /**
     * Initializes the buffer with records.
     */
    private void initBuffer()
        throws Exception
    {
        ListIterator<CloseableIterator<SAMRecord>> iter = null;
        int fileIndex = 0;

        this.buffer = new LinkedList<FlowAlignRecord>();

        iter = this.recordsIters.listIterator();
        while(iter.hasNext()) {
            CloseableIterator<SAMRecord> recordsIter = iter.next();
            if(recordsIter.hasNext()) {
                SAMRecord rec = recordsIter.next();
                this.addToBuffer(rec, fileIndex);
            }
            fileIndex++;
        }
    }

    /**
     * Adds the record in order genomic order to the buffer.
     * @param rec the SAM record to add.
     * @param the zero-based file index from which the record originated.
     */
    private void addToBuffer(SAMRecord rec, int fileIndex) 
        throws Exception
    {
        // search for correct position - implemented for simplicity
        // TODO: implement binary search for longer lists
        int i;
        SAMReadGroupRecord readGroup = null;
        FlowOrder flowOrder = null;
        for(i=0;i<this.buffer.size();i++) {
            FlowAlignRecord bufferRec = this.buffer.get(i);
            if(rec.getReferenceIndex() < bufferRec.record.getReferenceIndex() ||
                    (rec.getReferenceIndex() == bufferRec.record.getReferenceIndex() &&
                     rec.getAlignmentStart() <= bufferRec.record.getAlignmentStart())) 
            {
                break;
            }
        }

        try {
            readGroup = rec.getReadGroup();
        } catch(NullPointerException e) {
            throw new Exception("Read group ID not found in the SAM/BAM record");
        }
        if(null == this.flowOrders || null == this.flowOrders.get(fileIndex)) {
            throw new Exception("Could not find the flow order for the read with read group ID: " + readGroup.getId());
        }
        flowOrder = this.flowOrders.get(fileIndex).get(readGroup.getId());
        if(null == flowOrder) {
            throw new Exception("Could not find the flow order for the read with read group ID: " + readGroup.getId());
        }

        this.buffer.add(i, new FlowAlignRecord(rec, fileIndex, new FlowOrder(flowOrder)));
    }

    /**
     * @return true if the buffer is not empty, false otherwise.
     */
    public boolean hasNextAlignRecord()
    {
        if(0 == this.buffer.size()) {
            return false;
        }
        else {
            return true;
        }
    }

    /**
     * @return the next record in the buffer.
     */
    public FlowAlignRecord getNextAlignRecord()
        throws Exception
    {
        FlowAlignRecord ar = null;

        if(this.hasNextAlignRecord()) {
            ar = this.buffer.remove(0);
            CloseableIterator<SAMRecord> recordsIter = recordsIters.get(ar.fileIndex);
            if(recordsIter.hasNext()) {
                this.addToBuffer(recordsIter.next(), ar.fileIndex);
            }
        }
        return ar;
    }

    /**
     * Queries a range of SAM records in each input file, initializing the buffer with records from 
     * the new range.
     * @param sequenceName the sequence name, or null if we want the entire range.
     * @param  startPosition the start position (one-based).
     * @param endPosition the end position (one-based).
     */
    public void query(String sequenceName, int startPosition, int endPosition)
        throws Exception
    {
        ListIterator<SAMFileReader> readersIter = null;
        ListIterator<CloseableIterator<SAMRecord>> recordsItersIter = null;
        boolean indexed = true;

        if(null != sequenceName) {
            readersIter = this.readers.listIterator();
            recordsItersIter = this.recordsIters.listIterator();

            // Close all
            while(readersIter.hasNext()) {
                SAMFileReader reader = readersIter.next();
                if(reader.hasIndex()) {
                    CloseableIterator<SAMRecord> recordIter = recordsItersIter.next();
                    recordIter.close();
                    recordsItersIter.set(reader.query(sequenceName, startPosition, endPosition, false));
                }
                else {
                    indexed = false;
                }
            }
        }

        if(indexed || null == this.buffer) {
            this.initBuffer();
        }
    }

    /**
     * Adds the record for outputting.
     * @param rec the record to add to the output.
     */
    public void output(FlowAlignRecord rec)
    {
        if(1 == this.writers.size()) {
            if(this.headerCommentsAdded) {
                // Note: this does not check if we overwrite a previously set tag.
                rec.record.setAttribute(HEADER_FILE_INDEX_TAG, rec.fileIndex);
            }
            this.writers.get(0).addAlignment(rec.record);
        }
        else {
            this.writers.get(rec.fileIndex).addAlignment(rec.record);
        }
    }

    /**
     * Adds the unmapped record for outputting to the .
     * @param rec the record to add to the output.
     */
    public void outputUnmapped(FlowAlignRecord rec)
    {
        if(0 == this.writersUnmapped.size()) {
            // suppress
        }
        else if(1 == this.writersUnmapped.size()) {
            if(this.headerCommentsAdded) {
                // Note: this does not check if we overwrite a previously set tag.
                rec.record.setAttribute(HEADER_FILE_INDEX_TAG, rec.fileIndex);
            }
            this.writersUnmapped.get(0).addAlignment(rec.record);
        }
        else {
            this.writersUnmapped.get(rec.fileIndex).addAlignment(rec.record);
        }
    }

    /**
     * Closes all input(s) and output(s), and destroys the buffer.
     */
    public void closeAll()
    {
        int i;

        for(i=0;i<readers.size();i++) {
            readers.get(i).close();
        }
        for(i=0;i<writers.size();i++) {
            writers.get(i).close();
        }
        for(i=0;i<writersUnmapped.size();i++) {
            writersUnmapped.get(i).close();
        }
        for(i=0;i<recordsIters.size();i++) {
            recordsIters.get(i).close();
        }

        buffer = null;
    }
}
