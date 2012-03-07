/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2flowgram.util;

import java.util.*;

/**
 * Used to create align record queus or iterators from a single input queue.
 *
 * @author nils.homer@lifetech.com
 */
public class ThreadPoolLinkedList {

    // a buffer of align records
    private LinkedList<FlowAlignRecord> buffer; // internal use only

    /**
     * Initializes an empty buffer.
     */
    public ThreadPoolLinkedList()
    {
        this.buffer = new LinkedList<FlowAlignRecord>();
    }

    /**
     * @param e the record to add to the buffer.
     */
    public void add(FlowAlignRecord e)
    {
        this.buffer.add(e);
    }

    /**
     * @return the size of the buffer.
     */
    public int size() 
    {
        return this.buffer.size();
    }

    /**
     * @return the first record in the buffer, null if empty.
     */
    public FlowAlignRecord getFirst()
    {
        if(0 == this.size()) {
            return null;
        }
        return this.buffer.getFirst();
    }

    /**
     * @return the last record in the buffer, null if empty.
     */
    public FlowAlignRecord getLast()
    {
        if(0 == this.size()) {
            return null;
        }
        return this.buffer.getLast();
    }

    /**
     * @return removes and returns the first record in the buffer, null if empty.
     */
    public FlowAlignRecord removeFirst()
    {

        if(0 == this.size()) {
            return null;
        }
        return this.buffer.removeFirst();
    }

    /**
     * Creates a list of list of records for the given number of threads from the buffer.  This
     * will only return records with the matching reference index.
     * @param numThreads the number of threads.
     * @param referenceIndex the reference index of records to return, or -1 if all are to be returned.
     * @return a list of list of records, one for each thread.
     */
    public LinkedList<LinkedList<FlowAlignRecord>> getThreadLists(int numThreads, int referenceIndex)
    {
        int i = 0;
        int size = this.size();
        LinkedList<LinkedList<FlowAlignRecord>> threadLists = new LinkedList<LinkedList<FlowAlignRecord>>();

        for(i=0;i<numThreads;i++) {
            threadLists.add(new LinkedList<FlowAlignRecord>());
        }

        i=0;
        while(0 != size) {
            FlowAlignRecord rec = this.buffer.getFirst();
            if(0 <= referenceIndex && rec.record.getReferenceIndex()+1 != referenceIndex) {
                break;
            }
            threadLists.get(i).add(this.buffer.removeFirst());
            size--;
            i++;
            if(numThreads <= i) {
                i=0;
            }
        }

        return threadLists;
    }
    
    /**
     * Creates a list of list of records for the given number of threads from the buffer.  This
     * will only return records with the matching reference index and alignment end less than the
     * given upper bound.  Do not use this method if we are to flush all records.
     * @param numThreads the number of threads.
     * @param referenceIndex the reference index of records to return.
     * @param alignmentStartUpperBound the maximum alignment start position of any record to return.
     * @return a list of list of records, one for each thread, null if there were no entries found.
     */
    public LinkedList<LinkedList<FlowAlignRecord>> getThreadListsStartBound(int numThreads, int referenceIndex, int alignmentStartUpperBound)
    {
        int i = 0;
        int size, numAdded;
        LinkedList<LinkedList<FlowAlignRecord>> threadLists = null;

        size = this.size();
        if(size == 0) {
            return null;
        }
            
        threadLists = new LinkedList<LinkedList<FlowAlignRecord>>();

        for(i=0;i<numThreads;i++) {
            threadLists.add(new LinkedList<FlowAlignRecord>());
        }

        i = numAdded = 0;
        while(0 != size) {
            FlowAlignRecord rec = this.buffer.getFirst();
            if(rec.record.getReferenceIndex()+1 != referenceIndex ||
                    alignmentStartUpperBound <= rec.record.getAlignmentStart()) {
                break;
            }
            threadLists.get(i).add(this.buffer.removeFirst());
            size--;
            numAdded++;
            i++;
            if(numThreads <= i) {
                i=0;
            }
        }

        if(0 == numAdded) {
            return null;
        }
        else {
            return threadLists;
        }
    }

    /**
     * Creates a list of list of records for the given number of threads from the buffer.  This
     * will only return records with the matching reference index and alignment end less than the
     * given upper bound.  Do not use this method if we are to flush all records.
     * @param numThreads the number of threads.
     * @param referenceIndex the reference index of records to return.
     * @param alignmentEndUpperBound the maximum alignment end position of any record to return.
     * @return a list of list of records, one for each thread, null if there were no entries found.
     */
    public LinkedList<LinkedList<FlowAlignRecord>> getThreadListsEndBound(int numThreads, int referenceIndex, int alignmentEndUpperBound)
    {
        int i = 0;
        int size, numAdded;
        LinkedList<LinkedList<FlowAlignRecord>> threadLists = null;

        size = this.size();
        if(size == 0) {
            return null;
        }
            
        threadLists = new LinkedList<LinkedList<FlowAlignRecord>>();

        for(i=0;i<numThreads;i++) {
            threadLists.add(new LinkedList<FlowAlignRecord>());
        }

        i = numAdded = 0;
        while(0 != size) {
            FlowAlignRecord rec = this.buffer.getFirst();
            if(rec.record.getReferenceIndex()+1 != referenceIndex ||
                    alignmentEndUpperBound <= rec.record.getAlignmentEnd()) {
                break;
            }
            threadLists.get(i).add(this.buffer.removeFirst());
            size--;
            numAdded++;
            i++;
            if(numThreads <= i) {
                i=0;
            }
        }

        if(0 == numAdded) {
            return null;
        }
        else {
            return threadLists;
        }
    }
}
