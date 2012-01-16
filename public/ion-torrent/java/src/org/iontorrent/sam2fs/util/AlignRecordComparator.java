/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2fs.util;

import java.io.*;
import java.util.*;
import net.sf.samtools.*;

/**
 * Compares align records based on their genomic co-oordinates.
 *
 * @author nils.homer@lifetech.com
 */
public class AlignRecordComparator implements Comparator<AlignRecord> 
{
    /**
     * The SAM record comparator.
     */
    private SAMRecordComparator comp = null;

    /**
     * Creates a new genomic co-ordinate comparator.
     */
    public AlignRecordComparator()
    {
        comp = new SAMRecordCoordinateComparator();
    }

    /**
     * Compares two align records.
     * @param o1 the first record.
     * @param o2 the second record.
     * @return a negative value if the co-ordinate of o1 is less than o2, zero if equal, a positive value otherwise.
     */
    public int compare(AlignRecord o1, AlignRecord o2)
    {
        return this.comp.compare(o1.record, o2.record);
    }

    /**
     * @param o1 the first record.
     * @param o2 the second record.
     * @return true if both have the same co-ordinate, false otherwise.
     */
    public boolean equals(AlignRecord o1, AlignRecord o2)
    {
        return (0 == this.comp.compare(o1.record, o2.record)) ? true : false;
    }
}
