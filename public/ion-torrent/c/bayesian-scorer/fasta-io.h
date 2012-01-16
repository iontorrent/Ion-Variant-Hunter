/*
 * Copyright (C) 2010 Life Technologies Corporation. All rights reserved.
 */

#ifndef __fasta_io_h__
#define __fasta_io_h__

#include "util.h"

#define SEQLINE_LEN_DEFAULT 60   // default sequence characters per line
#define SEQLINE_LEN_MAX    120   // max sequence characters per line
#define SEQLINE_LEN_MIN     10   // min sequence characters per line

extern char chrValidAa[];
extern char chrValidNt[];

#define SEQTYPE_AA 1
#define SEQTYPE_NT 2


class FastaSeq
{
public:
	FastaSeq ();
	FastaSeq (char *def, char *seq);
	~FastaSeq ();

	const char * Label() const;
	const char * Title() const;
	const char * Defline() const;
	const char * Sequence() const;
	int Length() const;

	void Clear();

	bool  ParseText(char *text, const char *alphabet, int len);
	bool  ParseText(char *text, int seqtype, int len);

protected:
	char *m_tag;
	char *m_def;
	char *m_seq;
	int   m_len;
	void *m_bogus;

	void SetDefline (char *string);
	void SetSequence (char *string, int);

	friend class FastaFile;
};


class FastaFile 
{
public:
	FastaFile();
	FastaFile(int seqtype);
	~FastaFile();

	bool Open (const char *fname, const char *fmode);
	bool Close ();

	bool IsOpen () const
		{ return (m_file ==NULL) ? 0 : 1; }

	bool Read (FastaSeq &seq);
	bool Write (FastaSeq &seq);
	bool ReadHeader(FastaSeq &seq);
	bool ReadBody(FastaSeq &seq);

protected:
	FILE *m_file;
	int m_seqtype;
	long long  buf_total;
	char *buf;
	char title[20000];
};


extern int WriteSeqLines (FILE *fd, const char *seq, int len, int linelen=SEQLINE_LEN_DEFAULT);


#endif
