/*
 * Copyright (C) 2010 Life Technologies Corporation. All rights reserved.
 */

///////////////////////////////////////////////////////////////////
//
//		Electronic PCR (e-PCR) program
//
//		Gregory Schuler
//		Natonal Center for Biotechnology Information
//
//
//
//		Functions used by e-PCR to read FASTA files. 
//
///////////////////////////////////////////////////////////////////

#include "fasta-io.h"
#include "zutil.h"

char chrValidAa[] = "ABCDEFGHIKLMNPQRSTUVWXZ-*";
char chrValidNt[] = "GATCNBDHKMRSVWY-";  


///////////////////////////////////////////////////////////////
//
//
//		class FastaSeq
//
//

FastaSeq::FastaSeq ()
{
	m_tag = NULL;
	m_def = NULL;
	m_seq = NULL; 
	m_len = 0;
}


FastaSeq::FastaSeq (char *def, char *seq)
{
	m_tag = NULL;
	m_def = NULL;
	m_seq = NULL; 
	SetDefline(def);
	SetSequence(seq, strlen(seq));
}


FastaSeq::~FastaSeq ()
{
	Clear(); 
}


const char * FastaSeq::Label() const
{
	static char zippo[] = "XXXXXX";
	return (m_tag==NULL) ? zippo : m_tag; 
}


const char * FastaSeq::Defline () const
{
	static char zippo[] = "XXXXXX  No definition line found";
	return (m_def==NULL) ? zippo : m_def; 
}


const char * FastaSeq::Title () const
{
	static char zippo[] = "No definition line found";
	if (m_def==NULL) return zippo;
	char *p;
	if ((p = strchr(m_def, ' ')))
	{
		while (*p && isspace(*p))  p++;
		if (*p) return p;
	}
	return m_def; 
}


const char * FastaSeq::Sequence () const
{
	return m_seq; 
}


int FastaSeq::Length() const
{
	return m_len;
}


void FastaSeq::SetDefline (char *string)
{
	MemDealloc(m_tag);
	MemDealloc(m_def);
	m_def = string;
	if (string)
	{
		char *p = string;
		if (*p == ' ') p++;
		int n = strcspn(p," ");
		if (MemAlloc(m_tag,1+n))
		{
			memcpy((void*)m_tag,p,n);
			m_tag[n] = 0;
		}
	}
}


void FastaSeq::SetSequence (char *string, int len)
{
    //MemDealloc(m_seq); 
	m_seq = string; 
	m_len = string ? len : 0;
}


void FastaSeq::Clear()
{
	SetDefline(NULL); 
	SetSequence(NULL, 0); 
}


bool FastaSeq::ParseText (char *text, const char *alphabet, int len)
{
	if (text==NULL || *text==0)
		return 0;

	
	char *p2 = text;
	char *defline = NULL;
	if (*p2 == '>')
	{
		int n = strcspn((char*)text,"\r\n");
		p2 = text + n;
		n--;
		char *p = text+1;
		while (p<p2 && isspace(*p))
		{
			p++;
			n--;
		}
		if (n>0 && MemAlloc(defline,1+n))
		{
			memcpy((void*)defline,(void*)p,n);
			defline[n] =0;
		}
		
	}
	
	SetDefline(defline);
	*p2 = 0;
	SetSequence(p2+1, len-(p2+1-text));	
	
	return (m_len==0) ? 0 : 1;
}


bool FastaSeq::ParseText (char *text, int seqtype, int len)
{
	if (seqtype == SEQTYPE_AA)
		return ParseText(text,chrValidAa, len);

	if (seqtype == SEQTYPE_NT)
		return ParseText(text,chrValidNt, len);

	PrintError("ParseText(text,code);   ERROR: Invalid code");
	return 0;
}



///////////////////////////////////////////////////////////////
//
//
//		class FastaFile
//
//


FastaFile::FastaFile ()
{
	m_file = NULL;
	m_seqtype = 0;
	buf_total = 0;
}


FastaFile::FastaFile (int seqtype)
{
	m_file = NULL;
	m_seqtype = seqtype;
	buf_total = 0;
	buf = NULL;
	title[0] = 0;
}


FastaFile::~FastaFile ()
{
	if (m_file)  Close();
	MemDealloc(buf);
}


bool FastaFile::Open (const char *fname, const char *fmode)
{
	if (m_file != NULL)
	{
		PrintError("FastaFile::Open();  WARNING: file already open");
		return 0;
	}

	m_file = fopen(fname,fmode);
	if (m_file == NULL)
	{
		PrintError("FastaFile::Open();  ERROR: unable to open file");
		return 0;
	}

	return 1;
}


bool FastaFile::Close ()
{
	if (m_file == NULL)
	{
		PrintError("FastaFile::Close();  WARNING: file already closed");
		return 0;
	}

	fclose(m_file);
	m_file = NULL;
	return 1;
}


#define CHUNK 1000000

bool FastaFile::Read (FastaSeq &faseq)
{
	faseq.Clear();

	if (!IsOpen())
		return 0;

	if (title[0] == 0) {
	    do {
	        if (!fgets(title, sizeof title, m_file))  return 0;
		int len = strcspn(title, "\r\n");
		if (len == 0) continue;
		if (title[0] != '#' ) break;
	    } while (1);

	    //printf("%s\n", title);
	    if (title[0] != '>')
	    {
		PrintError("FastaFile::Read();  ERROR: was expecting '>'");
		return 0;
	    }
	}

	int len = strcspn(title, "\r\n")+1;
	if (buf_total == 0) {
	    buf_total = 1 + len + CHUNK;

	    if (!MemAlloc(buf,buf_total))
		return 0;
	}
	
	strcpy(buf,title);
	long long  buf_used = len;
	char *p = buf+len;
	title[0]= 0;

/*
#ifdef LINUX
	offset = ftello(m_file);
#else
	offset = ftello(m_file);
#endif
*/
	while (fgets(p, sizeof title, m_file))
	{
		if (p[0] == '#') continue;
		if (p[0] == '>')
		{
		    /*for (i = strlen(p)-1; i >=0; i--)
		        ungetc(p[i], m_file);
		    */
		    strcpy(title, p);
		    p[0] = '\0';
/*
#ifdef LINUX
		    fseeko(m_file,offset,SEEK_SET);
#else
		    fseek(m_file,offset,SEEK_SET);
#endif
*/
		    break;
		}
		len = strcspn(p, "\r\n");
		p[len] = 0;

		buf_used += len;
		p += len;
/*
#ifdef LINUX
		offset = ftello(m_file);
#else
		offset = ftello(m_file);
#endif
*/

		if (buf_used + (int) (sizeof(title))  > buf_total)
		{
			buf_total += buf_total;
			if (!MemResize(buf,buf_total))
			{
				PrintError("Out of memory");
				MemDealloc(buf);
				return 0;
			}
			p = buf + buf_used;
		}
	}

	bool result = faseq.ParseText(buf,m_seqtype, buf_used);

	//MemDealloc(buf);
	return result;
}

bool FastaFile::ReadHeader(FastaSeq &faseq)
{
	faseq.Clear();

	if (!IsOpen())
		return 0;

	char line[20000];
	while (fgets(line, sizeof line, m_file)) {
	    if (line[0] == '>') {

		int len = strlen(line);
		if (buf_total == 0) {
		    buf_total = 1 + len + CHUNK;
		    
		    if (!MemAlloc(buf,buf_total))
			return 0;
		}
		strcpy(buf,line);
		int buf_used = len;
		faseq.ParseText(buf,m_seqtype, buf_used);
		return 1;
	    }
	}
	return 0;
}

bool FastaFile::ReadBody(FastaSeq &faseq)
{
    strcpy(title, faseq.Defline());
    return Read(faseq);
}


bool FastaFile::Write (FastaSeq &seq)
{
	if (!IsOpen())
		return 0;

	fprintf(m_file,">%s\n", seq.Defline());
	WriteSeqLines(m_file,seq.Sequence(),seq.Length());
	return 0;
}


int WriteSeqLines (FILE *fd, const char *seq, int len, int linelen)
{
	char line[1+SEQLINE_LEN_MAX], *p1, ch;
	const char *p2 =seq;
	int i;

	///// force linelen to valid range
	if (linelen > SEQLINE_LEN_MAX)
		linelen = SEQLINE_LEN_MAX;
	else if (linelen < SEQLINE_LEN_MIN)
		linelen = SEQLINE_LEN_MIN;


	while (len >0)
	{
		p1 = line;
		for (i=0; i<len && i<linelen; ++i)
		{
			if ((ch = *p2++) ==0)
				break;
			*p1++ = ch;
		}
		*p1 = 0;
		if (fprintf(fd,"%s\n",line) <0)
			return 0;;
		len -= i;
	}
	return 1;
}
