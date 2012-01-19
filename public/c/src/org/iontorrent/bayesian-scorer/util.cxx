/*
 * Copyright (C) 2010 Life Technologies Corporation. All rights reserved.
 */

#include "util.h"

bool __MemAlloc (void **ptr, long long size)
{
	if (size > 0)
	{
		*ptr = malloc(size);
		if (*ptr != NULL)
			return 1;
	}
	*ptr = NULL;
	return 0;
}


bool __MemResize (void **ptr, long long new_size)
{
	void * old_ptr = *ptr;

	if (old_ptr == NULL)
		return MemAlloc(ptr,new_size);

	if (new_size == 0)
		return MemDealloc(*ptr);

	void *new_ptr = realloc(old_ptr,new_size);
	if (new_ptr == NULL)
	{
		return 0;
	}

	*ptr = new_ptr;
	return 1;
}


bool __MemDealloc (void **ptr)
{
	if (*ptr != NULL)
	{
		free(*ptr);
		*ptr = NULL;
		return 1;
	}
	return 0;
}


#ifdef GENERATE_HTML

void PrintError (const char *message)
{
	fprintf(stdout,"<P><B>%s</B>\n",message);
}


void FatalError (const char *message)
{
	PrintError(message);
	exit(0);
}

#else

void PrintError (const char *message)
{
	fprintf(stderr,"%s\n",message);
}


void FatalError (const char *message)
{
	PrintError(message);
	exit(1);
}

#endif
