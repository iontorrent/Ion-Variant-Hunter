/*
 * Copyright (C) 2010 Life Technologies Corporation. All rights reserved.
 */

#ifndef __util_h__
#define __util_h__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

/*
#ifndef bool
#define bool int
#warning defining bool as int
#endif
*/
#ifndef TRUE
#define TRUE true 
#endif

#ifndef FALSE
#define FALSE false
#endif 


bool __MemAlloc(void **ptr, long long size);
bool __MemResize(void **ptr, long long newsize);
bool __MemDealloc(void **ptr);

#define MemAlloc(x,y)     __MemAlloc((void**)&(x),(y))
#define MemResize(x,y)    __MemResize((void**)&(x),(y))
#define MemDealloc(x)     __MemDealloc((void**)&(x))

void PrintError(const char *message);
void FatalError(const char *message);

#endif /* __util_h__ */
