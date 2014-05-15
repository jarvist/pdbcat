// Generic routines used by a lot of programs I write

#include <stdio.h>
#include <stdlib.h>
#include "Common.h"


// print an exit message and leave
void die(char *s)
{
  fflush(stdout);
  fflush(stderr);
  fprintf( stderr, "%s\nBye.\n", s);
  fflush(stderr);
  exit(-1);
}


// print a warning and a new line
void warn(char *s)
{
  fflush(stdout);
  fflush(stderr);
  fprintf(stderr, "%s\n", s);
  fflush(stderr);
}

// in case there is no system strcasecmp -- for some reason
// at least one version of gcc complained about not having one
#ifdef NOSTRCASECMP
int strcasecmp(const char *s1, const char *s2)
{
 char c1, c2;
 while (*s1 && *s2) {
   c1 = *s1;
   c2 = *s2;
   if (c1 >='a' && c1<='z')
      c1 -= 'a';
   if (c2 >='a' && c2<='z')
      c2 -= 'a';
   if (c1 != c2)
     return c2 - c1;
 }
  return 0;
}
#endif
