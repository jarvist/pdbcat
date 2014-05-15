
// These are some type definitions used by the various atom programs

#ifndef COMMON_H
#define COMMON_H

enum Bool { FALSE, TRUE };

#ifndef NULL
#define NULL 0
#endif

// make a general definiton of floating point so it can be changed
// for different machines/ different floating point defs
typedef float Real;

// Print the string, print a "\n", and exit(-1)
void die(char *s);

// Print the string and a "\n", but keep on going
void warn(char *s);

#endif
