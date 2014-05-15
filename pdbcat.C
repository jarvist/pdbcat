// simple PBD cat program
// Input:  either the fields or column version of a PDB file
// Output: either the   "     "   "       "     "  "  "   "

// Reason:  I want an easy way to use awk to change the contents
//  of the PDB file.  I could do lots of substr tricks, but I
//  would rather just convert the column based PDB file (in either
//  PDB or XPLOR formats) into a white-spaced seperated field
//  format.
//    The general trick should be something like:
// pdbcat -fields file.pdb | awk '$2=="ATOM" || $2=="HETATM" {
//  $11+=5; } { print $0; }' | pdbcat -columns > file2.pdb
// This adds 5 Angstroms to the y coordinate

// rearranged to prevent fclose(null) calls - Sergei Izrailev 5/28/98

#include <stdlib.h> // for exit
#include <stdio.h>  // for feof() and other file manip stuff
#include <string.h> // for strcasecmp
#include <iostream.h>
#include "Common.h"
#include "PDBData.h"



int main(int argc, char *argv[])
{
  Bool good = TRUE;
  Bool fields = FALSE; 	 // default is to use columns
  Bool usefiles = FALSE;       // default is to read from stdin
  
  int i;				// parse the command line
  for (i=1; i<argc; i++) {
    if (argv[i][0]!='-')  // then I must be done and the rest are files
      { usefiles = TRUE; break; }
    if (!strcasecmp(argv[i], "-f"))  // then all the rest are files
      { usefiles = TRUE; i++; break;}
    if (!strcasecmp(argv[i], "-help") || !strcasecmp(argv[i], "-h"))
      { good = FALSE; break; }
    if (!strcasecmp(argv[i], "-fields") || !strcasecmp(argv[i], "-field"))
      { fields = TRUE; continue; }
    if (!strcasecmp(argv[i], "-columns") || !strcasecmp(argv[i], "-column") ||
	!strcasecmp(argv[i], "-col"))
      { fields = FALSE; continue; }
    good = FALSE;
  }
  if (!good) {
    cerr << "Usage:  " << argv[0] << " {-fields | -columns} [[-f] files]" << '\n';
    cerr << "  Read any pdb file from stdin or list of files and convert the\n";
    cerr << " data to either a column based or field based pdb file.  A '#'\n";
    cerr << " represents an empty field.  This is useful for field based \n";
    cerr << " tools like awk.  The default format is 'columns'.\n";
    cerr << "\n    Written by Andrew Dalke <dalke@ks.uiuc.edu>\n";
    exit(1);
  }
  
  PDBData *pdb;
  char buf[151];
  char prnt[200];
  FILE *infile; // I don't know the "real" C++ way of doing this :(
  char *s;
  
  if (usefiles)  // are there files on the command line?
    if ( i < argc) { // are there any files left?
      infile = fopen(argv[i], "r");  // read from a file
    }
    else
      infile = NULL;  // no files left on the command line
  else
    infile = stdin;                // read from stdin
  
  while( (usefiles && i<argc) || (!usefiles && i==argc) ) {
    if (!infile)
      cerr << "Cannot open " << argv[i] << '\n';
    else {
      while ( fgets(buf, 150, infile) ) {
	for (s=buf; *s && *s!='\n'; s++)  // chop off the '\n'
	  ;
	if (*s=='\n')
	  *s = 0;
	else
	  cerr << "Input line is too long -- ignoring the rest of the line\n";
	
	pdb = new_PDBData(buf);
	if (pdb) {
	  if (fields)
	    pdb->sprint(prnt, PDBData::FIELDS);
	  else
	    pdb->sprint(prnt, PDBData::COLUMNS);
	  cout << prnt << '\n';
	  if (!usefiles)    // This may be interactive, so flush the
	    cout.flush();   //  output to be "friendlier".
	  delete pdb;
	} /* if */   
      } /* while */
      fclose(infile);
    } /* else */
    if(++i<argc) 
      infile = fopen(argv[i], "r");
  } /* while */

  return 0;
}
