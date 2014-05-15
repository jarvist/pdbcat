/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile$
 *	$Author: dalke $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 94/10/05 17:05:06 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Read and parse a line of data from a PDB record.  There are many
 * different types of PDB records.  This version reads only the ATOM and
 * HETATM records and makes all fields accessible via the appropriate
 * member function.  In NAMD, this will be called only by the PDB class,
 * which reads PDB files, and a PDB writer class.
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log:	PDBData.h,v $
 * Revision 1.3  94/10/05  17:05:06  17:05:06  dalke (Andrew Dalke)
 * Changed residue length from 3 characters (PDB format) to 4
 * characters (XPLOR type PDB format)
 * 
 * Revision 1.2  94/08/12  15:12:10  15:12:10  nelson (Mark T. Nelson)
 * Moved destructor from protected to private
 * 
 * Revision 1.1  94/07/05  12:52:16  12:52:16  dalke (Andrew Dalke)
 * Initial revision
 * 
 ***************************************************************************/

#ifndef _PDBREADER_H_
#define _PDBREADER_H_

// These are added to the global namespace
//   the class PDBData
//   the class PDBUnknown, derived from PDBData
//   the class PDBAtom, derived from PDBData - contains ATOM and HETATM records
//   the class PDBAtomRecord, derived from PDBAtom - contains only ATOM records
//   the class PDBHetAtm, derived from PDBAtom - contains only HETATM records
//   the function new_PDBData, which creates the right pdb class given a string

#include <string.h>
#include "Common.h"

class PDBData {  // at the basic level, a PDB record only knows its type
   public:
// These data types come from the Protein Data Bank format
// description located via anon. ftp to pdb.pdb.bnl.gov
// in the file /pub/format.desc.ps
//  In addition, I define an extra type, UNKNOWN.  If I can't
// figure out what's going on, I store the complete string
// and return it when asked.
   enum PDBType {HEADER, OBSLTE, COMPND, SOURCE, EXPDTA, AUTHOR,
     REVDAT, SPRSDE, JRNL, REMARK, SEQRES, FTNOTE, HET, FORMUL,
     HELIX, SHEET, TURN, SSBOND, SITE, CRYST1, ORIGX, SCALE,
     MTRIX, TVECT, MODEL, ATOM, HETATM, SIGATM, ANISOU, SIGUIJ,
     TER, ENDMDL, CONECT, MASTER, END, UNKNOWN};

    static const char *PDBNames[UNKNOWN+1];   // string descriptors for each field

    enum PDBFormatStyle { COLUMNS, FIELDS };  // used to specify if the
       // output should be based on columns (FORTRAN style) or
       // fields (C/ awk style).
 // See, there are two different types of formats that this program
 // understands, one is the basic PDB w/ or w/o the XLPOR extension - these
 // are the column based versions.  The other is my own
 // field based version - each data element is seperated by a blank
 // and, if the element is empty, a pound sign ('#') is put in its place.
 // This type of record is denoted by a '#' in the first non-blank
 // character (hence, it is the first non-blank character of the first
 // field.  Basically, I'm a unix/ C/ awk/ yacc ... freak, I like field
 // based data rather than column based data.
 
 
   private:
    PDBType mytype;

   protected:
        // some parsing routines to get info from a line of text
    static void scan( const char *data, int len, int start, int size, 
                         int *ans, int defalt);
    static void scan( const char *data, int len, int start, int size,
                          Real *ans, Real defalt);
    static void scan( const char *data, int len, int start, int size,
                         char *ans);
    static void field( const char *data, int fld, char *result);
        // some routine to print to a specific column and width
    static void sprintcol( char *s, int start, int len, const char *val);
    static void sprintcol( char *s, int start, int len, int val);
    static void sprintcol( char *s, int start, int len, int prec, Real val);
    
   public:
     PDBData(PDBType newtype) {
       mytype = newtype;
     }
     virtual ~PDBData( void) {
     }
     PDBType type( void) {
       return mytype;
     }
                    // I know nothing, so I'll fake it and hope it works
     virtual void sprint( char *s, PDBFormatStyle usestyle = COLUMNS) {
       if (usestyle == COLUMNS)     // get rid of warning
         strcpy(s, "REMARK     (undefined remark - this is a bug)");
        else
         strcpy(s, "REMARK     (undefined remark - this is a bug)");
     }
};

//////******* the "UNKNOWN" class *****//////
class PDBUnknown : public PDBData {
  private:
    char *mystr;
  public:
   PDBUnknown(const char *data): PDBData(PDBData::UNKNOWN) {
     mystr = new char[strlen(data)+1];
     strcpy(mystr, data);
   }
   virtual ~PDBUnknown( void) {
     delete [] mystr;
   }  
   void sprint(char *s, PDBFormatStyle usestyle) {
     strcpy(s, mystr);
     if (usestyle == PDBData::COLUMNS)   // they are the same, but I won't
       strcpy( s, mystr);                //   get the stupid warning during
      else                               //   compilation
       strcpy( s, mystr);
   }
};

////************* routines used for ATOM and HETATM **********/////
class PDBAtom : public PDBData {
  private:
      // starting location for each record element
    enum Start {STYPE=1,SSERIAL=7, SNAME=13, SALT=17, SRESNAME=18, SCHAIN=22, 
                SRESSEQ=23, SINSERT=27, SX=31, SY=39, SZ=47,
                SOCC=55, STEMPF=61, SFOOT=68, SSEGNAME=73};
      // length of each element, the PREC is the number of digits
      // in the output after the decimal
// NOTE: The PDB says the length of the residue name is only 3 characters
//  whereas XPLOR allows 4 character names.  We choose 4 for compatability
//  with both systems (since we never change the length, we you give us is
//  what we use)
    enum Length {LTYPE=6, LSERIAL=5, LNAME=4, LALT=1, LRESNAME=4, LCHAIN=1, 
                 LRESSEQ=4, LINSERT=1, LCOOR=8,
                 LCOORPREC=3, LOCC=6, LOCCPREC=2, LTEMPF=6, 
                 LTEMPFPREC=2, LFOOT=3, LSEGNAME=4};

    static const int default_serial;         // some default values
    static const int default_residueseq;     // these are set in the .C file
    static const Real default_coor;
    static const Real default_occupancy;
    static const Real default_temperaturefactor;
    static const int no_footnote;

    int myserialnumber;                 // atom serial number
    char myname[LNAME+1];               // atom name
    char myalternatelocation[LALT+1];   // alternamte location identifier
    char myresiduename[LNAME+1];        // residue name
    char mychain[LCHAIN+1];             // chain indentifier
    int myresidueseq;                   // residue seq. no.
    char myinsertioncode[LINSERT+1];    // code for insertions of residues
    Real mycoor[3];                     // X, Y, and Z orthogonal A coordinates
    Real myoccupancy;                   // occupancy
    Real mytemperaturefactor;           // temperature factor
    int myfootnote;                     // footnote number
    char mysegmentname[LSEGNAME+1];     // XPLOR-type segment name

    void parse_field_data( const char *data);
    void parse_column_data( const char *data);
    void sprint_columns( char *outstr);
    void sprint_fields( char *outstr);

  protected:
    enum PDBPossibleAtoms {USE_ATOM = ATOM, USE_HETATM = HETATM};
    PDBAtom( const char *data,
           PDBPossibleAtoms whichatom);// parses a line from the PDB data file
    PDBAtom( void);        // makes a generic atom

  public:
    virtual ~PDBAtom( void);
    void parse( const char *s);  // reset to new input values
    void  sprint( char *s, PDBFormatStyle usestyle = COLUMNS);// write to string
    int serialnumber( void);
    void serialnumber( int newserialnumber);
    
    const char*name( void);
    void name( const char *newname);
    
    const char*alternatelocation( void);
    void alternatelocation( const char *newalternatelocation);
    
    const char*residuename( void);
    void residuename( const char *newresiduename);
    
    const char*chain( void);
    void chain( const char *newchain);
    
    int residueseq( void);
    void residueseq( int newresidueseq);
    
    const char*insertioncode( void);
    void insertioncode( const char *newinsertioncode);
    
    Real xcoor( void);
    void xcoor( Real newxcoor);
    Real ycoor( void);
    void ycoor( Real newycoor); 
    Real zcoor( void);
    void zcoor( Real newzcoor);
    
    const Real *coordinates( void);
    void coordinates(const Real *newcoordinates);
    
    Real occupancy( void);
    void occupancy( Real newoccupancy);

    Real temperaturefactor( void);
    void temperaturefactor( Real newtemperaturefactor);

    int footnote( void);
    void footnote( int newfootnote);
    
      // this is not part of the PDB format but is used by XPLOR instead of
      // the chain identifier (see XPLOR 3.1 manual, p 104)
    const char*segmentname( void);
    void segmentname( const char *newsegmentname);
};

// The two sub-classes of PDB Atom
class PDBAtomRecord : public PDBAtom{
   public:
     PDBAtomRecord( const char *data ) :
          PDBAtom( data, PDBAtom::USE_ATOM) {
     }
     virtual ~PDBAtomRecord( void) {
     }
};

class PDBHetatm : public PDBAtom {
  public:
    PDBHetatm( const char *data) :
         PDBAtom( data, PDBAtom::USE_HETATM) {
    }
    virtual ~PDBHetatm( void) {
    }
};


////********* Wrap up everything in one function call ***********//////
// somehow I need the base class to figure out which derived class
// to use to parse.   Since I don't know how to do that, I'll
// fake it with this.  Give it a string and it will create the
// correct PDB data type.
PDBData *new_PDBData(const char *data);  // nasty


#endif
