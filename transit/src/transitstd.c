// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <transit.h>

/* keeps tracks of number of errors that where allowed to continue. */
int verblevel;
int maxline=1000;

inline void transitdot(int thislevel,
                       int verblevel,
                       ...){
  if(thislevel <= verblevel)
    fwrite(".", 1, 1, stderr);
}

void
tr_output_fcn (int flags,
               const char *file,
               const long line,
               const char *str,
               ...) {

  va_list format;
  va_start(format, str);
  tr_output_vfcn(flags, file, line, str, format);
  va_end(format);
}

void
tr_output_vfcn (int flags,
                const char *file,
                const long line,
                const char *str,
                va_list format) {

  // Obtain the level of output from the flag data.
  int level = flags & TOUT_VERBMASK;

  // Choose output location: stdout or stderr.
  FILE *output = (level == TOUT_ERROR) ? stderr : stdout;

  // If the caller has chosen to use a banner, print the preceeding line.
  if (flags & TOUT_BANNER)
    fprintf(output, "\n--------------------------------------------------\n");

  /* Always print debugging information (file and line number) for errors and
   * warnings. Selectively print this information for other messages based on
   * the TOUT_LOCATE flag.
   */
  if (level == TOUT_ERROR)
    fprintf(output, "Transit ERROR :: (%s, line %lu)\n", file, line);

  else if (level == TOUT_WARN)
    fprintf(output, "Transit WARNING :: (%s, line %lu)\n", file, line);

  else if (flags & TOUT_LOCATE) {

    switch(level) {
      case TOUT_INFO:
        fprintf(output, "Transit INFO :: (%s, line %lu)\n", file, line);
        break;

      case TOUT_RESULT:
        fprintf(output, "Transit RESULT :: (%s, line %lu)\n", file, line);
        break;

      case TOUT_DEBUG:
        fprintf(output, "Transit DEBUG :: (%s, line %lu)\n", file, line);
        break;

      default:
        break;
    }
  }

  // Print the output message itself.
  vfprintf(output, str, format);

  // If the caller has chosen to use a banner, print the succeeding line.
  if (flags & TOUT_BANNER)
    fprintf(output, "--------------------------------------------------\n");
}


/*\fcnfh
  Check whether the file 'in' exists and is openable. If so, return
  the opened file pointer 'fp'. Otherwise a NULL is returned and a
  status of why the opening failed is returned.

  Return: 1 on success in opening
          0 if no file was given
         -1 File doesn't exist
         -2 File is not of a valid kind (it is a dir or device)
         -3 File is not openable (permissions?)
         -4 Some error happened, stat returned -1                           */
int
fileexistopen(char *in,    /* Input filename                                */
              FILE **fp){  /* Opened file pointer if successful             */
  struct stat st;

  /* Return if no file was given:                                           */
  if (!in)
    return 0;

  /* Check if the suggested file exists. If it doesn't, use defaults:       */
  if (stat(in, &st) == -1){
    if (errno == ENOENT)
      return -1;
    else
      return -4;
  }

  /* Check if suggested file is of a valid type:                            */
  if (!(S_ISREG(st.st_mode) || S_ISFIFO(st.st_mode)))
    return -2;

  /* If no file pointer variable is given, we've done all we can:           */
  if (fp == NULL)
    return 1;

  /* If a file pointer is given, attempt to open the file:                  */
  *fp = NULL;

  if (((*fp)=fopen(in,"r")) == NULL)
    return -3;

  /* Success:                                                               */
  return 1;
}


/*\fcnfh
  Output for the different cases. of fileexistopen()

  @return fp of opened file on success
          NULL on error (doesn't always returns though */
FILE *
verbfileopen(char *in,     /* Input filename              */
             char *desc){  /* Comment on the kind of file */
  FILE *fp;

  switch(fileexistopen(in, &fp)){
  /* Success in opening or user don't want to use atmosphere file: */
  case 1:
    return fp;
  case 0:
    tr_output(TOUT_ERROR, "No file was given to open.\n");
    return NULL;
  /* File doesn't exist: */
  case -1:
    tr_output(TOUT_ERROR, "%s file '%s' doesn't exist.\n", desc, in);
    return NULL;
  /* Filetype not valid: */
  case -2:
    tr_output(TOUT_ERROR, "%s file '%s' is not of a valid kind "
                               "(it is a dir or device)\n", desc, in);
    return NULL;
  /* File not openable: */
  case -3:
    tr_output(TOUT_ERROR, "%s file '%s' is not openable. Probably "
                               "because of permissions.\n", desc, in);
    return NULL;
  /* stat returned -1: */
  case -4:
    tr_output(TOUT_ERROR,
                 "Error happened for %s file '%s', stat() returned -1, "
                 "but file exists.\n", desc, in);
    return NULL;
  default:
    tr_output(TOUT_ERROR,
                 "Something weird in file %s, line %i.\n", __FILE__, __LINE__);
  }
  return NULL;
}

/*
    Check that the 'n' functions 'variable_argument' have been called
    before 'fcn'. If any of them have not been yet called, throw an error
    message detailed all functions not yet called.                          */
void
transitcheckcalled(const long pi,   /* Progress indicator variable          */
                   const char *fcn, /* Name of function being checked       */
                   const int n,     /* Number of functions that have to be
                                       called before fcn                    */
                   ...){            /* Function's name and flag pairs       */
  /* TD: make 'mess' dynamic size */
  va_list ap;      /* Variable-argument pointer        */
  int i;           /* Auxiliary for index              */
  char mess[1000], /* Output message                   */
       *name;      /* Function's name                  */
  long flag;       /* Function's flag                  */
  _Bool stop=0;    /* Any function has not been called */

  *mess = '\0';
  va_start(ap, n);
  /* Check the pair of arguments: */
  for(i=0; i<n; i++){
    name = (char *)va_arg(ap, char *);
    flag = (long)va_arg(ap, long);
    /* Check if it was called before: */
    if(!(pi&flag)){
      /* Append text to mess: */
      if(!stop){ /* Append only once at first function not called: */
        strcpy(mess, "The following function(s) were not executed before "
                     "this execution of '%s()':\n");
      }
      strcat(mess, "  ");
      strcat(mess, name);
      strcat(mess, "()\n");
      stop = 1;
    }
  }
  va_end(ap);
  /* Print out the error: */
  if(stop) {
    tr_output(TOUT_ERROR, mess, fcn);
    exit(EXIT_FAILURE);
  }
}


void
freemem_molecules(struct molecules *mol, long *pi){
  /* Free structures: */
  for(int i=0; i<mol->nmol;  i++)
    free_mol(mol->molec+i);
  /* Free arrays:     */
  free(mol->name[0]);
  free(mol->name);
  free(mol->molec);
  free(mol->mass);
  free(mol->radius);
  free(mol->pol);
  free(mol->ID);
  /* FINDME: Define a pi for molec */
}

/* \fcnfh
   Frees array in prop_isov, this should be called only once for all the
   isotopes.                                                                */
void
free_isov(prop_isov *isov){
  free(isov->z);
}


/* \fcnfh
   Frees array in prop_isof, this should be called for each of the
   isotopes.                                                       */
void
free_isof(prop_isof *isof){
  free(isof->n);
}


void
free_mol(prop_mol *molec){
  free(molec->d);
  free(molec->q);
}


/* \fcnfh
   Frees array in prop_db, this should be called for each of the
   isotopes.                                                       */
void
free_db(prop_db *db){
  free(db->n);
}


/* \fcnfh
   Frees array in prop_dbnoext, this should be called once per each
   database.                                                       */
void
free_dbnoext(prop_dbnoext *db){
  free(db->T);
}

/* \fcnfh
   Frees array in prop_samp */
void
free_samp(prop_samp *samp){
  free(samp->v);
}


/* \fcnfh
   Frees array in prop_atm */
void
free_atm(prop_atm *atm){
  free(atm->p);
  free(atm->t);
  free(atm->mm);
}

/* \fcnfh
   Saves a string in binary in open file */
void
savestr(FILE *out,
        char *str){
  long len = strlen(str)+1;

  fwrite(&len, sizeof(long), 1, out);
  fwrite(str, 1, len, out);
}

/* \fcnfh
   Restores a string from a binary open file */
int
reststr(FILE *in,
        char **str){
  long len;

  if(fwrite(&len, sizeof(long), 1, in) != 1)
    return -1;
  if(len<0)
    return -2;
  if(len>1000)
    return 1;
  if((*str=(char *)calloc(len, sizeof(char)))==NULL)
    return -3;
  if(fread(*str, 1, len, in) != len)
    return -1;

  return 0;
}


/* \fcnfh
   This function is called if a line of 'file' was longer than 'max'
   characters
*/
void
linetoolong(int max,     /* Maxiumum length of an accepted line */
            char *file,  /* File from which we were reading     */
            int line){   /* Line who was being read             */
  tr_output(TOUT_ERROR,
    "Line %i of file '%s' has more than %i characters, "
    "that is not allowed.\n", file, max);
  exit(EXIT_FAILURE);
}


double
timestart(struct timeval tv,  /* timeval structure                          */
          char *str){         /* Time stamp                                 */
  /* Get current time:                                                      */
  gettimeofday(&tv, NULL);
  /* Calculate time in seconds:                                             */
  double sec = tv.tv_sec + 1e-6*tv.tv_usec;
  tr_output(TOUT_INFO, "%s\n", str);
  return sec;
}

/* Print to screen the time elapsed since time t0
   Return: current time in seconds                       */
double
timecheck(int verblevel,      /* Verbosity level         */
          long iter,          /* Iteration index         */
          long index,         /* Sequencial index        */
          char *str,          /* Time stamp description  */
          struct timeval tv,  /* timeval structure       */
          double t0){         /* Time in seconds         */
  /* Get current time:          */
  gettimeofday(&tv, NULL);
  /* Calculate time in seconds: */
  double sec = tv.tv_sec + 1e-6*tv.tv_usec;
  /* Print time stamp:          */
  tr_output(TOUT_RESULT, "Check point: %02li - %02li %s:  dt = %.4f "
    "sec.\n\n", iter, index, str, sec-t0);
  return sec;
}
