struct mafAli
/* A multiple alignment. */
{
  struct mafAli *next;
  double score;
  struct mafComp *components;	/* List of components of alignment */
  int textSize;	 /* Size of text in each component. */
  int chain_len;
  int label;
  char orient; /* Relative orientation of the reference */
  
};
struct mafComp
/* A component of a multiple alignment. */
    {
      struct mafComp *next;
      char *name;        /* comman name of sequence source. */
      char *src;	 /* Name of sequence source.  */
      char *text;        /* The sequence including dashes. */
      char* contig;
      int* mafPosMap;
      int srcSize;       /* Size of sequence source.  */
      int start;	 /* Start within sequence. Zero based. If strand is - is relative to src end. */
      int size;	         /* Size in sequence (does not include dashes).  */
      short nameID;
      char strand;       /* Strand of sequence.  Either + or -*/
      char paralog;
};

extern "C" void parseSrcName(char* srcName, char* name, char* src);
extern "C" struct mafFile *mafOpen(const char *fileName, int verbose);
extern "C" struct mafAli *mafNext(struct mafFile *mafFile);
extern "C" void mafWrite(FILE *f, struct mafAli *maf);
extern "C" void mafWriteStart(FILE *f, char *scoring);
extern "C" void mafFileFree(struct mafFile **pObj);
extern "C" void mafAliFree(struct mafAli **pObj);
