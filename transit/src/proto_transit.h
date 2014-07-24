#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/transit.c */
/* FINDME: Main commented out to make MPI to work: */
//extern int main P_((int argc, char **argv));
extern void freemem_transit P_((struct transit *tr));

#undef P_
