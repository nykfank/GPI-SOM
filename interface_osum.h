#define MAXSEQ 40000                    // max number of seqs that can be loaded
#define CLN 32                          // length of the C terminal sequence
#define LOCI 22 	                // number of loci to be looked at
#define SLN 100000
#define PLEN 8         // max omega pattern length
#define PNUM 3          // number of omega patterns
#define INP LOCI+20+2   	          // input neurons = loci+AAzentiole+omega
#define WORDSIZE 20
int load_fasta(char* fn,int start);     // loads fasta in seq and desc array
void interface(int k);                  // computes input for NN from seqs in x array
void omega_search(char sq[100000],int sn);
void ohm(char os[SLN],int osp,int scor,int oi,int moi);
char *seq[MAXSEQ],*desc[MAXSEQ];	// string and string arrays
float x[MAXSEQ][INP];                   // 2D array with NN inputs for all seqs
int ome[MAXSEQ];			// array of omega values
int ope[MAXSEQ];			// array of omega quality
int numseq;                             // number of seqs read
int mx;long mxp;

void interface(int k) {
	// input: LOCI hydro - 20 Zentiole - 2 omega - 1 pkr - 1 mass
	int pol['Z'];int pos['Z'];long sp,i;
	int loc[LOCI]={2,9,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
	// kyte & doolittle hydrophobicity grouping, data from M. Solioz
	pol['G']=-24;pol['A']=-19;pol['V']=-20;pol['L']=-23;pol['I']=-22;
	pol['C']=-12;pol['M']=-15;pol['F']=8;pol['Y']=61;pol['W']=59;
	pol['P']=60;pol['S']=51;pol['T']=49;pol['N']=97;pol['Q']=94;
	pol['D']=110;pol['E']=102;pol['H']=103;pol['K']=150;pol['R']=200;
	// positions of AAs in input field
	pos['A']=0;pos['R']=1;pos['N']=2;pos['D']=3;pos['C']=4;pos['Q']=5;
	pos['E']=6;pos['G']=7;pos['H']=8;pos['I']=9;pos['L']=10;pos['K']=11;
	pos['M']=12;pos['F']=13;pos['P']=14;pos['S']=15;pos['T']=16;
	pos['W']=17;pos['Y']=18;pos['V']=19;
	pos['X']=0;pos['Z']=0;pos['O']=0;pos['J']=0;pos['B']=0;
	printf(".");
	sp=strlen(seq[k])-CLN;  // calculate beginning of C terminus
	x[k][i+LOCI+22]=0;x[k][i+LOCI+23]=0;
	if (sp<1) {printf("WARNING: seq%d < %dAAs.\n",k,CLN);return;}
	for (i=LOCI;i<LOCI+20;i++) {x[k][i]=0;} // reset pos value
	for (i=0;i<LOCI;i++) {          // loop through loci
	// calc normalized hydrophobicity of C terminus
		x[k][i]=(pol[seq[k][sp+loc[i]]]+25)/(float)225;
	// calc mean positions AAs
		x[k][LOCI+pos[seq[k][sp+loc[i]]]]=(x[k][LOCI+pos[seq[k][sp+loc[i]]]]+loc[i])/2;
	}
	// normalize mean positions was 20
	for (i=LOCI;i<LOCI+20;i++) {x[k][i]/=20;}
	omega_search(seq[k],k);
}

void omega_search(char sq[SLN],int sn) {
        long k;float perc,nopo;
        mx=0;mxp=0;
        for (k=strlen(sq)-CLN;k<strlen(sq);k++) {ohm(sq,k,0,0,PNUM-1);}
        mxp-=strlen(sq)-32;if (mx==0) {mxp=0;}
        perc=(float)mx/12;nopo=(float)mxp/32;
	ome[sn]=32-mxp;ope[sn]=(int)((float)mx/12*100);
        x[sn][LOCI+20]=nopo;x[sn][LOCI+21]=perc;
}

void ohm(char os[SLN],int osp,int scor,int oi,int moi) {
	char w[PNUM][PLEN]={
        "SNDGACXX",
        "ASDTRCMW",
        "AGSTVXXX"};
        int wp[PNUM][PLEN]={
        4,3,2,1,1,1,0,0,
        4,3,2,1,1,1,1,1,
        4,3,2,1,1,0,0,0};
        int ps[PNUM]={0,1,2};int ok;
        for (ok=0;ok<PLEN;ok++) {
                if (os[osp+ps[oi]]==w[oi][ok]) {
                        scor+=wp[oi][ok];
                        if (oi<moi) {ohm(os,osp,scor,oi+1,moi);}
                        else {if (scor>mx) {mx=scor;mxp=osp;}}
                }
        }
}

int load_fasta(char* fn,int n) {
        FILE *ifp;int once=0;char line[SLN];
        if ((ifp=fopen(fn,"r"))==NULL) {printf("%s not found\n",fn);return(-1);}
        numseq=0; // n is pos in the seq/desc array, numseq a counter
        while (!feof(ifp)) {  // while there are lines..
                strcpy(line,"");fgets(line,sizeof(line),ifp);line[strlen(line)-1]='\0';
                if (line[0]=='>') {  // if its a desc
                        if (once>1) {printf("One seq per line!\n");return(-1);}
                        n++;once=0;numseq++; // inc counters
                        desc[n]=(char*)malloc(sizeof(char)*(strlen(line)+1));
                        strcpy(desc[n],line);   // copy desc into desc array
                }
                else {
                        if (strlen(line)>1) { // because of the final newline
                                seq[n]=(char*)malloc(sizeof(char)*(strlen(line)+1));
                                strcpy(seq[n],line);once++;
                        }
                }
        }
        fclose(ifp);printf ("- %s: %d sequence(s)\n",fn,numseq);return(0);
}
