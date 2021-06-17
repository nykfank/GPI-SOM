/*
KohGPI v1.5 by Nick Fankhauser. Identifies GPI anchored proteins by C terminus.
v1.0 written 9.7.2003
v1.1 update 15.8.2003: Map also made from training set, error percent display.
v1.2 update 17.7.2003:	no more maybe.
			undescidable mapped as undef.
			improved radial search.
			train/plot: DIM+1 -> =>DIM. ga,nga size.
v1.3 update 19.2.2004:	support for unequal training sets
v1.4 update 2.3.2004:	Add [unmapGPI:NOGPI] to unmapped seqs when evaluating
v1.5 update 3.3.2004:	replaced radial search by radius2 quadrant search
*/
#include "gd.h"		/* C library for the dynamic creation of images */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
/* /System/Library/Frameworks/Kernel.framework/Versions/A/Headers/sys/ */
#include <sprlib.h>	/* StatisticalPatterRecognitionArtificialNeuralNetworkLibrary */
#include "interface_osum.h"	/* contains interface and load_fasta functions */
#define MAXCYC 10000	/* number of cycles to train the NN with the training set */
#define DIM1 40		/* first dimension of the kohonen self organizing map */
#define DIM2 40		/* second dimension of the kohonen self organizing map */
#define LEARN 0.001	/* the learning rate */
#define PIX 10		/* size of a pixel in the NN visualisation */
#define PIZ 5		/* half of pix */
#define ZD 2		/* distance of cross from pixel border */
#define MINC 100	/* minimal color brightness */
#define INFILE "input.txt"	/* input seq file for the web daemon */
void summon_dataset(void);	/* converts the input vector into sprannlib format */
void train(void);		/* trains the NN */
void search(void);		/* evaluates input sequences */
void parse(char uname[80],int *xp,int *yp);	/* get coordinates from unit desc */
void plot(long c);	/* creates a PNG image of the NN */
void resplot(void);	/* creates a PNG image of results */
void webdaemon(void);	/* loads the NN&map and then waits for INFILE to be created by CGI */
void mapnet(void);	/* loads the map and the net */
void interactive(void);	/* loads map&net and then asks for filenames to evaluate */
void outfiles(char filename[100]);		/* adds extensions for output files */
DATASET *set;				/* input vectors in sprannlib format */
char logn[100],posn[100],negn[100],undn[100],pngn[100];	/* for output filenames */
int tn1,tn2; 			/* number of sequences in the pos/neg training set */
int ga[DIM1+2][DIM2+2],nga[DIM1+2][DIM2+2];	/* GPIarray / NonGPIarray to map valid evalultion */
int mode;		/* mode of operation (search, train, web, interactive) */
NET *net;		/* the neural net pointer */
char posi[DIM1+2][DIM2+2];	/* array to hold the map */
long rG[DIM1+2][DIM2+2];	/* array to hold the G values */
long rN[DIM1+2][DIM2+2];	/* array to hold the N values */
long rU[DIM1+2][DIM2+2];	/* array to hold the U values */

int main(int argc, char *argv[]) {
	int i; /* loop counter */
    if (argc<2) {	/* less than one command line arg -> help/exit */
       	printf("\nSyntax:\tkohgpi t\n\tkohgpi <prot.fasta>\n\tkohgpi i\n\n");
		exit(-1);
	}
	if ((argv[1][1]==0) && (argv[1][0]=='t')) {mode=1;} /* train */
	else {mode=0;} /* search/evaluate */
	printf("Kohonen gpAI v1.5 - ");
	if (mode==1) {printf("training mode ");} else {printf("search mode ");}
	printf("- 3.3.2004 Nick Fankhauser\nDimensions: %d\n",INP);
	if ((argv[1][1]==0) && (argv[1][0]=='i')) {mode=2;interactive();} /* interactiv */
	if ((argv[1][1]==0) && (argv[1][0]=='d')) {mode=2;webdaemon();}	/* webdaemon */
	if (mode==1) {	/* training mode */
		if (load_fasta("trainp.txt",-1)==-1) {exit(-1);} /* load pos train in array */
		tn1=numseq;	/* store number of seq */
		if (load_fasta("trainn.txt",tn1-1)==-1) {exit(-1);} /* load neg train in array */
		tn2=numseq;	/* store number of seq */
		numseq=tn1+tn2; /* numseq is the size of train set, used in summon_dataset */
	}
	else { /* search mode */
		if (load_fasta(argv[1],-1)==-1) {exit(-1);};	/* load seqs in array */
		outfiles(argv[1]);				/* create output file names */
	}
	for (i=0;i<numseq;i++) {interface(i);} /* pass train pos/neg seqs through interface */
	printf("\n- Successfuly passed the Interface %d!\n",LOCI);
	if (mode==1) {train();} else {search();}
	return(0);
}

void outfiles(char filename[100]) {
	strcpy(logn,filename);strcat(logn,".log");	/* file for all map-pos&results */
	strcpy(posn,filename);strcat(posn,".pos");	/* file for positive seqs */
	strcpy(negn,filename);strcat(negn,".neg");	/* file for negative seqs */
	strcpy(undn,filename);strcat(undn,".und");	/* file for undecided seqs */
	strcpy(pngn,filename);strcat(pngn,".png");	/* file for result PNG */
}

void webdaemon(void) {
	FILE *stream;int i;
	mapnet(); /* loads map and neural net */
	printf("waiting for input.txt to appear...\n");
	while (0==0) {	/* loop forever (or until you kill it...) */
		sleep(1);	/* sleep for a second */
		if ((stream=fopen(INFILE,"r"))!=NULL) {	/* does infile exist? */
			close(stream);	/* if yes: close it. load_fasta opens it. */
			load_fasta(INFILE,-1);	/* load the infile in array */
			for (i=0;i<numseq;i++) {interface(i);} /* pass seqs through interface */
	    		printf("\n- Successfuly passed the Interface %d!\n",LOCI);
			outfiles(INFILE); /* create the output filenames */
			search();	/* evaluate these seqs */
			remove(INFILE);	/* del file, so it's only analysed once */
		}
	}
}

void interactive(void) {
	char fn[100];	/* string to hold the user input */
	int i; 		/* loop counter */
	mapnet();	/* loads map and neural net */
	while (0==0) {	/* loop forever, or until you enter q as filename */
		printf("Enter filename: ");
		scanf("%s",&fn); 	/* get filename from user */
		if ((fn[0]=='q') && (fn[1]==0)) {exit(0);} /* file is just a q? -> quit! */
		if (load_fasta(fn,-1)==0) {
			for (i=0;i<numseq;i++) {interface(i);}	/* pass seqs through interface */
		    	printf("\n- Successfuly passed the Interface %d!\n",LOCI);
			outfiles(fn);	/* create the output filenames */
			search();		/* evaluate these seqs */
    	}
    }
}

void train(void) {
	FILE *stream,*log,*nmap,*pmap;	/* file handle to save map/net and log */
	int *NetVect;		/* vector to describe neural network strucure */
	double *InpVect;	/* input vector */
	int i,j,k;		/* loop counters */
	int trn;		/* number of seqs in training set pos and neg */
	int n1,n2; 		/* number of sequences in the pos/neg validation set */
	int xp,yp;		/* x,y position of winning unit */
	int mxs,min,minp=0;	/* mix sum, minimal mix sum, min mix sum at cycle */
	struct unit *winner;	/* unit (neuron) pointer */
	long options,done=FALSE;	/* options flags for learning, done flag */
	long cycles=1,cyc=0;	/* nb of cycles per learn step, current cycle  */
	char sv[10];		/* string to store whether it saved this cycle */
	summon_dataset();	/* converts seqs from input vector to sprannlib dataset */
	printf("- Converted into dataset.\n");trn=numseq;
	if (load_fasta("validp.txt",trn-1)==-1) {exit(-1);} /* load valid pos in array */
	n1=numseq;	/* save number of seqs */
	if (load_fasta("validn.txt",trn+n1-1)==-1) {exit(-1);}  /* load valid neg in array */
	n2=numseq;	/* save number of seqs */
	for (i=trn;i<(trn+n1+n2);i++) {interface(i);} /* pass valid pos/neg though interface */
	printf("\n- Successfuly passed the Interface %d!\n",LOCI);
	printf("Switching to daemon mode - training neural network...\n");
	daemon(1,0);
	sprinit(INIT_EXIT|INIT_SILENT);		/* sprannlib init */
	min=trn+n1+n2;	/* set min mix sum to maximum possible */
  	NetVect=ivector(1,2);	/* create vector for two integers */
  	NetVect[1]=DIM1;NetVect[2]=DIM2;	/* put net dimensions into vector */
	net=create_koh_net(1L,INP,2,NetVect);	/* creates 2D kohonen net with ID=1 */
	check_all_network(net);		/* check the network for errors */
	unif_rand_net(1,-.001,0.001,net);	/* randomize network weights */
	InpVect=vector(1,INP);		/* create input vector */
	options=(NO_HISTU|NO_HISTT|NO_HISTW|KOH_KERNEL_UPD);	/* no history and no limit */
	log=fopen("kohgpi.log","w");fclose(log);		/* create logfile */
	while (cyc<MAXCYC) {					/* loop MAXCYC times */
		if (koh_learn(net,set,LEARN,5,cycles,options)) { /* train net with dataset */
			sprerror("Error2");exit(-1);	/* quit on error. very unlikely */
		}
		cyc+=cycles;mxs=0;		/* increase cycle counter and reset mix sum */
		for (i=0;i<DIM1+1;i++) {	/* resets pos/neg winner array */
			for (j=0;j<DIM2+1;j++) {ga[i][j]=0;nga[i][j]=0;}
		}
   		for (k=0;k<trn;k++) {		/* for all seqs in the training set */
			for (i=1;i<INP+1;i++) {InpVect[i]=x[k][i-1];}	/* copy input in vector */
			eval_koh_net(net,InpVect);		/* evaluate net with this input */
			winner=best_matching_unit(net);	/* get the winning unit for this input */
			parse(winner->UnitName,&xp,&yp);	/* get winning units coordinates */
			if (k<tn1) {ga[xp][yp]++;}		/* put pos of positive seqs into pos map */
			else {nga[xp][yp]++;}		/* put pos of negative seqs into neg map */
		}
   		for (k=trn;k<trn+n1+n2;k++) {		/* for all seqs in the validation set */
			for (i=1;i<INP+1;i++) {InpVect[i]=x[k][i-1];}	/* copy input in vector */
			eval_koh_net(net,InpVect);		/* evaluate net with this input */
			winner=best_matching_unit(net);	/* get the winning unit for this input */
			parse(winner->UnitName,&xp,&yp);	/* get winning units coordinates */
			if (k<trn+n1) {ga[xp][yp]++;}		/* put pos of positive seqs into pos map */
			else {nga[xp][yp]++;}		/* put pos of negative seqs into neg map */
		}
		for (i=1;i<DIM1+1;i++) {			/* loop first dimention */
			for (j=1;j<DIM2+1;j++) {		/* loop second dimention */
				if ((ga[i][j]>0) && (nga[i][j]>0)) { /* unit won on pos AND neg seqs */
					mxs+=ga[i][j]+nga[i][j];	/* then sum up the nb of seqs than won */
				}
			}
		}
		if (mxs<min) {		/* do we have a new minimum? */
			stream=fopen("kohgpi.net","w");	/* open NN save file */
			fprintf_network(stream,net,(NO_HISTU | NO_HISTT | NO_HISTW));	/* write NN */
			fclose(stream);	/* close stream */
			strcpy(sv,"Yes");	/* Yes, we saved */
			min=mxs;minp=cyc;	/* set new minimum and minimum cycle position */
			pmap=fopen("kohgpi.pmp","w");	/* open pos map save file */
			nmap=fopen("kohgpi.nmp","w");	/* open neg map save file */
			for (i=1;i<=DIM1;i++) {		/* loop first dimention */
				for (j=1;j<=DIM2;j++) {	/* loop second dimention */
					fprintf(pmap,"%d ",ga[i][j]);	/* write char to stream */
					fprintf(nmap,"%d ",nga[i][j]);	/* write char to stream */
				}
				fprintf(nmap,"\n");fprintf(pmap,"\n");
			}
			fclose(pmap);fclose(nmap);	/* close stream */
			plot(cyc);			/* save network as PNG image */
		} else {strcpy(sv,"No");}		/* if mo new minimum, we don't save */
		log=fopen("kohgpi.log","a");	/* open the log to append and write status */
		fprintf(log,"cyc:%d mix:%d min:%d n:%d minp:%d save:%s\n",cyc,mxs,min,trn+n1+n2,minp,sv);
		printf("cyc:%d mix:%d min:%d n:%d minp:%d save:%s\n",cyc,mxs,min,trn+n1+n2,minp,sv);
		fclose(log);	/* close the logfile */
	}
}

void parse(char uname[80],int *xp,int *yp) {
	int i,flg=0,xcp=0,ycp=0; /* loop counter, flag, x/y string pointer */
	char xc[10],yc[10];		/* x/y string */
	for (i=0;i<strlen(uname);i++) {	/* go through the unit name */
		if (uname[i]==')') {flg=0;}	/* after the final paranthese, no more digits follow */
		if (flg==2) {		/* if we're the y number */
			yc[ycp]=uname[i];	/* copy y digits into y string */
			ycp++;		/* increase y string pointer */
			yc[ycp]='\0';	/* null terminate the y string */
		}	
		if (uname[i]==',') {flg=2;}	/* after the comma, there is the y number */
		if (flg==1) {		/* if we're the x number */
			xc[xcp]=uname[i];	/* copy x digits into x string */
			xcp++;		/* increase x string pointer */
			xc[xcp]='\0';	/* null terminate the x string */
		}
		if (uname[i]=='(') {flg=1;}	/* after the first paranthese, there is the x number */
	}
	*xp=atoi(xc);*yp=atoi(yc);	/* return x and y */
}

void plot(long c) {
	gdImagePtr im;	/* GD library image pointer */
	FILE *pngout;	/* file handle for png out */
	char fn[30],zahl[10];	/* string var for filename and int->str conversions */
	int i,j,k;		/* loop counters */
	int black,white,col,farb;	/* color variables */
	mkdir("png",0777);
	strcpy(fn,"png/kohgpi");	/* filename starts with kohgpi */
	sprintf(zahl,"%d",DIM1);strcat(fn,zahl);	/* add DIM1 to filename */
	strcat(fn,"x");				/* seperate dimensions with an x */
	sprintf(zahl,"%d",DIM2);strcat(fn,zahl);	/* add DIM2 to filename */
	strcat(fn,"_");				/* seperate with an underscore */
	sprintf(zahl,"%d",(int)(1/LEARN));strcat(fn,zahl); /* add reciprocal of learn rate */ 
	strcat(fn,"_");				/* seperate with an underscore */
	sprintf(zahl,"%d",c);strcat(fn,zahl);		 /* add nb of current cycle */
	strcat(fn,".png");				/* add extension */
	im=gdImageCreate(DIM1*PIX,DIM2*PIX);		/* create image with dim*pix size */
	black=gdImageColorAllocate(im,255,255,255);		/* allocate white as background */
	for (i=1;i<=DIM1;i++) {					
		for (j=1;j<=DIM2;j++) {		/* loop second dimention */
			if ((ga[i][j]==0) && (nga[i][j]==0)) {col=black;}	/* nothing -> black */
			if ((ga[i][j]>0) && (nga[i][j]==0)) {	/* GPI hits */
				farb=ga[i][j]*15+100; /* set intensity according to nb of seqs on unit */
				if (farb>255) {farb=255;}	/* but not bigger than 255 */
				col=gdImageColorAllocate(im,0,farb,0);	/* set this as green value */
			}
			if ((ga[i][j]==0) && (nga[i][j]>0)) {	/* non GPI hits */
				farb=nga[i][j]*15+100; /* set intensity according to nb of seqs on unit */
				if (farb>255) {farb=255;}	/* but not bigger than 255 */
				col=gdImageColorAllocate(im,farb,0,0);	/* set this as red value */
			}
			if ((ga[i][j]>0) && (nga[i][j]>0)) {	/* mixed hits */
				farb=(ga[i][j]+nga[i][j])*15+100;	/* set intensity prop to nb of seqs */
				if (farb>255) {farb=255;}		/* but not bigger than 255 */
				col=gdImageColorAllocate(im,farb,farb,0);	/* set this as yellow value */
			}
			for (k=0;k<PIX;k++) {	/* loop PIX times */
			/* creates the point for a neuron with color selected earlier */
			gdImageLine(im,(i-1)*PIX+k,(j-1)*PIX,(i-1)*PIX+k,(j-1)*PIX+PIX,col);
			}
			/* lines to separate the point, so it looks better */
			gdImageLine(im,(i-1)*PIX,(j-1)*PIX,(i-1)*PIX,(j-1)*PIX+PIX,black);
			gdImageLine(im,(i-1)*PIX,(j-1)*PIX,(i-1)*PIX+PIX,(j-1)*PIX,black);
		}
	}
	pngout=fopen(fn,"wb");gdImagePng(im,pngout);	/* writes image to file */
	fclose(pngout);gdImageDestroy(im);		/* close file and destroy image in memory */
}

void mapnet(void) {
	FILE *stream;	/* file handle */
	int i,k;		/* loop counters */
	printf("- loading map\n");
	if ((stream=fopen("kohgpi.map","r"))==NULL) {	/* open saved map */
		printf("map not found\n");exit(-1);	/* quit if it's not here */
	}
	i=0;
	while (!feof(stream)) {			/* read file till end */
		fgets(posi[i],sizeof(posi[i]),stream);	/* get a line of the map */
		posi[i][strlen(posi[i])-1]='\0';	/* null terminate it */
		for (k=0;k<strlen(posi[i]);k++) {	/* for all chars of the line */
			printf("%c ",posi[i][k]);	/* print the char */
		}
		printf("\n");i++;			/* print new line */
	}
	fclose(stream);				/* close the file */
	sprinit(INIT_EXIT|INIT_SILENT);		/* init sprannlib */
	printf("- loading network...\n");
	if ((stream=fopen("kohgpi.net","r"))==NULL) {	/* open saved neural net */
		printf("net not found\n");exit(-1);	/* exit if it's not here */
	}
	net=fscanf_network(stream);fclose(stream);	/* load the net and close file */
	printf("- Network loaded\n");
	check_all_network(net);			/* check network for errors  */
}

void search(void) {
	FILE *log,*pos,*neg,*und;	/*  file handles for output files */
	int i,j,k,rad;		/* loop counters */
	int xp,yp,gc,nc;	/* x/y position, pos/neg counters in map area */
	double *InpVect;	/* input vector for evaluation */
	float error;
	struct unit *winner;	/* unit (neuron) pointer */
	long gpi=0;		/* pos counter */
	char result;
	char qual[100];
	if (mode==0) {mapnet();}	/* load map/net when in search mode */
   	InpVect=vector(1,INP);		/* create input vector */
	log=fopen(logn,"w");pos=fopen(posn,"w"); /* open log and positive files for writing */
	neg=fopen(negn,"w");und=fopen(undn,"w"); /* open negative files for writing */
	for (i=1;i<=DIM1;i++) {
		for (j=1;j<=DIM2;j++) {		/* loop second dimention */
			rG[i][j]=0;rN[i][j]=0;rU[i][j]=0;
		}
	}
   	for (k=0;k<numseq;k++) {			/* loop through all seqs to be avaluated */
		for (i=1;i<INP+1;i++) {		/* for all input values */
			InpVect[i]=x[k][i-1];	/* copy value into input vector */
		}
		eval_koh_net(net,InpVect);		/* evaluate the network */
		winner=best_matching_unit(net);	/* get a pointer to the best matching unit */
		parse(winner->UnitName,&xp,&yp);	/* get its coordinates from its name */
		sprintf(qual," [Cleavage at C-%d of %d%% quality]",ome[k],ope[k]);
		if (ope[k]==0) {strcpy(qual,"");}
		fprintf(log,"%s%s\n",desc[k],qual);	/* print sequence description to log */
		switch (posi[xp-1][yp-1]) {
			case 'G':
				gpi++;result='G';
				/* write seq to positive file */
				fprintf(pos,"%s%s\n%s\n",desc[k],qual,seq[k]);
				break;
			case 'N':
				result='N';
				/* write seq to negative file */
				fprintf(neg,"%s\n%s\n",desc[k],seq[k]);
				break;
			case 'U':
				result='U';
				/* write seq to undecided file */
				fprintf(und,"%s\n%s\n",desc[k],seq[k]);
				break;
			default:
				/* made obsolete by mapfill: better weighting */
				fprintf(log,"Unmapped Sequence.\n");
				gc=0;nc=0;rad=2;	/* reset counters */
				/* loop through the radius around the hit */
				for (i=xp-1-rad;i<xp+rad;i++) {
					for (j=yp-1-rad;j<yp+rad;j++) {
						if ((i>-1) && (j>-1) && (i<DIM1) && (j<DIM2)) {
							fprintf(log,"%c ",posi[i][j]);
							/* count the Gs around the hit */
							if (posi[i][j]=='G') {gc++;}
							/* count the Ns around the hit */
							if (posi[i][j]=='N') {nc++;}
						}
					}
					fprintf(log,"\n");
				}
				fprintf(log,"rad: %d gc: %d nc: %d\n",rad,gc,nc);
				if (gc>nc) {
					/* increase the GPI count */
					gpi++;result='G';
					/* write seq to positive file */
					fprintf(pos,"%s [unmap%d:%d]\n%s\n",desc[k],gc,nc,seq[k]);
				}
				if (gc<=nc) {
				 	/* write seq to negative file */
					fprintf(neg,"%s [unmap%d:%d]\n%s\n",desc[k],gc,nc,seq[k]);
					result='N';
				}
		}
		/* write position and result to log */
		if (result=='G') {rG[xp][yp]++;}
		if (result=='N') {rN[xp][yp]++;}
		if (result=='U') {rU[xp][yp]++;}
		fprintf(log,"seq%d x:%d y:%d map:%c omega:%d oq:%d\n\n",k,xp,yp,result,ome[k],ope[k]);
	}
	fclose(log);fclose(pos);fclose(neg);fclose(und); /* close all output files */
	resplot();
	if (gpi>numseq-gpi) {error=(float)(numseq-gpi)/numseq*100;}
	else {error=(float)gpi/numseq*100;}
	printf("GPI: %d None: %d -> %2.2f %%\n",gpi,numseq-gpi,error);
}

void resplot(void) {
	gdImagePtr im;	/* GD library image pointer */
	FILE *pngout;	/* file handle for png out */
	int i,j,k;		/* loop counters */
	int hit;
	int black,white,col,farb;	/* color variables */
	im=gdImageCreate(DIM1*PIX,DIM2*PIX);		/* create image with dim*pix size */
	black=gdImageColorAllocate(im,0,0,0);		/* allocate white as background */
	for (i=1;i<=DIM1;i++) {					
		for (j=1;j<=DIM2;j++) {		/* loop second dimention */
			col=black;hit=0;	/* nothing -> black */
			if (posi[i-1][j-1]=='G') {col=gdImageColorResolve(im,0,250,0);}
			if (posi[i-1][j-1]=='N') {col=gdImageColorResolve(im,0,0,250);}
			if (posi[i-1][j-1]=='U') {col=gdImageColorResolve(im,250,0,0);}
			if ((rG[i][j]>0)) {	/* GPI hits */
				farb=250-rG[i][j]*10; /* set intensity according to nb of seqs on unit */
				if (farb<MINC) {farb=MINC;}	/* but not smaller than 20 */
				col=gdImageColorResolve(im,0,farb,0);	/* set this as green value */
				hit=1;
			}
			if ((rN[i][j]>0)) {	/* non GPI hits */
				farb=250-rN[i][j]*10; /* set intensity according to nb of seqs on unit */
				if (farb<MINC) {farb=MINC;}	/* but not smaller than 20 */
				col=gdImageColorResolve(im,0,0,farb);	/* set this as red value */
				hit=1;
			}
			if ((rU[i][j]>0)) {	/* mixed hits */
				farb=250-rU[i][j]*10;	/* set intensity prop to nb of seqs */
				if (farb<MINC) {farb=MINC;}	/* but not smaller than 20 */
				col=gdImageColorResolve(im,farb,0,0);	/* set this as yellow value */
				hit=1;
			}
			gdImageFilledRectangle(im,(i-1)*PIX,(j-1)*PIX,
			(i-1)*PIX+PIX,(j-1)*PIX+PIX,col);
			if (hit==1) {
				gdImageLine(im,(i-1)*PIX+ZD,(j-1)*PIX+ZD,
				(i-1)*PIX+PIX-ZD,(j-1)*PIX+PIX-ZD,black);
				gdImageLine(im,(i-1)*PIX+PIX-ZD,(j-1)*PIX+ZD,
				(i-1)*PIX+ZD,(j-1)*PIX+PIX-ZD,black);
				gdImageLine(im,(i-1)*PIX+ZD,(j-1)*PIX+PIZ,
				(i-1)*PIX+PIX-ZD,(j-1)*PIX+PIZ,black);
				gdImageLine(im,(i-1)*PIX+PIZ,(j-1)*PIX+ZD,
				(i-1)*PIX+PIZ,(j-1)*PIX+PIX-ZD,black);
			}
			/* lines to separate the point, so it looks better */
			gdImageLine(im,(i-1)*PIX,(j-1)*PIX,(i-1)*PIX,(j-1)*PIX+PIX,black);
			gdImageLine(im,(i-1)*PIX,(j-1)*PIX,(i-1)*PIX+PIX,(j-1)*PIX,black);
		}
	}
	pngout=fopen(pngn,"wb");gdImagePng(im,pngout);	/* writes image to file */
	fclose(pngout);gdImageDestroy(im);		/* close file and destroy image in memory */
}

void summon_dataset(void) {
	SAMPLE  *sample_ptr = NULL;	/* pointer to allocate samples */
	int i,k;			/* loop counters */
	FILE *ofp;		/* output stream for dataset. output is just needed for debuging */
	float f;			/* var to temporaly hold the output for a sample */
	char str[50];		/* used to temporaly hold strings. direct copy gave a segfault... */
	set = malloc_dataset();	/* creates the new dataset */
	strcpy(str,"GPI Training Sequence Data");	/* direct copy leads to segfault */
	strcpy(set->SetText,str);	/* describes the set */
	strcpy(str,"today");
	strcpy(set->CreationDate,str);	/* dataset creation date */
	set->SetId=1L;			/* dataset identification number. dont know what for... */
	set->SetFlag=(IOSET|LEARNSET|NOSTATIST);/* set type of dataset */
	set->NumOutputs=2L;	/* number of output values. not needed for kohonen net. */
	set->OutputDescr=(char**)malloc(sizeof(char*)*2);	/* allocate pointer to 2 pointers */
	set->OutputDescr[0]=malloc(sizeof(char)*8);	/* allocate string on the first pointer */
	set->OutputDescr[1]=malloc(sizeof(char)*8);	/* allocate string on the second pointer */
	strcpy(set->OutputDescr[0],"GPI");		/* first output is GPI */
	strcpy(set->OutputDescr[1],"NO GPI");		/* first output is not GPI */
	set->NumSamples=(long)numseq;			/* set number of samples */
	set->NumInputs=(long)INP;			/* set number of inputs */
	set->FirstSample=NULL;			/* doesn't make sense, but it's in the example... */
	for(i=0;i<numseq;i++) {			/* for all sequences */
		/* if it's the first sample, allocate & save pointer in FirstSample and sample_ptr*/
	   	if (sample_ptr==NULL) {set->FirstSample=sample_ptr=malloc_sample();}
		/* if not, allocate next sample and save the pointer in Next and then sample_ptr */
		else {sample_ptr->Next=malloc_sample();sample_ptr=sample_ptr->Next;}
      		sample_ptr->SampleId=i;		/* set sample number */
		sample_ptr->SampleFlag=SAMPLE_ENABLED;	/* enable sample */
      		sample_ptr->Input=vector(1,INP);	/* define Input as vector */
		for(k=1;k<INP+1;k++) {					/* for all samples */
			sample_ptr->Input[k]=x[i][k-1];		/* copy samples into Input vector */
		}
		sample_ptr->Output.Vector=vector(1,2);	/* define Output.vector as vector */
     		if (i<tn1) {f=0.9;} else {f=0.1;}	/* GPI -> 0.9, not GPI -> 0.1 */
		sample_ptr->Output.Vector[1]=f;	/* set the GPI output node */
		sample_ptr->Output.Vector[2]=1-f;	/* set the not GPI output node */
   	}
	random_order_dataset(set);	/* shuffle the dataset */
	ofp=fopen("kohgpi.dat","w");fprintf_dataset(ofp,set);fclose(ofp);	/* save dataset */
}
