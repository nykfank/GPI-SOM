#!/usr/bin/perl

# GPI-webserver control daemon - 6.10.2003 Nick Fankhauser

$bd="/var/www/cgi-bin";
&setlog;
$seqlimit=39000;        # kohgpi max is 40000
$w="waiting for unprocessed files...\n";print $w;
while (0==0) {
        sleep(2);
        $i=1;
        while ((-r "$bd/data/job$i.fasta.pos.gpi") and (-r "$bd/data/job$i.fasta")) {$i++}
        if (-r "$bd/data/job$i.fasta.go") {process("$bd/data/job$i.fasta")}
}

sub setlog {$logfile="$bd/logs/gpid.log"}

sub process {
        $fn=shift;$number=$i;
        if (exe("sreformat fasta $fn >$fn.srf")) {
                logg("$fn problem with input\n");
                open(O,">$fn.pos.gpi") or die;close O;
                goto fin;
        }
        $total=count("$fn.srf");
	if ($total<0) {logg("$fn.srf not found\n");&jobterm;goto fin}
        logg("$fn total: $total seqs\n");
        if (pregpi("$fn.srf")) {
                logg("$fn problem with input\n");
                open(O,">$fn.pos.gpi") or die;close O;
                goto fin;
        }
        $total=count("$fn.srf");
	if ($total<0) {logg("$fn.srf not found\n");&jobterm;goto fin}
        logg("$fn pregpi: $total seqs\n");
        if ($total>$seqlimit) {
                logg("$fn too many seqs\n");
                open(O,">$fn.pos.gpi") or die;close O;
                goto fin;
        }
        exe("cp $fn.srf input.txt");
        $pos=0;$neg=0;$und=0;$lp=0;
        until ($total==$pos+$neg+$und) {
                sleep(2);
                $pos=count("$bd/input.txt.pos");
		if ($pos<0) {logg("$bd/input.txt.pos not found\n");&jobterm;goto fin}
		$neg=count("$bd/input.txt.neg");
		if ($neg<0) {logg("$bd/input.txt.neg not found\n");&jobterm;goto fin}
		$und=count("$bd/input.txt.und");
		if ($und<0) {logg("$bd/input.txt.und not found\n");&jobterm;goto fin}
                $lp++;
                if ($lp>5) {    # quits, when kohgpi doesnt complete all seqs
                        $lp=0;
                        if ($lsum==$pos+$neg+$und) {
                                logg("$fn problem with kohgpi\n");
                                open(O,">$fn.pos.gpi") or die;
                                close O;goto fin;
                        }
                        $lsum=$pos+$neg+$und;
                }

        }
        logg("$fn pos: $pos seqs\n");
	logg("$fn neg: $neg seqs\n");
	logg("$fn und: $und seqs\n");
        exe("mv input.txt.pos $fn.pos");
        exe("mv input.txt.neg $fn.neg");
	exe("mv input.txt.und $fn.und");
        exe("mv input.txt.log $fn.log");
	exe("mv input.txt.png $fn.png");
	exe("mv $fn.png ../html/gif");
	&imgmap;
        if (sig("$fn.pos")) {
		&setlog;
                logg("$fn problem with signalp\n");
                &jobterm;goto fin;
        }
        $gpi=count("$fn.pos.nt.log");
	if ($gpi<0) {logg("$fn.pos.nt.log not found\n");&jobterm;goto fin}
        logg("$fn gpi: $gpi seqs\n");
        lookup("$fn.pos.nt.log","$fn.pos");
        exe("sreformat fasta $fn.pos.gpi >$fn.pos.gpi.srf");
        exe("mv $fn.pos.gpi.srf $fn.pos.gpi");
        fin:
        print $w;
}

sub imgmap {
	my $px=10;
	# build hash-reference of used units
	%hr=();
	open(F,"$fn.log") or die;
	while (<F>) {
		if (/x:(\d+) y:(\d+)/) {$hr{"$1-$2"}=1}
	}
	close F;
	# create image map of used units
	open(O,">$fn.imp") or die;
	for ($i=0;$i<40;$i++) {
                for ($j=0;$j<40;$j++) {
			$ii=$i+1;$jj=$j+1;
			if ($hr{"$ii-$jj"}) {
				print O "<area shape=rect coords=";
        	                print O $i*$px,",",$j*$px,",",$i*$px+$px,",",$j*$px+$px;
                	        $ii=$i;$jj=$j;
				unless ($ii) {$ii="O"}
				unless ($jj) {$jj="O"}
				print O " href=?id=$number&x=$ii&y=$jj>\n";
			}
                }
        }
	close O;
}

sub jobterm {open(O,">$fn.pos.gpi");close O}

sub pregpi {
	my $f=shift;open(F,$f) or return(1);open(O,">$f.cln") or return(1);
	my $i=-1;my $minlen=32;my $maxlen=6000;my @s=();my @d=();my $n=0;
	while (<F>) {chomp;if (/^>/) {$i++;$d[$i]=$_;} else {$s[$i]=$s[$i].$_}}
	my $numseq=$i+1;
	for ($i=0;$i<$numseq;$i++) {
		$s[$i]=~tr/a-z/A-Z/;$s[$i]=~s/[^A-Z]//g;
		if ((length($s[$i])>$minlen) and (length($s[$i])<$maxlen)) {
			print O "$d[$i]\n$s[$i]\n";$n++;
		}
	}
	close F;close O;logg("$n of $numseq seqs written\n");
	system("mv $f $f.old");system("mv $f.cln $f");
	return(0);
}

sub sig {
	my $f=shift;$fna="$f.nt";
	open(F,$f) or return(1);open(O,">$fna") or return(1);
	my $i=-1;my $minlen=30;my $cterm=70;my @d=();my @s=();
	while (<F>) {chomp;if (/^>/) {$i++;$d[$i]=$_;} else {$s[$i]=$s[$i].$_}}
	my $numseq=$i+1;my $n=0;
	for ($i=0;$i<$numseq;$i++) {
		$s[$i]=~s/\*//g;
		if (length($s[$i])>$minlen) {
			my @sar=split("",$s[$i]);my $beg=0;my $seq="";
			for ($j=$beg;$j<$beg+$cterm;$j++) {$seq=$seq.$sar[$j]}
			print O "$d[$i]\n$seq\n";$n++;
		}
	}
	close F;close O;logg("$n of $numseq seqs written\n");
	open(F,$fna) or return(1);$flg=1;$nb=1;
	
	while (<F>) {
		if (/^>/) {$cnt++}
		if ($cnt>599) {&proc}
		if ($flg) {open(O,">$fna.$nb") or return(1);$flg=0;}
		print O $_;
	}
	if ($cnt) {&proc}
	exe("cat $fna.*.sig >$fna.sig");
	logg("total GPI in $fna: $sum\n",1);
	&setlog;
	return(0);
}

sub proc {  # used by sig; needs global $fna,$nb, $cnt, $flg, $sum
	close O;
	exe("/usr/local/signalp-2.0/signalp -t euk $fna.$nb >$fna.$nb.sig");
	open(S,"$fna.$nb.sig") or return;my $c=0;
	$logfile="$fna.log";
	while (<S>) {
	        if (/^>/) {$d5=$d4;$d4=$d3;$d3=$d2;$d2=$d1;$d1=$_}
	        if (/Prediction: Signal/) {logg($d5,1);$c++}
	}
	close S;&setlog;
	logg("Proteins with signal sequence in $fna.$nb: $c\n");
	$cnt=0;$nb++;$flg=1;$sum+=$c;
}

sub lookup {
	my $df=shift;my $sf=shift;my $of=$sf.".gpi";
	open(F,$sf) or return(1);my @s=();my @d=();
	my $i=-1;while (<F>) {chomp;if (/^>/) {$i++;$d[$i]=$_;} else {$s[$i]=$s[$i].$_}}
	my $ns=$i+1;close F;
	open(F,$df) or return(1);open(O,">$of") or return(1);
	logg("$ns seqs read from $sf. searching...\n");
	my $lc=0;my $put=0;
	while (<F>) {
		chomp;
		if (/^>/) {
			$lc++;
			for ($i=0;$i<$ns;$i++) {
				if (/\Q$d[$i]\E/) {print O "$d[$i]\n$s[$i]\n";$put++}
			}
		}
	}
	logg("$put of $lc seqs written to $of\n");
}

sub logg {
        if (length($_[0])>1) {
                print $_[0];my $t=zeit();
		if ($_[1]) {$t=""}
                unless (-r $logfile) {open(L,">$logfile");close L}
                open(L,">>$logfile");print L $t,$_[0];close L;
        }
}
   
sub zeit {      # returns current time like 18:37:12 in a string
        ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
         return "$mday.",$mon+1,".",$year+1900," $hour:$min ";
}
   
sub exe {my $ret=1;logg("Run: $_[0] ");$ret=system($_[0]);logg("-> $ret\n");return $ret}
  
sub count {
        open(C,shift) or return(-1);my $cnt=0;
        while (<C>) {chomp;if (/^>/) {$cnt++}}
        close C;return $cnt;
}
