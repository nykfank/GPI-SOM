$f=shift;$fn=$f.".nt";
open(F,$f) or die;open(O,">$fn") or die;
$i=-1;$minlen=30;$cterm=70;
while (<F>) {chomp;if (/^>/) {$i++;$d[$i]=$_;} else {$s[$i]=$s[$i].$_}}
$numseq=$i+1;
for ($i=0;$i<$numseq;$i++) {
	$s[$i]=~s/\*//g;
	if (length($s[$i])>$minlen) {
		@sar=split("",$s[$i]);$beg=0;$seq="";
		for ($j=$beg;$j<$beg+$cterm;$j++) {$seq=$seq.$sar[$j]}
		print O "$d[$i]\n$seq\n";$n++;
	}
}
close F;close O;print "$n of $numseq seqs written\n";
open(F,$fn) or die;$flg=1;$nb=1;$logfile="$fn.log";
while (<F>) {
	if (/^>/) {$cnt++}
	if ($cnt>599) {&proc}
	if ($flg) {open(O,">$fn.$nb") or die;$flg=0;}
	print O $_;
}
if ($cnt) {&proc}
exe("cat $fn.*.sig >$fn.sig");
logg("total GPI in $fn: $sum\n");

sub proc {
	close O;
	exe("/usr/local/signalp-2.0/signalp -t euk $fn.$nb >$fn.$nb.sig");
	open(S,"$fn.$nb.sig") or die;$c=0;
	while (<S>) {
	        if (/^>/) {$d5=$d4;$d4=$d3;$d3=$d2;$d2=$d1;$d1=$_}
	        if (/Prediction: Signal/) {logg($d5);$c++}
	}
	close S;logg("Proteins with signal sequence in $fn.$nb: $c\n");
	$cnt=0;$nb++;$flg=1;$sum+=$c;
}
sub exe {
	my $ret=1;
	logg("Run: $_[0] ");
	$ret=system($_[0]);
	logg("-> $ret\n");
	return $ret;
}

sub logg {
	if (length($_[0])>1) {
		print $_[0];
		unless (-r $logfile) {open(L,">$logfile");close L;}
		open(L,">>$logfile");print L $_[0];close L;
	}
}
