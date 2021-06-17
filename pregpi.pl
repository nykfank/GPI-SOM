$f=shift;open(F,$f) or die;open(O,">$f.cln") or die;
$i=-1;$minlen=32;$maxlen=6000;
while (<F>) {chomp;if (/^>/) {$i++;$d[$i]=$_;} else {$s[$i]=$s[$i].$_}}
$numseq=$i+1;
for ($i=0;$i<$numseq;$i++) {
	$s[$i]=~tr/a-z/A-Z/;$s[$i]=~s/[^A-Z]//g;
	if ((length($s[$i])>$minlen) and (length($s[$i])<$maxlen)) {print O "$d[$i]\n$s[$i]\n";$n++}
}
close F;close O;print "$n of $numseq seqs written\n";
system("mv $f $f.old");system("mv $f.cln $f");
