$df=shift;$sf=shift;$of=$sf.".gpi";
unless ($sf) {print "lookup.pl <desc.file> <seq.file>\n\n";exit}
open(F,$sf) or die;
$i=-1;while (<F>) {chomp;if (/^>/) {$i++;$d[$i]=$_;} else {$s[$i]=$s[$i].$_}}
$ns=$i+1;close F;open(F,$df) or die;open(O,">$of") or die;
print "$ns seqs read from $sf. searching...\n";
while (<F>) {
	chomp;
	if (/^>/) {
		$lc++;
		for ($i=0;$i<$ns;$i++) {
			if (/\Q$d[$i]\E/) {print O "$d[$i]\n$s[$i]\n";$put++}
		}
	}
}
print "$put of $lc seqs written to $of\n";
