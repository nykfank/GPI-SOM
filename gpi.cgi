#!/usr/bin/perl
# GPI-webserver user interface - 6.10.2003 Nick Fankhauser
 
use CGI qw(:standard);
$ofn="input.fasta";
$linelimit=10000;
$mail="niklaus.fankhauser\@izb.unibe.ch";
$logfile="logs/gpid.log";
$url="http://sublx.unibe.ch";
$title="KohGPI: Identification of GPI-anchor signals<br>by a Kohonen Self Organizing Map";
$ch=param("ch");$id=param("id");$seq=param("seq");$fasta=param("fasta");$ref=param("ref");
$arch=param("arch");$jobnam=param("jobnam");
print header(-type => 'text/html');print "\n<!-- code by nick fankhauser -->\n\n";
$t=$title;$t=~s/\<br\>/ /g;
print "<HTML><HEAD><TITLE>$t</TITLE>\n";
if (-r "data/job$id.fasta.pos.gpi") {$ref=0}
if (($seq) or ($fasta)) {&kohonen}
if ($ref) {print "<meta http-equiv=\"refresh\" content=\"2; URL=gpi.cgi?ref=$ref&id=$id\">\n"}
print "</head><body>\n";
if (($seq) or ($fasta)) {
        print "Job #$id started! -> <a href=$url/cgi-bin/gpi.cgi?id=$i&ref=1>Results</a><br>\n";
        goto fuss;
}
if ($arch) {&archiv;goto fuss}
if ($ch) {&show;goto fuss}
if ($id) {&result;goto fuss}
&iform;&refs;
fuss:
print "</body></html>\n";
 
sub archiv {
        print "<a href=http://gpi.unibe.ch><img src=$url/gpis.gif width=100 border=0>
	</a><p><h2>Archive of analyzed sequences</h2>
        <table border=1><tr><td><b>Job</b></td><td><b>Submitted by</b></td><td><b>Date</b></td></tr>\n";
        open(L,"logs/gpi.log") or die;
        while (<L>) {
                $nam="";
                if (/Job(\d+?):/) {$id=$1}
                if (/^(.+?) Job/) {$dat=$1}
                if (/Job$id: (.+?)\(/) {$host=$1}
                if (/Jobnam: (.+?)$/) {$nam=$1}
                if ($nam) {print "<tr><td><a href=$url/cgi-bin/gpi.cgi?id=$id>$nam</a></td><td>$host</td><td>$dat</td></tr>\n"}
        }
        close L;print "</table>\n";
}
 
sub show {
        print "<pre>\n";$c=0;
        open(F,"data/job$id.fasta$ch") or die;
        while (<F>) {
                if ($c<$linelimit) {print $_}
                $c++;
                if ($c>$linelimit) {
                        print "</pre><h2>Only $linelimit lines shown. Need more? <a href=mailto:$mail>$mail</a></h2>";
                        goto fuss;
                }
        }
        close F;print "</pre>";
}
 
sub result {if (-r "data/job$id.fasta.pos.gpi") {&finished} else {&waiting}}
 
sub finished {
	print "<table><tr><td><a href=http://gpi.unibe.ch>
        <img src=$url/gpis.gif width=100 border=0></a><p>";
        open(L,"logs/gpi.log") or die;
        while (<L>) {chomp;if (/Job$id:/) {if (/Jobnam: (.+?)$/) {$nam=$1}}}
        close L;if ($nam) {$n=" ($nam)"} else {$n=""}
        print "<h2>Results for Job #$id$n:</h2><p>\n";
        open(L,$logfile) or die;
        while (<L>) {
                if (/job$id\.fasta total: (\d+) /) {$total=$1}
                if (/job$id\.fasta pregpi: (\d+) /) {$pregpi=$1}
                if (/job$id\.fasta pos: (\d+) /) {$pos=$1}
                if (/job$id\.fasta neg: (\d+) /) {$neg=$1}
                if (/job$id\.fasta gpi: (\d+) /) {$gpi=$1}
		if (/job$id\.fasta und: (\d+) /) {$und=$1}
                if (/job$id\.fasta problem with input/) {$prob=1}
                if (/job$id\.fasta problem with kohgpi/) {$prob=2}
                if (/job$id\.fasta too many seqs/) {$prob=3}
        }
        close L;$anf="<a href=$url/cgi-bin/gpi.cgi?id=$id&ch=";
        if ($prob==1) {print "<h2>Invalid input sequence!</h2>\n";goto fuss}
        if ($prob==2) {
                print "<h2>kohGPI wasnt able to process all input sequences.</h2>
                Something unexpected is contained within your fasta file!<br>\n
                Please report this problem along with the file causing it to 
                <a href=mailto:$mail>$mail</a>.\n";
                goto fuss;
        }
        if ($prob==3) {print "<h2>Max 39000 sequences per job allowed.</h2>\n";goto fuss}
        
	print "<table bgcolor=#a4d4d6 border=1><tr><td>Seqs submitted:</td><td align=center>$total</td></tr>
        <tr><td>Ignored Seqs (<32AAs):</td><td align=center>",$total-$pregpi,"</td></tr>
        <tr><td>Seqs with C-terminal signal ($anf.log>kohGPI</a>):</td><td align=center>$pos</td></tr>
	<tr><td>$anf.und>Undecidable seqs</a>:</td><td align=center>$und</td></tr>
        <tr><td>$anf.pos.gpi><b>GPI anchored</b></a> (C&N-term signal) 
	($anf.pos.nt.sig>SignalP</a>):</td>\n";
	
        unless ($gpi) {$gpi="O"}
	print "<td align=center><font color=red>$gpi</font></td></tr></table>
	</td><td width=50></td><td>\n";
	&help;
	print "</td></tr></table><map name=koh>\n";
	$px=10;
	open(F,"data/job$id.fasta.imp");
	while (<F>) {print}
	close F;
	
	print "</map><br><table><tr><td>
	<img src=$url/gif/job$id.fasta.png usemap=#koh border=0>\n
	<p><table border=1 bgcolor=#a4d6c1><tr><td><b>Legend:</b></td></tr>
	<tr><td><font color=green>green</font>=GPI anchored, 
	<font color=blue>blue</font>=not GPI anchored,<br>
	<font color=red>red</font>=uncertain, cross=activated by sequence(s)
	</td></tr></table></td>\n";
	
	$px=param("x");$py=param("y");
	if (($px) and ($py)) {
		if ($px eq "O") {$px=0}
		if ($py eq "O") {$py=0}
		print "<td width=50></td><td valign=top><table bgcolor=#a4d6c1 border=1><tr>
		<tr><td><b>Sequences activating unit $px / $py</b></td></tr><td><ul>\n";
		open(G,"data/job$id.fasta.pos.gpi") or die;
		@gpis=<G>;close G;
		open(F,"data/job$id.fasta.log") or die;
		while (<F>) {
			if (/x:(\d+) y:(\d+)/) {
				$x=$1;$y=$2;
				if (($x==$px+1) and ($y==$py+1)) {
					$fnd="";
					for ($i=0;$i<@gpis;$i++) {
						$a=$des;$b=$gpis[$i];
						$a=~s/\s//g;$b=~s/\s//g;
						if ($a eq $b) {$fnd="<font color=green>GPI</font>"}
					}
					$des=~s/^>//;
					print "<li>$des $fnd</li>\n";
				}
			} else {$des=$_}
		}
		close F;print "</ul></td></tr></table></td>\n";
	}
	print "</tr></table>\n";
	&netinfo;
}

sub help {
	print "<table border=1 bgcolor=#49b0f4><tr><td><b><font face=arial><center>
	HELP</center></font></b></td></tr>
	<tr><td width=650><font size=-1 face=arial>
	The <b>table to the left</b> shows the results of your query.
	The <b>first</b> column shows how many sequences you submitted.
	The <b>second</b> tells you how many sequences had to be rejected because they
	were too short to possibly contain even one signal sequence.
	The <b>third</b> column show
	how many sequences have a C-terminal GPI-anchor signal sequence.
	You can click there to access the KohGPI logfile. <b>Fourth</b>,
	there is the number of sequences
	that could not be classified because they could belong to both classes. You can view these
	sequences by clicking on them.
	The <b>fifth</b> and last colums shows how many sequences will
	really be GPI anchored, because they
	have both C- and N-terminal signal sequences. Click there to view a list of these sequences.
	The N-terminal signal sequence was determined by SignalP, whose output can also be seen
	by clicking the link.<br>
	Below this table there is the a <b>graphical representation</b> of the Kohonen
	Self-Organizing Map. Each square symbolizes a unit or simulated neuron.
	According to legend, when a green unit gets activated by a sequence, it has
	a GPI signal sequence. Blue units mean no GPI signal and red ones cannot decide.
	The <b>darker the colour</b> of a square is, the more sequences activated it.
	The squares with crosses in them were activated by one or more of your input sequences.
	You can click on them to see what sequences activated this unit.
	In <b>table to the right of the image</b>,
	that will be opened by this, you see the sequence descriptions with the position
	and quality of the cleavage site added. This information is also contained in the other
	sequence output lists. If the quality is zero, there was no cleavage site found.
	If the is a green <font color=green>GPI</font> next to the description,
	this sequence has a C- as well as a N-terminal
	signal sequence. If not, it only has a C-terminal signal sequence, and this protein
	would not get a GPI anchor in a cell.
	</font></td></tr></table>"
}

sub netinfo {
        open(F,"kohgpi.net") or die;
        $dat=<F>;$dat=<F>;$dat=<F>;$dat=<F>;close F;
        print "<p><font size=-1>Neural Network training date: $dat</font><p>\n";
}

sub navi {  # makes the prev / next links
        $next=$id+1;$prev=$id-1;$src=1;
        open(L,"logs/gpi.log") or die;
        while (<L>) {if (/Job(\d+):.+Jobnam: (.+)$/) {$job[$1]=$2}}
        close L;
        while ((-r "data/job$next.fasta") and ($src)) {
                if ((-r "../html/gif/job$next.fasta.log.gif") and ($job[$next])) {$src=0} else {$next++}
        }
        $src=1;
        while ((-r "data/job$prev.fasta") and ($src)) {
                if ((-r "../html/gif/job$prev.fasta.log.gif") and ($job[$prev])) {$src=0} else {$prev--}
        }
        print "<a href=$url/cgi-bin/gpi.cgi?id=$prev>previous</a> - 
	<a href=$url/cgi-bin/gpi.cgi?id=$next>next</a>\n";
}

sub waiting {
        unless (-r "data/job$id.fasta") {print "Job #$id does not exist!";goto fuss}
	print "<h2>Waiting for Job #$id:</h2><p>\n";
        if (-r "data/job$id.fasta.srf.old") {$phase=1}
        if (-r "data/job$id.fasta.pos") {$phase=2}
        if (-r "data/job$id.fasta.nt.log") {$phase=3}
        unless ($phase) {print "Your Job is on queue... Please be patient.<p>\n"}
        print "<table border=1><tr><td>Sequence PreProcessing</td><td>\n";
        if ($phase>0) {print "<img src=$url/gut.gif>\n"}
        print "</td></tr><tr><td>KohGPI</td><td>\n";
        if ($phase==1) {
                $total=count("data/job$id.fasta.srf");
                $n=count("input.txt.pos")+count("input.txt.neg")+count("input.txt.und");
                print int($n/$total*100),"%";
        }
        if ($phase>1) {print "<img src=$url/gut.gif>\n"}
        print "</td></tr><tr><td>SignalP</td><td>\n";
        if ($phase==2) {print "run"}
        if ($phase>2) {print "<img src=$url/gut.gif>\n"}
        print "</td></tr></table>\n";
}

sub kohonen {
        $i=1;while (-r "data/job$i.fasta") {$i++};$ofn="job$i.fasta";
        $logfile="logs/gpi.log";my $adr=ip();logg("Job$i: $adr Jobnam: $jobnam\n");
        open(O,">data/$ofn") or die;
        if ($fasta) {while(read $fasta,$data,1024) {print O $data}}
        else {if ($seq) {print O $seq}}
        close O;open(O,">data/$ofn.go") or die;print O "Go!";close O;
        $id=$i;$ref=1;
}

sub iform {
        print "<table><tr><td><img src=$url/gpis.gif></td>
        <td width=10></td><td><font face=arial size=+3><b>
        $title</b></font>
        </td></tr><tr><td></td><td></td><td>
	<form action=$url/cgi-bin/gpi.cgi method=post enctype=multipart/form-data>
        <font face=arial>Protein Sequence(s):
        <br><textarea cols=100 rows=15 name=seq wrap=off>
        $seq</textarea><p>
        Upload seqence file: <input type=file name=fasta><p>
        Job name: <input type=text size=32 name=jobnam maxlength=200><p>
        <input type=submit value=GO></form></font>
        <p><font face=arial size=-1>Sequences may be in fasta,
	embl, genbank, SWISS-PROT, gcg, gcgdata, 
        pir or raw format.</font></td></tr></table>\n";
}
   
sub refs {
        print "<a href=$url/cgi-bin/gpi.cgi?arch=1>Job Archive</a><br>
        <a href=$url/gpi_server.png>Server Architecture</a><p>
        <p><b><u>References:</u></b><ul><li><font size=-1>
        $t,<br><a href=http://www.nyk.ch>Nick Fankhauser</a> and
        <a href=http://www.izb.unibe.ch/res/maeser/index.php>Pascal Maeser</a>, Publication
        in preparation</li><li><font size=-1>GPI anchor picture by Haematological Malignancy
        Diagnostic Service (<a href=http://www.hmds.org.uk>www.hmds.org.uk</a>)</font></li></ul>\n";
}
 
sub count {
        open(C,shift) or die;my $cnt=0;
        while (<C>) {if (/^>/) {$cnt++}}
        close C;return $cnt;
}
 
sub ip {
        my $ip=$ENV{'REMOTE_ADDR'};
        system("host $ip >ip.txt");
        open(I,"ip.txt");my $host=<I>;close I;system("rm ip.txt");chomp($host);
        if ($host=~/domain name pointer (.+)\.$/) {$host=$1}
        if ($host eq "") {$host=$ENV{'REMOTE_ADDR'}}
        my $ret="$host ($ip)";return $ret;
}
 
sub logg {
        if (length($_[0])>1) {
                unless (-r $logfile) {open(L,">$logfile");close L;}
                open(L,">>$logfile");print L zeit(),$_[0];close L;
        }
}
 
sub zeit {      # returns current time like 18:37:12 in a string
        ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
         return "$mday.",$mon+1,".",$year+1900," $hour:$min ";
}
