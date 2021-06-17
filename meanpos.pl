use GD;$px=10;
$fn=shift;open(F,$fn) or die;print "file opened, counting...\n";@mapp=();
while (<F>) {if (/^seq.+map:G$/) {if (/x:(\d+) y:(\d+) /) {$mapp[int($1)][int($2)]++}}}
for ($x=0;$x<40;$x++) {for ($y=0;$y<40;$y++) {if ($mapp[$x][$y]>$max) {$max=$mapp[$x][$y]}}}
print "max: $max\n";$im=new GD::Image($px*40,$px*40) or die;
$black = $im->colorAllocate(0,0,0);
for ($i=51;$i<256;$i++)  {$red[$i]=$im->colorAllocate($i,0,0)}
for ($x=1;$x<41;$x++) {
	for ($y=1;$y<41;$y++) {
		$c=int(($mapp[$x][$y]/$max)*200)+50;
		for ($i=0;$i<$px-1;$i++) {for ($j=0;$j<$px-1;$j++) {$im->setPixel(($y-1)*$px+$j,($x-1)*$px+$i,$red[$c])}}
	}
}

open PNG,">$fn.gif" or die;binmode PNG;
print PNG $im->gif;close PNG;print "$fn.gif written\n";
