@files = <*.xml>;

foreach $file_model (@files)
{
	open(RP, "pcacoeff.bin");

	open(R, "$file_model");
	open(W, ">${file_model}_");

	binmode RP;
	binmode W;

	while(<R>){
		$tmp = $_;
		print W $tmp;
		if($tmp =~ /<P>/)
		{
		  ($s, $p, $f) = split(/\<\/?P\>/, $tmp);
		  print "$p";
		  last;
		}
	}

	print W "<PCAcoeff>";

	read (RP, $tmp, 961* 8);
	print W $tmp;
	print W "\n</PCAcoeff>\n";

	binmode R;

	while(read (R, $tmp, 1024)){
		print W $tmp;
	}
	close(W);
	close(R);
	close(RP);
}