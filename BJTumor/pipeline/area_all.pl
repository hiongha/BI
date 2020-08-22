my %hash;
my $seq="";

while(<>){
	chomp;
	my @a=split;
	$seq=$a[4];
	my $cigar=$a[3];#150M90S20M
	my $XF=$a[-1];
	$XF=~s/XF:i://g;
	my $cigar_tmp=$cigar;
	$cigar_tmp=~s/^(\d+)S//;
	$cigar_tmp=~s/(\d+)S//;
	my @chai1=split(/[A-Za-z]/,$cigar_tmp);#150 90	20
	$cigar_tmp=~s/(\d+)//g;
	my @chai2=split(//,$cigar_tmp);#M S M
	my $chr=$a[1];
	my $start=$a[2];
	my $end=0;
	my $check=1;
#	print $_."\n";
	for(my $m=0;$m<@chai2;$m++){
	  if($chai2[$m] eq "M"){
	    if($check==1){
		  $end=$start+$chai1[$m]-1;
		  $check++;
		}
		else{
		  $start=$end+1;
		  $end=$start+$chai1[$m]-1;
		}
		
		
#		for(my $p=$start;$p<=$end;$p++){
#		  my $key=$a[0]."\t".$chr."\t".$p."\t".$XF;
#		  if(exists $hash{$key}){}
#		  else{$hash{$key}=1;}
#		}

		print $chr."\t".$start."\t".$end."\t".$XF."\n";
	  }
	  elsif($chai2[$m] eq "S" or $chai2[$m] eq "D"){
	    $end=$end+$chai1[$m];
	  }
	}
}
#foreach(keys %hash){
#  print $_."\n";
#}
