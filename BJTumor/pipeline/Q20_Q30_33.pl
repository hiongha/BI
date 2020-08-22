print "length\tQ20\tQ30\n";
my $length=0;
my $Q20_all=$Q30_all=0;
my $n=1;
my $N20=$N30=0;
my ($file1)=@ARGV;
open IN1,"zcat $file1 |" or die "can not open IN1";
while(<IN1>){
  if($n==1 or $n==2 or $n==3){$n++;}
  elsif($n==4){
    chomp;
    $length=length($_);
    my @a=split;
    my @b=split(//,$a[0]);
    for(my $m=0;$m<@b;$m++){
      my $score=ord($b[$m])-33;
      if($score>=20){$N20++;$Q20_all++;}
      if($score>=30){$N30++;$Q30_all++;}
    }
    $n=1;
    print length($_)."\t".$N20."\t".$N30."\n";
    $N20=$N30=0;
  }
}
close IN1;
print $length."\t".$Q20_all."\t".$Q30_all."\n";
