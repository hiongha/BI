my %hash;
my ($file1,$file2)=@ARGV;
open IN1,"<$file1" or die "can not open IN1";
while(<IN1>){
  chomp;
  my @a=split(/:/,$_);
  my @b=split(/-/,$a[1]);
  for(my $n=$b[0];$n<=$b[1];$n++){
    my $key=$a[0]."\t".$n;
	$hash{$key}=1;
  }
}
close IN1;

open IN2,"<$file2" or die "can not open IN2";
while(<IN2>){
  chomp;
  my @a=split;
  my $key=$a[1]."\t".$a[2];
  if(exists $hash{$key}){
    print $_."\n";
  }
}
close IN2;