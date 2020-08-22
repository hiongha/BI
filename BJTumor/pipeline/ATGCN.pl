my $n=1;
my %hash;
my %quality;
while(<>){
  chomp;
  my @a=split;
  if($n==1 or $n==3){$n++;}
  elsif($n==2){
    my @b=split(//,$a[0]);
    for(my $m=0;$m<@b;$m++){
      if(exists $hash{$m}){
        if($b[$m] eq "A"){
          $hash{$m}{"A"}++;
        }
        elsif($b[$m] eq "T"){
          $hash{$m}{"T"}++;
        }
        elsif($b[$m] eq "G"){
          $hash{$m}{"G"}++;
        }
        elsif($b[$m] eq "C"){
          $hash{$m}{"C"}++;
        }
        else{
          $hash{$m}{"N"}++;
        }
      }
      else{
        if($b[$m] eq "A"){
          $hash{$m}{"A"}=1;
        }
        elsif($b[$m] eq "T"){
          $hash{$m}{"T"}=1;
        }
        elsif($b[$m] eq "G"){
          $hash{$m}{"G"}=1;
        }
        elsif($b[$m] eq "C"){
          $hash{$m}{"C"}=1;
        }
        else{
          $hash{$m}{"N"}=1;
        }
      }
    }
    $n++;
  }
  elsif($n==4){
    my @b=split(//,$a[0]);
    for(my $m=0;$m<@b;$m++){
      my $quality_tmp=ord($b[$m])-33;
#      print $quality_tmp."\t";
      if(exists $quality{$quality_tmp}){
        if(exists $quality{$quality_tmp}{$m}){
          $quality{$quality_tmp}{$m}++;
        }
        else{
          $quality{$quality_tmp}{$m}=1;
        }
      }
      else{
        if(exists $quality{$quality_tmp}{$m}){
          $quality{$quality_tmp}{$m}++;
        }
        else{
          $quality{$quality_tmp}{$m}=1;
        }
      }
    }
#    print "\n";
    $n=1;
  }
}

open OUT1,">base.txt" or die "can not output OUT1";
my @a=sort {$a<=>$b} keys %hash;
print OUT1 "pos\t";print OUT1 join("\t",@a);print OUT1 "\n";
my @base=("A","T","G","C","N");
for(my $p=0;$p<@base;$p++){
  print OUT1 $base[$p];
  for(my $m=0;$m<@a;$m++){
    if(exists $hash{$a[$m]}{$base[$p]}){
      print OUT1 "\t".$hash{$a[$m]}{$base[$p]};
    }
    else{
      print OUT1 "\t0";
    }
  }
  print OUT1 "\n";
}


open OUT2,">quality.txt" or die "can not output OUT2";
print OUT2 "pos";
for(my $m=0;$m<@a;$m++){
  print OUT2 "\t".$m;
}
print OUT2 "\n";

for(my $p=0;$p<=50;$p++){
  print OUT2 $p;
  for(my $m=0;$m<@a;$m++){
    if(exists $quality{$p}{$m}){
      print OUT2 "\t".$quality{$p}{$m};
    }
    else{
      print OUT2 "\t0";
    }
  }
  print OUT2 "\n";
}

