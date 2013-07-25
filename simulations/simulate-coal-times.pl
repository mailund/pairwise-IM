
use strict;

sub simulate_coal_times;

my $pop1size = 1;
my $pop2size = 2;
my $m12 = 0.1;
my $m21 = 0.2;

my $change = 0.5;
my $m12t = 0.4;
my $m21t = 0.8;
my $pop1size2 = 2;
my $pop2size2 = 1;

# no change...
#my $m12t = 0.2;
#my $m21t = 0.2;
#my $pop1size2 = 1;
#my $pop2size2 = 2;

my $repeats = 500;

my $sample_conf;
my $mscmd;

for $sample_conf ('2 0', '1 1', '0 2') {
  $mscmd = "./ms 2 $repeats -T -I 2 $sample_conf " .
                 " -m 1 2 $m12 -m 2 1 $m21 -n 1 $pop1size -n 2 $pop2size " .
                 " -em $change 1 2 $m12t -em $change 2 1 $m21t " .
	         " -en $change 2 $pop2size2 -en $change 1 $pop1size2 ";
  print $mscmd, "\n";
  my $printable_conf = $sample_conf;
  $printable_conf =~ tr/ /-/;

  open OUT, ">coal.$printable_conf.txt";

  for my $ct (simulate_coal_times()) {
    # Multiply by 2 to get to a time unit of 2N rather than 4N
    print OUT $ct * 2, "\n";
  }
}



sub simulate_coal_times {
  my @coal_times;
  for my $line (`$mscmd`) {
    # This heavily depends on there being a tree with two samples!
    if ($line =~ /\(1:(\d+.\d+),/) {
      push @coal_times, $1; 
    }
  }

  return @coal_times;
}
