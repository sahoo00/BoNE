#!/usr/bin/perl -I /booleanfs/sahoo/scripts

if (scalar(@ARGV) <= 0) {
    print "perl analyze.pl <cmd> <args>\n";
    exit(1);
}

use U;
use Hegemon;
use JSON;
use Data::Dumper;

my $cmd = shift(@ARGV);
if ($cmd eq "tindex") {
  &analyzeTherapeuticIndex(@ARGV);
}
if ($cmd eq "uc-map") {
  &analyzeUCmap(@ARGV);
}
if ($cmd eq "cd-map") {
  &analyzeCDmap(@ARGV);
}

sub loadJSON {
  my $filename = shift;

  my $json_text = do {
    open(my $json_fh, "<:encoding(UTF-8)", $filename)
      or die("Can't open \$filename\": $!\n");
    local $/;
    <$json_fh>
  };

  my $json = JSON->new;
  my $data = $json->decode($json_text);

  return $data;
}

sub loadGenes {
  my $file = shift;
  my $res = [];
  open(my $fh, "<$file");
  while (<$fh>) {
    chomp;
    my @list = split(/\s/);
    push @$res, @list;
  }
  close($fh);
  return $res;
}

sub analyzeTherapeuticIndex {
  my @list = @_;
  my $shash = &U::getHash("ibd-targets-score.txt");
  my @ids = sort { $shash->{$a}->[1] <=> $shash->{$b}->[1] } keys %{$shash};
  my $ihash = {};
  my $thr = 0.1;
  my $thrNum;
  foreach my $i (0 .. $#ids) {
    $ihash->{$ids[$i]} = $i;
    if (!defined $thrNum && $shash->{$ids[$i]}->[1] > $thr) {
      $thrNum = $i;
    }
  }
  my $m1 =  $shash->{$ids[0]}->[1];
  my $m2 =  $shash->{$ids[$#ids]}->[1];
  my $num = scalar(@ids);
  my $obj = {"Minimum" => $m1, "Maximum" => $m2, "Total Num" => $num,
    "thr" => $thr, "thrNum" => $thrNum};
  foreach my $id (@list) {
    $id = uc($id);
    if (defined $shash->{$id}) {
      $obj->{$id} = [@{$shash->{$id}}, $ihash->{$id}];
    }
  }
  my $fhash = &U::getHash("targets-2.txt");
  foreach my $k (keys %{$fhash}) {
    push @{$obj->{"FDA Approved"}}, 
    [@{$shash->{$k}}, $ihash->{$k}, $fhash->{$k}->[1]];
  }
  my $json = encode_json $obj;
  print "--data--\n$json\n";
}

sub analyzeUCmap {
  my @list = @_;
  my $hash = {};
  foreach my $i (1 .. 5, 7 .. 9) {
    my $file = "UC-info/node-$i\.txt";
    my $l = &U::getList($file);
    foreach my $id (@$l) {
      $hash->{$id} = $i;
    }
  }
  my $obj = {};
  foreach my $id (@list) {
    $id = uc($id);
    if (defined $hash->{$id}) {
      $obj->{$id} = $hash->{$id};
    }
  }
  my $json = encode_json $obj;
  print "--data--\n$json\n";
}

sub analyzeCDmap {
  my @list = @_;
  my $hash = {};
  foreach my $i (1 .. 8) {
    my $file = "CD-info/node-$i\.txt";
    my $l = &U::getList($file);
    foreach my $id (@$l) {
      $hash->{$id} = $i;
    }
  }
  my $obj = {};
  foreach my $id (@list) {
    $id = uc($id);
    if (defined $hash->{$id}) {
      $obj->{$id} = $hash->{$id};
    }
  }
  my $json = encode_json $obj;
  print "--data--\n$json\n";
}

sub analyzeTest3 {
  my $shash = &U::getHash("ibd-targets-score.txt");
  my $fhash = &U::getHash("targets-2.txt");
  foreach my $k (keys %{$fhash}) {
    print join("\t", @{$fhash->{$k}}, @{$shash->{$k}}), "\n";
  }
}
