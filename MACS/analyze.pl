#!/usr/bin/perl -I /booleanfs/sahoo/scripts

if (scalar(@ARGV) <= 0) {
    print "perl analyze.pl <cmd> <args>\n";
    exit(1);
}

use IPC::Open2;
use U;
use Hegemon;
use Data::Dumper;
use JSON;

my $cmd = shift(@ARGV);
if ($cmd eq "dynrange") {
  &analyzeDynRange(@ARGV);
}
if ($cmd eq "select-noad") {
  &analyzeNOAD(@ARGV);
}
if ($cmd eq "network") {
  &analyzeNetwork(@ARGV);
}
if ($cmd eq "figure-1") {
  &analyzeFigure1(@ARGV);
}
if ($cmd eq "surv-1") {
  &analyzeSurv1(@ARGV);
}
if ($cmd eq "surv-2") {
  &analyzeSurv2(@ARGV);
}

sub analyzeDynRange {
  my $pre = "/booleanfs2/sahoo/Data/Colon/Jung/gpl570-colon";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
    survival => "$pre\-survival.txt");
  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $carcinoma = $h->getSurvivalArray(9, "Carcinoma");
  my $adenoma = $h->getSurvivalArray(9, "Adenoma");
  my $normal = $h->getSurvivalArray(9, "Normal");
  my $series = $h->getSurvivalData(4);
  my $type = $h->getSurvivalName('c Histology');

  my $chash = &U::getHash("cec.txt");
  my $cec = [ map { $h->{'hhash'}->{$_} } keys %{$chash} ];
  $normal = &U::diff($normal, $cec);

  my $pG = [ ["All", "red", [ $start .. $end ] ] ];
  my $pG = [ ["Carcinoma", "#00EEEE", $carcinoma],
    ["Adenoma", "#EE1000", $adenoma],
    ["Normal", "#000000", $normal]];
  my $f = "$pre\-vinfo.txt";
  open(my $fh, "<$f") || die "Can't open $f\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  my $hhash = {};
  for (my $i = 0; $i < scalar(@headers); $i++) {
    $hhash->{$headers[$i]} = $i;
  }
  my $dhash = {};
  my $data = [];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $dyn = $list[$hhash->{"Perc 0.90"}] - $list[$hhash->{"Perc -0.90"}];
    $dhash->{$list[0]} = $dyn;
    push @$data, $dyn;
  }
  close($fh);
  $data = [sort { $a <=> $b } @$data];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  print join("\t", @$thr), "\n";
  #my $ddata = [grep { $_ > $thr->[0] } @$data];
  #$ddata = [sort { $a <=> $b } @$ddata];
  #my $res = &U::fitStep($ddata, 0, scalar(@$ddata) - 1);
  #my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  #print join("\t", @$thr), "\n";
  my $res = &histogramArray("figures/dynrange-stepminer.png", $data, 
    thr => $thr->[0]);

  my @ids = keys(%{$h->{"idHash"}});
  my $nhash = {};
  foreach my $id (@ids) {
    next if (defined $dhash->{$id} && $dhash->{$id} < $thr->[0]);
    my $name = $h->getName($id);
    next if ($name eq "---");
    $name = (split(" /// ", $name))[0];
    push @{$nhash->{$name}}, [$id, $dhash->{$id}];
  }
  my $finalIDs = [];
  open(my $ofh, ">colon-dyn-select.txt");
  print $ofh join("\t", "ArrayID", "DynRange", "Name"), "\n";
  foreach my $name (keys %{$nhash}) {
    my $l = [sort { $b->[1] <=> $a->[1] } @{$nhash->{$name}}];
    my ($id, $dyn) = @{$l->[0]};
    push @$finalIDs, $id;
    print $ofh join("\t", $id, $dyn, $name, map { $_->[1] } @$l), "\n";
  }
  close($ofh);

  open(my $ofh, ">colon-dyn-expr.txt");
  my $hdr = $h->{'headers'};
  print $ofh join("\t", @$hdr), "\n";
  my $rows = [];
  my $expr = [];
  my $columns = [ map { $hdr->[$_] } 2 .. $#{$hdr}];
  foreach my $id (@$finalIDs) {
    next if (defined $dhash->{$id} && $dhash->{$id} < $thr->[0]);
    my $name = $h->getName($id);
    $name = (split(" /// ", $name))[0];
    my $e = $h->getExprData($id);
    $e->[0] = $name;
    $e->[1] = $id;
    print $ofh join("\t", @$e), "\n";
    push @$rows, $id;
    push @$expr, [ map { $e->[$_] } 2 .. $#{$hdr}];
  }
  close($ofh);
  open(my $ofh, ">colon-noad-expr.txt");
  my $hdr = $h->{'headers'};
  print $ofh join("\t", map { $hdr->[$_] } (0, 1, @$normal, @$adenoma)), "\n";
  my $rows = [];
  my $expr = [];
  my $columns = [ map { $hdr->[$_] } 2 .. $#{$hdr}];
  foreach my $id (@$finalIDs) {
    next if (defined $dhash->{$id} && $dhash->{$id} < $thr->[0]);
    my $name = $h->getName($id);
    $name = (split(" /// ", $name))[0];
    my $e = $h->getExprData($id);
    $e->[0] = $name;
    $e->[1] = $id;
    print $ofh join("\t", map { $e->[$_] } (0, 1, @$normal, @$adenoma)), "\n";
    push @$rows, $id;
    push @$expr, [ map { $e->[$_] } 2 .. $#{$hdr}];
  }
  close($ofh);
  $h->subset("normal-adenoma", [@$normal, @$adenoma]);
}

sub analyzeNOAD {
  my $pre = "normal-adenoma";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
    survival => "$pre\-survival.txt");
  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $f = "$pre\-vinfo.txt";
  open(my $fh, "<$f") || die "Can't open $f\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  my $hhash = {};
  for (my $i = 0; $i < scalar(@headers); $i++) {
    $hhash->{$headers[$i]} = $i;
  }
  my $dhash = {};
  my $data = [];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $dyn = $list[$hhash->{"Perc 0.90"}] - $list[$hhash->{"Perc -0.90"}];
    my $name = $h->getName($list[0]);
    if ($name =~ /CCDC88A/) {
      print join("\t", $list[0], $dyn, $name), "\n";
    }
    if ($name =~ /CDX/) {
      print join("\t", $list[0], $dyn, $name), "\n";
    }
    $dhash->{$list[0]} = $dyn;
    push @$data, $dyn;
  }
  close($fh);
  $data = [sort { $a <=> $b } @$data];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  print join("\t", @$thr), "\n";
  return;
  #my $ddata = [grep { $_ > $thr->[0] } @$data];
  #$ddata = [sort { $a <=> $b } @$ddata];
  #my $res = &U::fitStep($ddata, 0, scalar(@$ddata) - 1);
  #my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  #print join("\t", @$thr), "\n";
  my $res = &histogramArray("figures/noad-stepminer.png", $data, 
    thr => $thr->[0]);

  my @ids = keys(%{$h->{"idHash"}});
  my $nhash = {};
  foreach my $id (@ids) {
    next if (defined $dhash->{$id} && $dhash->{$id} < $thr->[2]);
    my $name = $h->getName($id);
    next if ($name eq "---");
    $name = (split(" /// ", $name))[0];
    push @{$nhash->{$name}}, [$id, $dhash->{$id}];
  }
  my $finalIDs = [];
  open(my $ofh, ">colon-noad-select.txt");
  print $ofh join("\t", "ArrayID", "DynRange", "Name"), "\n";
  foreach my $name (keys %{$nhash}) {
    my $l = [sort { $b->[1] <=> $a->[1] } @{$nhash->{$name}}];
    my ($id, $dyn) = @{$l->[0]};
    push @$finalIDs, $id;
    print $ofh join("\t", $id, $dyn, $name, map { $_->[1] } @$l), "\n";
  }
  close($ofh);

  open(my $ofh, ">colon-noad-expr.txt");
  my $hdr = $h->{'headers'};
  print $ofh join("\t", map { $hdr->[$_] } (0 .. $end)), "\n";
  my $rows = [];
  my $expr = [];
  my $columns = [ map { $hdr->[$_] } 2 .. $#{$hdr}];
  foreach my $id (@$finalIDs) {
    next if (defined $dhash->{$id} && $dhash->{$id} < $thr->[2]);
    my $name = $h->getName($id);
    $name = (split(" /// ", $name))[0];
    my $e = $h->getExprData($id);
    $e->[0] = $name;
    $e->[1] = $id;
    print $ofh join("\t", map { $e->[$_] } (0 .. $end)), "\n";
    push @$rows, $id;
    push @$expr, [ map { $e->[$_] } 2 .. $#{$hdr}];
  }
  close($ofh);
}

sub histogramArray {
  my ($ofile, $arr, %params) = @_;
use IPC::Open2;
  my $break = 20;
  $break = $params{'break'} if (defined $params{'break'});
  my $width = 640;
  my $height = 480;
  my $dpi = 100;
  $width = $params{'width'} if (defined $params{'width'});
  $height = $params{'height'} if (defined $params{'height'});
  $dpi = $params{'dpi'} if (defined $params{'dpi'});
  my $wd = $width/$dpi;
  my $ht = $height/$dpi;
  $wd = $params{'wd'} if (defined $params{'wd'});
  $ht = $params{'ht'} if (defined $params{'ht'});
  my ($phr, $phw);
  if (defined $params{'ph'}) {
    ($phr, $phw) = @{$params{'ph'}};
  }
  else {
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  my $devstr = "";
  if (defined $ofile) {
    my @l = split("\\.", $ofile);
    if ($l[$#l] eq "png") {
      $devstr = "png(filename=\"$ofile\", width=$width, height=$height, pointsize=15)";
    }
    if ($l[$#l] eq "ps") {
      $devstr = "postscript(file=\"$ofile\", width=$wd, height=$ht, pointsize=11)";
    }
    if ($l[$#l] eq "pdf") {
      $devstr = "pdf(file=\"$ofile\", width=$wd, height=$ht, pointsize=11)";
    }
  }
  print $phw "l <- c(", join(",", @$arr), ");\n";
  print $phw "$devstr\n";
  my $main = "$ofile";
  if (defined $params{'title'}) {
    $main = $params{'title'};
  }
  my $pa = "";
  if (defined $params{'pa'}) {
    $pa = $params{'pa'};
  }
  if ($params{'logy'}) {
print $phw <<END
h <- hist(l, breaks=$break, freq=TRUE, main="$main", xlab="value", ylab="count")
n <- length(h\$counts)
plot(h\$breaks[1:n], log10(h\$counts), main="$main", xlab="value", ylab="log10(count)"$pa)
cat(h\$breaks, "\\n")
cat(h\$counts, "\\n")
dev.off()
END
;
  }
  else {
print $phw <<END
h <- hist(l, breaks=$break, freq=TRUE, main="$main", xlab="value", ylab="count")
n <- length(h\$counts)
plot(h\$breaks[1:n], h\$counts, main="$main", xlab="value", ylab="count"$pa)
cat(h\$breaks, "\\n")
cat(h\$counts, "\\n")
dev.off()
END
;
  }
  close($phw);
  my $breaks = <$phr>;
  my $counts = <$phr>;
  close($phr);
  $breaks =~ s/[\r\n]//g;
  $counts =~ s/[\r\n]//g;
  my @breaks = split(" ", $breaks);
  my @counts = split(" ", $counts);
  for (my $i = 0; $i < scalar(@counts); $i++) {
    my $l = "";
    if ($counts[$i] > 0) {
      $l = log($counts[$i])/log(10);
    }
    print join("\t", $breaks[$i], $breaks[$i+1], $counts[$i], $l), "\n";
  }
  return [\@breaks, \@counts];
}

sub gamma {
  my ($u, $edges, $hash) = @_;
  if (defined $hash && defined $hash->{$u}) {
    return $hash->{$u};
  }
  my $res = [$u, keys %{$edges->{$u}}];
  if (defined $hash) {
    $hash->{$u} = $res;
  }
  return $res;
}
sub rho {
  my ($u, $v, $edges, $hash, $shash) = @_;
  if (defined $shash && defined $shash->{$u}->{$v}) {
    return $shash->{$u}->{$v};
  }
  my $gu = &gamma($u, $edges, $hash);
  my $gv = &gamma($v, $edges, $hash);
  my $g_union = &U::union($gu, $gv);
  my $g_int = &U::intersection($gu, $gv);
  my $res = scalar(@$g_int) / scalar(@$g_union);
  if (defined $shash) {
    $shash->{$u}->{$v} = $res;
  }
  return $res;
}

sub convertCodes {
  my $vcode = shift;
  my $b0 = $vcode->bit_test(1);
  my $b1 = $vcode->bit_test(0);
  my $b2 = $vcode->bit_test(3);
  my $b3 = $vcode->bit_test(2);
  my $total = $b0 + $b1 + $b2 + $b3;
  if ($total == 1) {
    if ($b0 == 1) { return 1; }
    if ($b1 == 1) { return 2; }
    if ($b2 == 1) { return 3; }
    if ($b3 == 1) { return 4; }
  }
  if ($total == 2) {
    if ($b1 == 1 && $b2 == 1) { return 5; }
    if ($b0 == 1 && $b3 == 1) { return 6; }
  }
  return 0;
}

sub getNCodes{
  my ($n, $net, $l, $balanced, $balancedhash) = @_;
  my $num = $n->getNum();
  my $relations = {};
  my $nb = $n->{'numBits'};
  my $vcode = Bit::Vector->new($nb*2*2);
  foreach my $u (@$l) {
    my $i = $balancedhash->{$u};
    for (my $j = 0; $j < $num/2; $j++) {
      $n->readNCode($net, $i, $j, $vcode);
      my $code = &convertCodes($vcode);
      next if ($code <= 0);
      if (!defined $relations->{$code}) {
        $relations->{$code} = {};
      }
      $relations->{$code}->{$balanced->[$j]}++;
    }
  }
  return $relations;
}

sub getNCodesPair {
  my ($n, $net, $u, $v, $balanced, $balancedhash) = @_;
  my $num = $n->getNum();
  my $relations = {};
  my $nb = $n->{'numBits'};
  my $vcode = Bit::Vector->new($nb*2*2);
  my $i = $balancedhash->{$u};
  my $j = $balancedhash->{$v};
  $n->readNCode($net, $i, $j, $vcode);
  my $code = &convertCodes($vcode);
  print STDERR "$u $i $v $j $code ". $vcode->to_Bin() . "\n";
  if ($code > 0) {
    if (!defined $relations->{$code}) {
      $relations->{$code} = {};
    }
    $relations->{$code}->{$balanced->[$j]} = 1;
  }
  my $code = $n->readCode_1_1($i, $j);
  print STDERR "$u $i $v $j $code\n";
  return $relations;
}

sub getBooleanPath {
  my ($g, $nodes, $edges, $start, $rank, $ehash, $pathway, $nodelist) = @_;
  my $path = &Graph::getDeepPath($nodes, $edges, $start, $rank);
  my $index = 0;
  my $ihash = {};
  my $ids = [];
  use UnionFind;
  my $uf = UnionFind->new();
  foreach my $id (map { $g->getNode($_) } @$path) {
    last if (defined $pathway->{$id});
    $pathway->{$id} = 1;
    my @k = keys %{$ehash->{$id}};
    foreach my $i (@k) {
      push @{$ihash->{$i}}, $index;
    }
    push @$ids, $id;
    $uf->makeset($id);
    $index++;
  }
  foreach my $i (keys %{$ihash}) {
    my $s = &U::min($ihash->{$i});
    my $e = &U::max($ihash->{$i});
    for (my $j = $s + 1; $j <= $e; $j ++) {
      $uf->union($ids->[$j], $ids->[$j-1]);
    }
  }
  my $khash = {};
  foreach my $i (keys %{$ihash}) {
    my $s = &U::min($ihash->{$i});
    my $e = &U::max($ihash->{$i});
    for (my $j = $s; $j <= $e; $j ++) {
      my $p = $uf->findset($ids->[$j]);
      $khash->{$p}->{$i} = 1;
    }
  }
  #print Dumper($ids);
  #print Dumper($ihash);
  #print Dumper($khash);
  my $clusters = {};
  my $bpath = [];
  foreach my $n (@$ids) {
    my $p = $uf->findset($n);
    if (!defined $clusters->{$p}) {
      push @$bpath, $p;
    }
    push @{$clusters->{$p}}, [$n, scalar(@{$nodelist->{$n}})];
  }
  my $res = [map { [$clusters->{$_}, [keys %{$khash->{$_}}]] } @$bpath];
  #print Dumper($res);
  return $res;
}

sub cleanupBooleanPath {
  my ($g, $path, $nodelist, $ehash) = @_;
  my ($nodes, $edges) = $g->getNodesEdges("5");
  my $kh = {};
  my $res = [];
  my $ekh = {};
  foreach my $k (@$path) {
    my $arr2 = [];
    foreach my $k2 (@{$k->[1]}) {
      next if (defined $ekh->{$k2});
      $ekh->{$k2} = 1;
      push @$arr2, $k2;
    }
    my $nl = [];
    foreach my $k1 (@{$k->[0]}) {
      my $ki = $g->getIndex($k1->[0]);
      next if (defined $kh->{$ki});
      push @$nl, $k1;
      $kh->{$ki} = 1;
    }
    next if (scalar(@$nl) == 0);
    my $extra = [];
    foreach my $k1 (@$nl) {
      my $i = $g->getIndex($k1->[0]);
      foreach my $ki (keys %{$edges->{$i}}) {
        next if (defined $kh->{$ki});
        my $n = $g->getNode($ki);
        next if (!defined $nodelist->{$n});
        if (defined $ehash->{$n}) {
          foreach my $k2 (keys %{$ehash->{$n}}) {
            next if (defined $ekh->{$k2});
            $ekh->{$k2} = 1;
            push @$arr2, $k2;
          }
        }
        push @$extra, [$n, scalar(@{$nodelist->{$n}})];
        $kh->{$ki} = 1;
      }
    }
    push @$res, [ [@$nl, @$extra], $arr2 ];
  }
  return $res;
}

sub analyzeNetwork {
  my @params = @_;
  if (defined $params[0] && $params[0] eq "eq") {
    use Graph;
    my $g = Graph->new();
    my $file = "colon-network-res.txt";
    my $count = 0;
    open(my $fh, "<$file");
    while (<$fh>) {
      next if (!/^Found : 5/);
      chomp;
      my @list = split("\t");
      my $id1 = $list[3];
      my $id2 = $list[4];
      $g->addEdge($id1, $id2, "eq");
      $count++;
      #if ($count > 10000) {
      #  last;
      #}
    }
    close($fh);
    my ($nodes, $edges) = $g->getNodesEdges("eq");
    my $hash = {};
    my $shash = {};
    use UnionFind;
    my $uf = UnionFind->new();
    my $count = 0;
    my $num = scalar(keys %{$edges});
    foreach my $u (keys %{$edges}) {
      print STDERR "$count/$num\n";
      my $data = [];
      my $scores = [];
      foreach my $v (keys %{$edges->{$u}}) {
        my $r = &rho($u, $v, $edges, $hash, $shash);
        push @$data, $r;
        push @$scores, [$r, $v];
        $uf->makeset($u);
        $uf->makeset($v);
      }
      my $filter = [sort { $b->[0] <=> $a->[0] } @$scores];
      foreach my $s (@$filter) {
        if ($uf->findset($u) != $uf->findset($s->[1])) {
          $uf->union($u, $s->[1]);
          print join("\t", $g->getNode($u), $s->[0], 
            $g->getNode($s->[1])), "\n";
          last;
        }
      }
      $count++;
    }
  }
  if (defined $params[0] && $params[0] eq "eq-rho-thr") {
    my $file = "colon-network-eq.txt";
    my $index = 0;
    open(my $fh, "<$file");
    my $data = [];
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      if ($list[1] > 0) {
        push @$data, $list[1];
      }
      $index++;
    }
    close($fh);
    $data = [sort { $a <=> $b } @$data];
    my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
    my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
    print STDERR join("\t", "DR", @$thr), "\n";
  }
  if (defined $params[0] && $params[0] eq "cls") {
    use UnionFind;
    my $uf = UnionFind->new();
    my $nodes = {};
    use Network;
    my $file = "colon-network-1.rl";
    my $n = Network->new(file => $file, mode => "r");
    $n->readNetworkFile();
    my $balanced = $n->getBalanced();
    my $balancedhash = {};
    for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
      $balancedhash->{$balanced->[$i]} = $i;
      my $id = $balanced->[$i];
      $uf->makeset($id);
      $nodes->{$id} = {} if (!defined $nodes->{$id});
      $nodes->{$id}->{$id} = 1;
    }
    my $num = $n->getNum();
    my $rhash = {};
    my $net = $n->readAll();
    my $file = "colon-network-eq.txt";
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id1 = $list[0];
      my $id2 = $list[2];
      $uf->makeset($id1);
      $nodes->{$id1} = {} if (!defined $nodes->{$id1});
      $nodes->{$id1}->{$id1} = $list[1];
      $uf->makeset($id2);
      $nodes->{$id2} = {} if (!defined $nodes->{$id2});
      $nodes->{$id2}->{$id2} = $list[1];
      if ($list[1] > 0.5) {
        $uf->union($id1, $id2);
        $nodes->{$id1}->{$id2} = $list[1];
        $nodes->{$id2}->{$id1} = $list[1];
      }
    }
    close($fh);
    my $rank = {};
    foreach my $n (keys %{$nodes}) {
      $rank->{$n} = scalar(keys %{$nodes->{$n}});
    }
    my $clusters = {};
    foreach my $n (keys %{$nodes}) {
      push @{$clusters->{$uf->findset($n)}}, $n;
    }
    if (defined $params[1]) {
      my $n = $uf->findset($params[1]);
      my @l = sort { $rank->{$b} <=> $rank->{$a} } @{$clusters->{$n}};
      print scalar(@l), "\n";
      foreach my $id (@l) {
        print join("\t", $id, $rank->{$id}), "\n";
      }
    }
    else {
      foreach my $n (keys %{$clusters}) {
        my @l = sort { $rank->{$b} <=> $rank->{$a} } @{$clusters->{$n}};
        print join("\t", $l[0], scalar(@l), @l), "\n";
      }
    }
  }
  if (defined $params[0] && $params[0] eq "g-cls") {
    my $file = "colon-network-cls.txt";
    my $nodes = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id = $list[0];
      my $l = [@list[2 .. $#list]];
      $nodes->{$id} = $l;
      push @$ids, $id;
    }
    close($fh);
    use Network;
    my $file = "colon-network-1.rl";
    my $n = Network->new(file => $file, mode => "r");
    $n->readNetworkFile();
    my $balanced = $n->getBalanced();
    my $balancedhash = {};
    for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
      $balancedhash->{$balanced->[$i]} = $i;
    }
    my $num = $n->getNum();
    my $rhash = {};
    my $net = $n->readAll();
    foreach my $i (0 .. $#{$ids}) {
      my $u = $ids->[$i];
      if (defined $params[1]) {
        $u = $params[1];
      }
      print STDERR "$i $u\n";
      my $n1 = scalar(@{$nodes->{$u}});
      my $mid = int($n1 / 2);
      my $l = [$u];
      if ($n1 > 2) {
        $l = [ $u, map { $nodes->{$u}->[$_] } (1, $mid)];
      }
      if ($n1 > 10) {
        push @$l, map { $nodes->{$u}->[$_] } (int($mid/4), int($mid/2), $mid - 1);
      }
      my $ru = &getNCodes($n, $net, $l, $balanced, $balancedhash);
      if (defined $params[2]) {
        my $v = $params[2];
        print STDERR $v, " ", scalar(@$l), " $n1\n";
        foreach my $u (@$l) {
          print STDERR $u, "\n";
          &getNCodesPair($n, $net, $u, $v, $balanced, $balancedhash);
        }
      }
      foreach my $c (1 .. 6) {
        foreach my $v (keys %{$ru->{$c}}) {
          next if (!defined $nodes->{$v});
          next if ($n1 > 10 && $ru->{$c}->{$v} < 3);
          next if (defined $params[2] && $params[2] ne $v);
          print join("\t", $u, $c, $v, $ru->{$c}->{$v}), "\n";
        }
      }
      if (defined $params[1]) {
        last;
      }
    }
  }
  if (defined $params[0] && $params[0] eq "genes") {
    use Graph;
    my $g = Graph->new();
    my $file = "colon-network-eq-g.txt";
    my $nodelist = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      $g->addEdge($list[0], $list[2], $list[1]);
    }
    close($fh);
    my $file = "colon-network-cls.txt";
    my $nhash = {};
    my $rank = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id = $list[0];
      my $l = [@list[2 .. $#list]];
      $nodelist->{$id} = $l;
      $nhash->{$id} = $id;
      foreach my $i (@$l) {
        $nhash->{$i} = $id;
      }
      $rank->{$g->getIndex($id)} = scalar(@$l);
      push @$ids, $g->getIndex($id);
    }
    close($fh);
    my $ehash = {};
    my $ehash1 = {};
    my ($nodes, $edges) = $g->getNodesEdges("5");
    if (defined $params[1]) {
      foreach my $id (@params) {
        next if (!defined $nhash->{$id});
        my $i = $g->getIndex($nhash->{$id});
        if (defined $ehash1->{$g->getNode($i)}) {
          $rank->{$i} += 500;
        }
        $ehash1->{$g->getNode($i)}->{$id} = 1;
        $ehash->{$g->getNode($i)}->{$id} = 1;
        if ($rank->{$i} < 2000) {
          $rank->{$i} += 4000;
        }
        foreach my $k (keys %{$edges->{$i}}) {
          if ($rank->{$k} > 5 && $rank->{$k} < 2000) {
            $rank->{$k} += 2000;
          }
          $ehash->{$g->getNode($k)}->{$id} = 1;
        }
      }
    }
    #print Dumper($ehash);
    $ids = [sort { $rank->{$b} <=>  $rank->{$a} } @$ids];
    my $res = [];
    push @$res, [0, scalar(@$ids), 
      [map { [$_, $rank->{$g->getIndex($_)}] } 
      map { $g->getNode($ids->[$_]) } 0 .. 5], []];
    my $pathway = {};
    my ($nodes, $edges) = $g->getNodesEdges("2");
    my $start = $ids->[0];
    my $path = &getBooleanPath($g, $nodes, $edges, $start, $rank, $ehash, $pathway, $nodelist);
    my $path = &cleanupBooleanPath($g, $path, $nodelist, $ehash);
    my $index = 1;
    foreach my $k (@$path) {
      push @$res, [-$index, scalar(@{$k->[0]}), @$k];
      $index++;
    }
    my ($nodes, $edges6) = $g->getNodesEdges("6");
    my ($nodes, $edges4) = $g->getNodesEdges("4");
    my $edges = &Graph::doUnion($edges6, $edges4);
    my $ids = [sort { $rank->{$b} <=> $rank->{$a} } keys %{$edges->{$start}}];
    push @$res, [0, scalar(@$ids), 
      [map { [$_, $rank->{$g->getIndex($_)}] } 
      map { $g->getNode($ids->[$_]) } 0 .. 5], []];
    my ($nodes, $edges) = $g->getNodesEdges("2");
    my $start = $ids->[0];
    my $path = &getBooleanPath($g, $nodes, $edges, $start, $rank, $ehash, $pathway, $nodelist);
    my $path = &cleanupBooleanPath($g, $path, $nodelist, $ehash);
    my $index = 1;
    foreach my $k (@$path) {
      push @$res, [$index, scalar(@{$k->[0]}), @$k];
      $index++;
    }

    my $json = encode_json $res;
    print $json, "\n";
    #foreach (@$res) {
    #  print join("\t", @$_), "\n";
    #}
  }
  if (defined $params[0] && $params[0] eq "query") {
    my $file = "colon-network-cls.txt";
    my $nodelist = {};
    my $nhash = {};
    my $rank = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id = $list[0];
      my $l = [@list[2 .. $#list]];
      $nodelist->{$id} = $l;
      $nhash->{$id} = $id;
      foreach my $i (@$l) {
        $nhash->{$i} = $id;
      }
      $rank->{$id} = scalar(@$l);
      push @$ids, $id;
    }
    close($fh);
    use Network;
    my $file = "colon-network-1.rl";
    my $n = Network->new(file => $file, mode => "r");
    $n->readNetworkFile();
    my $balanced = $n->getBalanced();
    my $balancedhash = {};
    for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
      $balancedhash->{$balanced->[$i]} = $i;
    }
    my $num = $n->getNum();
    my $rhash = {};
    my $net = $n->readAll();
    my $u = $nhash->{$params[1]};
    my $n1 = scalar(@{$nodelist->{$u}});
    print "$u ($n1)\n";
    my $l = $nodelist->{$u};
    my $v = $nhash->{$params[2]};
    my $n2 = scalar(@{$nodelist->{$v}});
    print "$v ($n2)\n";
    my $ru = &getNCodes($n, $net, $l, $balanced, $balancedhash);
    my $vhash = {};
    foreach my $v1 (@{$nodelist->{$v}}) {
      $vhash->{$v1} = 1;
    }
    foreach my $c (1 .. 6) {
      foreach my $v (keys %{$ru->{$c}}) {
        next if (!defined $vhash->{$v});
        print join("\t", $u, $c, $v, $ru->{$c}->{$v}), "\n";
      }
    }
  }
  if (defined $params[0] && $params[0] eq "path") {
    use Graph;
    my $g = Graph->new();
    my $file = "colon-network-eq-g.txt";
    my $nodelist = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      $g->addEdge($list[0], $list[2], $list[1]);
    }
    close($fh);
    my $file = "colon-network-cls.txt";
    my $nhash = {};
    my $rank = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id = $list[0];
      my $l = [@list[2 .. $#list]];
      $nodelist->{$id} = $l;
      $nhash->{$id} = $id;
      foreach my $i (@$l) {
        $nhash->{$i} = $id;
      }
      $rank->{$g->getIndex($id)} = scalar(@$l);
      push @$ids, $g->getIndex($id);
    }
    close($fh);
    my $u = $nhash->{$params[1]};
    my $n1 = scalar(@{$nodelist->{$u}});
    print "$u ($n1)\n";
    my $l = $nodelist->{$u};
    my $v = $nhash->{$params[2]};
    my $n2 = scalar(@{$nodelist->{$v}});
    print "$v ($n2)\n";
    my $res = $g->Path("2", $u, $v, 0, 100);
    print join("\t", map { $g->getNode($res->[$_])." ".
     $g->getType($res->[$_], $res->[$_ + 1]) } 0 .. $#{$res}), "\n";
  }
  if (defined $params[0] && $params[0] eq "state") {
    my $file = "colon-network-cls.txt";
    my $nodelist = {};
    my $nhash = {};
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id = $list[0];
      my $l = [@list[2 .. $#list]];
      $nodelist->{$id} = $l;
      $nhash->{$id} = $id;
      foreach my $i (@$l) {
        $nhash->{$i} = $id;
      }
    }
    close($fh);
    my $pre = "/booleanfs2/sahoo/Data/Colon/Jung/gpl570-colon";
    my $pre = "colon-noad";
    my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => "$pre\-survival.txt");
    my $start = $h->getStart();
    my $end = $h->getEnd();
    my $carcinoma = $h->getSurvivalArray(9, "Carcinoma");
    my $adenoma = $h->getSurvivalArray(9, "Adenoma");
    my $normal = $h->getSurvivalArray(9, "Normal");
    my $series = $h->getSurvivalData(4);
    my $type = $h->getSurvivalName('c Histology');

    my $pG = [ ["All", "red", [ $start .. $end ] ] ];
    my $pG = [
      ["Adenoma", "#EE1000", $adenoma],
      ["Normal", "#000000", $normal]];
    #foreach my $g (@$pG) {
    #  print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
    #}
    foreach my $id (keys %{$nodelist}) {
      my $e = $h->getExprData($id);
      my $eG = [];
      foreach my $g (@$pG) {
        my $m = &U::mean([map { $e->[$_] } @{$g->[2]}]);
        push @$eG, $m;
      }
      my $l = [sort { $eG->[$b] <=> $eG->[$a] } 0 .. 1];
      print join("\t", $id, scalar(@{$nodelist->{$id}}), 
        $pG->[$l->[0]]->[0], @$eG), "\n";
    }
  }
}

sub analyzeFigure1 {
  my @params = @_;
  if (defined $params[0] && $params[0] eq "clusters") {
    my $file = "colon-network-cls.txt";
    my $data = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      next if ($list[1] == 0);
      push @$data, $list[1];
    }
    close($fh);
    $data = [ sort { $b <=> $a } @$data ];
    use Statistics::R;
    my $R = Statistics::R->new();
    $devstr = "pdf(file=\"figures/cluster-sizes.pdf\", width=6, height=4, pointsize=10)";
    $R->run($devstr);
    $R->run("y <- c(". join(",", @$data). ")");
    $R->run("x <- c(". join(",", 1 .. scalar(@$data)). ")");
      my $str = <<END
lo <- data.frame(x=x, y=y)
plot(lo, type="l", col=1, main="Cluster Sizes", xlab="Index", 
    ylab="Count", log="xy")
#points(lo, pch=19, col="cyan")
#with(lo, symbols(x=x, y=y, circles=log(y)/log(2), inches=1/10,
#                      add=T, bg="steelblue2", fg=NULL))
END
;
      $R->run($str);
  }
}

sub analyzeSurv1 {
  my $pre = "/booleanfs2/sahoo/Data/Colon/CDX2/Analysis/jstom-colon";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $hash = &U::getHash("Supplementary/jstom-cdata-1.txt");
  my $ct = 250;
  my $maxt = 300;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } $start .. $end ];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } $start .. $end ];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high:", "green", $g1],
    ["MACS low:", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/jstom-surv-1.png";
  my $res = $h->plotSurvival($outfile, $pG, 120, 130, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Disease-free Survival');

  my $stage = $h->getSurvivalName('c Stage');
  my $stage = &Hegemon::convertValues($stage,
      'I' => 1, 'II'=> 2, 'III' => 3, 'IV' => 4);
  my $grade = $h->getSurvivalName('c Grade');
  my $age = $h->getSurvivalName('n Age');
  my $sex1 = $h->getSurvivalName('c Sex');
  my $sex = &Hegemon::convertValues($sex1, 'M' => 1, 'F' => 3, "" => 2);
  my $res1 = $h->univariate(60, 'groups' => $groups, 'stage' => $stage,
      'age' => $age, 'sex' => $sex);
  my $res2 = $h->multivariate(60, 'groups' => $groups, 'stage' => $stage,
      'age' => $age, 'sex' => $sex);
  &Hegemon::printUMStatsAll($res1, $res2, 'groups', 'stage',
      'age', 'sex');
  my $ofile = "Supplementary/jstom-dfs-table.tex";
  &Hegemon::printUMStatsTexAll($ofile, $res1, $res2, 'groups', 'stage',
      'age', 'sex');

  print "With Grade\n";
  my $res1 = $h->univariate(60, 'groups' => $groups, 'stage' => $stage,
      'grade' => $grade, 'age' => $age, 'sex' => $sex);
  my $res2 = $h->multivariate(60, 'groups' => $groups, 'stage' => $stage,
      'grade' => $grade, 'age' => $age, 'sex' => $sex);
  &Hegemon::printUMStatsAll($res1, $res2, 'groups', 'stage',
      'grade', 'age', 'sex');
  my $ofile = "Supplementary/jstom-grade-table.tex";
  &Hegemon::printUMStatsTexAll($ofile, $res1, $res2, 'groups', 'stage',
      'grade', 'age', 'sex');

  print "With Stage\n";
  my $res1 = $h->univariate(60, 'groups' => $groups, 'stage' => $stage);
  my $res2 = $h->multivariate(60, 'groups' => $groups, 'stage' => $stage);
  &Hegemon::printUMStatsAll($res1, $res2, 'groups', 'stage');
  my $ofile = "Supplementary/jstom-stage-table.tex";
  &Hegemon::printUMStatsTexAll($ofile, $res1, $res2, 'groups', 'stage');

  my $pre = "/booleanfs2/sahoo/Data/Colon/Kirzin-Sheffer/kirzin-sheffer-colon";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $tissue = $h->getSurvivalName('c Tissue');
  my $all = [ grep { $tissue->[$_] eq "Primary Tumor" } $start .. $end ];
  my $hash = &U::getHash("Supplementary/sheffer-cdata-1.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/sheffer-surv-1.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Disease-free Survival');

  my $pre = "/booleanfs2/sahoo/Data/Colon/Gaedcke/gaedcke-2014-colon";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $all = [ $start .. $end ];
  my $hash = &U::getHash("Supplementary/gaedcke-cdata-1.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/gaedcke-surv-1.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Disease-free Survival');

  my $pre = "/booleanfs2/sahoo/Data/Piero/Colon/tcga-2017-m";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $tissue = $h->getSurvivalName('c Histology');
  my $all = [ grep { $tissue->[$_] eq "Primary Tumor" } $start .. $end ];
  my $hash = &U::getHash("Supplementary/tcga-cdata-1.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/tcga-surv-1.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Overall Survival');

  my $pre = "/booleanfs2/sahoo/Data/Colon/Cancer/del-rio-2017-crc";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $all = [ $start .. $end ];
  my $hash = &U::getHash("Supplementary/del-rio-1-cdata-1.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/del-rio-1-surv-1.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Overall Survival');

  my $pre = "/booleanfs2/sahoo/Data/Colon/Cancer/hu-2018-crc";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $tissue = $h->getSurvivalName('c tissue');
  my $all = [ grep { $tissue->[$_] eq "rectal tumor" } $start .. $end ];
  my $hash = &U::getHash("Supplementary/hu-cdata-1.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/hu-surv-1.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Disease-free Survival');

  my $pre = "/booleanfs2/sahoo/Data/Colon/Cancer/allen-2018-crc";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $all = [ $start .. $end ];
  my $hash = &U::getHash("Supplementary/allen-cdata-1.txt");
  my $ct = 250;
  my $maxt = 300;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/allen-surv-1.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Overall Survival');

  my $pre = "/booleanfs/sahoo/Networks/Survival/khambata-ford-2007-colon";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $all = [ $start .. $end ];
  my $hash = &U::getHash("Supplementary/kf-cdata-1.txt");
  my $ct = 400;
  my $maxt = 450;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 80 } 0 .. 5];
  my $outfile = "Supplementary/kf-surv-1.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Days', 'ylabel' => 'Progression-free Survival');

}

sub analyzeSurv2 {
  my $pre = "/booleanfs2/sahoo/Data/Colon/CDX2/Analysis/jstom-colon";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $hash = &U::getHash("Supplementary/jstom-cdata-2.txt");
  my $ct = 250;
  my $maxt = 300;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } $start .. $end ];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } $start .. $end ];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high:", "green", $g1],
    ["MACS low:", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/jstom-surv-2.png";
  my $res = $h->plotSurvival($outfile, $pG, 120, 130, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Disease-free Survival');

  my $stage = $h->getSurvivalName('c Stage');
  my $stage = &Hegemon::convertValues($stage,
      'I' => 1, 'II'=> 2, 'III' => 3, 'IV' => 4);
  my $grade = $h->getSurvivalName('c Grade');
  my $age = $h->getSurvivalName('n Age');
  my $sex1 = $h->getSurvivalName('c Sex');
  my $sex = &Hegemon::convertValues($sex1, 'M' => 1, 'F' => 3, "" => 2);
  my $res1 = $h->univariate(60, 'groups' => $groups, 'stage' => $stage,
      'age' => $age, 'sex' => $sex);
  my $res2 = $h->multivariate(60, 'groups' => $groups, 'stage' => $stage,
      'age' => $age, 'sex' => $sex);
  &Hegemon::printUMStatsAll($res1, $res2, 'groups', 'stage',
      'age', 'sex');
  my $ofile = "Supplementary/jstom-dfs-table.tex";
  &Hegemon::printUMStatsTexAll($ofile, $res1, $res2, 'groups', 'stage',
      'age', 'sex');

  print "With Grade\n";
  my $res1 = $h->univariate(60, 'groups' => $groups, 'stage' => $stage,
      'grade' => $grade, 'age' => $age, 'sex' => $sex);
  my $res2 = $h->multivariate(60, 'groups' => $groups, 'stage' => $stage,
      'grade' => $grade, 'age' => $age, 'sex' => $sex);
  &Hegemon::printUMStatsAll($res1, $res2, 'groups', 'stage',
      'grade', 'age', 'sex');
  my $ofile = "Supplementary/jstom-grade-table.tex";
  &Hegemon::printUMStatsTexAll($ofile, $res1, $res2, 'groups', 'stage',
      'grade', 'age', 'sex');

  print "With Stage\n";
  my $res1 = $h->univariate(60, 'groups' => $groups, 'stage' => $stage);
  my $res2 = $h->multivariate(60, 'groups' => $groups, 'stage' => $stage);
  &Hegemon::printUMStatsAll($res1, $res2, 'groups', 'stage');
  my $ofile = "Supplementary/jstom-stage-table.tex";
  &Hegemon::printUMStatsTexAll($ofile, $res1, $res2, 'groups', 'stage');

  my $pre = "/booleanfs2/sahoo/Data/Colon/Kirzin-Sheffer/kirzin-sheffer-colon";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $tissue = $h->getSurvivalName('c Tissue');
  my $all = [ grep { $tissue->[$_] eq "Primary Tumor" } $start .. $end ];
  my $hash = &U::getHash("Supplementary/sheffer-cdata-2.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/5;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] + $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] + $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/sheffer-surv-2.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Disease-free Survival');

  return;
  my $pre = "/booleanfs2/sahoo/Data/Colon/Gaedcke/gaedcke-2014-colon";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $all = [ $start .. $end ];
  my $hash = &U::getHash("Supplementary/gaedcke-cdata-2.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/gaedcke-surv-2.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Disease-free Survival');

  my $pre = "/booleanfs2/sahoo/Data/Piero/Colon/tcga-2017-m";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $tissue = $h->getSurvivalName('c Histology');
  my $all = [ grep { $tissue->[$_] eq "Primary Tumor" } $start .. $end ];
  my $hash = &U::getHash("Supplementary/tcga-cdata-2.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/tcga-surv-2.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Overall Survival');

  my $pre = "/booleanfs2/sahoo/Data/Colon/Cancer/del-rio-2017-crc";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $all = [ $start .. $end ];
  my $hash = &U::getHash("Supplementary/del-rio-1-cdata-2.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/del-rio-1-surv-2.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Overall Survival');

  my $pre = "/booleanfs2/sahoo/Data/Colon/Cancer/hu-2018-crc";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $tissue = $h->getSurvivalName('c tissue');
  my $all = [ grep { $tissue->[$_] eq "rectal tumor" } $start .. $end ];
  my $hash = &U::getHash("Supplementary/hu-cdata-2.txt");
  my $ct = 150;
  my $maxt = 160;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] + $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] + $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/hu-surv-2.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Disease-free Survival');

  my $pre = "/booleanfs2/sahoo/Data/Colon/Cancer/allen-2018-crc";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $all = [ $start .. $end ];
  my $hash = &U::getHash("Supplementary/allen-cdata-2.txt");
  my $ct = 250;
  my $maxt = 300;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 24 } 0 .. 5];
  my $outfile = "Supplementary/allen-surv-2.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Months', 'ylabel' => 'Overall Survival');

  my $pre = "/booleanfs/sahoo/Networks/Survival/khambata-ford-2007-colon";
  my $thrf = "$pre\-thr.txt";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => $thrf,
    survival => "$pre\-survival.txt");

  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $all = [ $start .. $end ];
  my $hash = &U::getHash("Supplementary/kf-cdata-2.txt");
  my $ct = 400;
  my $maxt = 450;

  my $f_ranks = [ map { $hash->{$h->{'headers'}->[$_]}->[2] } $start .. $end];
  my $data = [sort { $a <=> $b } @$f_ranks];
  my $res = &U::fitStep($data, 0, scalar(@$data) - 1);
  my $thr = [$res->[6], $res->[3], $res->[6]-0.5, $res->[6]+0.5];
  my $nm = (&U::max($f_ranks) - &U::min($f_ranks))/15;
  my $nm = 0;
  my $g1 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] > ($thr->[0] - $nm) } @$all];
  my $g2 = [ grep { $hash->{$h->{'headers'}->[$_]}->[2] <= ($thr->[0] - $nm) } @$all];
  print STDERR join("\t", "DR", @$thr), "\n";

  my $groups = [ map { "" } 0 .. $end];
  $groups->[$_] = 1 foreach (@$g1);
  $groups->[$_] = 2 foreach (@$g2);

  my $pG = [ ["MACS high :", "green", $g1],
    ["MACS low  :", "red",   $g2] ];
  foreach my $g (@$pG) {
    print join("\t", $g->[0], scalar(@{$g->[2]})), "\n";
  }
  my $times = [map { $_ * 80 } 0 .. 5];
  my $outfile = "Supplementary/kf-surv-2.png";
  my $res = $h->plotSurvival($outfile, $pG, $ct, $maxt, undef, $times);
  &Hegemon::printSurvStats($res, $times);
  $outfile =~ s/png/tex/g;
  &Hegemon::printSurvTikz($outfile, $res, $times, 
      'xlabel' => 'Days', 'ylabel' => 'Progression-free Survival');

}

