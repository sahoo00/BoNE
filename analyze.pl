#!/usr/bin/perl -I .
# Author: Debashis Sahoo <dsahoo@ucsd.edu>

if (scalar(@ARGV) <= 0) {
    print "perl analyze.pl <cmd> <args>\n";
    exit(1);
}

use U;
use Network;
use Graph;
use Hegemon;
use LWP::UserAgent;
use Data::Dumper;
use JSON;

my $cmd = shift(@ARGV);
if ($cmd eq "toidx") {
  &U::convertidx(@ARGV);
}
if ($cmd eq "thr") {
  &U::printThr(@ARGV);
}
if ($cmd eq "Info") {
  &convertInfo(@ARGV);
}
if ($cmd eq "VInfo") {
  &convertVInfo(@ARGV);
}
if ($cmd eq "bv") {
  &convertBv(@ARGV);
}
if ($cmd eq "data-download") {
  &downloadData(@ARGV);
}
if ($cmd eq "ibd") {
  &analyzeIBD(@ARGV);
}

sub convertInfo {
  my ($dataset, $thr, @rest) = @_;
  my $pre = $dataset;
  $thr = "$pre\-thr.txt" if (!defined $thr);
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => $thr,
      survival => undef);
  $h->printInfo();
}

sub convertVInfo {
  my ($dataset, $thr, @rest) = @_;
  my $pre = $dataset;
  $thr = "$pre\-thr.txt" if (!defined $thr);
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => $thr,
      survival => undef);
  $h->printVInfo();
}

sub convertBv {
  my ($dataset, @rest) = @_;
  my $pre = $dataset;
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => undef);
  $h->printBvFile();
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
#https://arxiv.org/pdf/1606.00950.pdf
#https://link.springer.com/chapter/10.1007/978-3-319-91452-7_13
sub DCT {
  my ($nodes, $edges) = @_;
  my $T = [];
  my $S = [];
  my $V = [keys %{$nodes}];
  my $VB = {};
  my $checked = {};
  my $density = {};
  my $connect = {};
  my $count = 0;
  my $num = scalar(@$V);
  if (scalar(@$V) > 0) {
    foreach my $u (@$V) {
      if (scalar(keys %{$edges->{$u}}) == 0) {
        $checked->{$u} = 1;
        push @$S, $u;
        $num --;
      }
      else {
        $VB->{$u} = 1;
      }
    }
    foreach my $u (@$V) {
      if (!defined $checked->{$u} && defined $VB->{$u}) {
        $checked->{$u} = 1;
        push @$T, $u;
        $count ++;
        last;
      }
    }
  }
  my $hash = {};
  while ($count < $num) {
    my $maxv = -1;
    my ($p, $q);
    foreach my $u (@$T) {
      my $gu = &gamma($u, $edges, $hash);
      foreach my $v (@$gu) {
        if (!defined $checked->{$v} && defined $VB->{$v}) {
          my $r = &rho($u, $v, $edges, $hash);
          if ($r > $maxv) {
            $maxv = $r;
            $p = $v;
            $q = $u;
          }
        }
      }
    }
    if (defined $p) {
      $checked->{$p} = 1;
      $connect->{$p} = $q;
      $density->{$p} = $maxv;
      push @$T, $p;
      $count++;
    }
    else {
      foreach my $u (@$V) {
        if (!defined $checked->{$u} && defined $VB->{$u}) {
          $checked->{$u} = 1;
          push @$T, $u;
          $count++;
          $p = $u;
          my $gu = &gamma($u, $edges);
          $maxv = scalar(@$gu);
          last;
        }
      }
    }
    print STDERR "$count/$num $p $q $maxv\n";
  }
  return [$T, $density, $connect];
}

sub DCT1 {
  my ($nodes, $edges) = @_;
  print STDERR scalar(keys %{$nodes}), "\n";
  my $V = [keys %{$nodes}];
  my $VB = {};
  foreach my $u (@$V) {
    if (scalar(keys %{$edges->{$u}}) > 0) {
      $VB->{$u} = 1;
    }
  }
  my ($nodes_2, $edges_2) = &Graph::subGraph1($nodes, $edges, $VB);
  print STDERR scalar(keys %{$nodes_2}), "\n";
  return &DCT2($nodes_2, $edges_2);
}

sub DCT2 {
  my ($nodes, $edges) = @_;
  my $T = [];
  my $V = [keys %{$nodes}];
  my $checked = {};
  my $density = {};
  my $connect = {};
  my $count = 0;
  my $num = scalar(@$V);
  if (scalar(@$V) > 0) {
    foreach my $u (@$V) {
      if (!defined $checked->{$u}) {
        $checked->{$u} = 1;
        push @$T, $u;
        $count ++;
        last;
      }
    }
  }
  my $hash = {};
  my $shash = {};
  while ($count < $num) {
    my $maxv = -1;
    my ($p, $q);
    foreach my $u ($T->[$#{$T}]) {
      my $gu = &gamma($u, $edges);
      foreach my $v (@$gu) {
        if (!defined $checked->{$v}) {
          my $r = &rho($u, $v, $edges);
          if ($r > $maxv) {
            $maxv = $r;
            $p = $v;
            $q = $u;
          }
        }
      }
    }
    if (defined $p) {
      $checked->{$p} = 1;
      $connect->{$p} = $q;
      $density->{$p} = $maxv;
      push @$T, $p;
      $count++;
    }
    else {
      foreach my $u (@$V) {
        if (!defined $checked->{$u}) {
          $checked->{$u} = 1;
          push @$T, $u;
          $count++;
          $p = $u;
          my $gu = &gamma($u, $edges, $hash);
          $maxv = scalar(@$gu);
          last;
        }
      }
    }
    print STDERR "$count/$num $p $q $maxv\n";
  }
  return [$T, $density, $connect];
}

sub getCodes {
  my ($n, $l, $balanced, $balancedhash) = @_;
  my $num = $n->getNum();
  my $relations = {};
  foreach my $u (@$l) {
    my $i = $balancedhash->{$u};
    for (my $j = 0; $j < $num/2; $j++) {
      my $code = $n->readCode_1_1($i, $j);
      if ($code > 0) {
        if (!defined $relations->{$code}) {
          $relations->{$code} = {};
        }
        $relations->{$code}->{$balanced->[$j]} = 1;
      }
    }
  }
  return $relations;
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

sub getNCodes {
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
      if ($code > 0) {
        if (!defined $relations->{$code}) {
          $relations->{$code} = {};
        }
        $relations->{$code}->{$balanced->[$j]}++;
      }
    }
  }
  return $relations;
}

sub getNCodesArijs {
  my ($arijs, $n, $net, $l, $balanced, $balancedhash) = @_;
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
      my $res = &getArijsCode($arijs, $u, $balanced->[$j]);
      if (&isCodeCompatible($code, $res)) {
        if (!defined $relations->{$code}) {
          $relations->{$code} = {};
        }
        $relations->{$code}->{$balanced->[$j]}++;
      }
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

sub getNCodesPair2 {
  my ($n, $net, $u, $v, $balanced, $balancedhash) = @_;
  my $num = $n->getNum();
  my $relations = {};
  my $nb = $n->{'numBits'};
  my $vcode = Bit::Vector->new($nb*2*2);
  my $i = $balancedhash->{$u};
  my $j = $balancedhash->{$v};
  $n->readNCode($net, $i, $j, $vcode);
  my $code = &convertCodes($vcode);
  return $code;
}

sub getRelationship {
  my $counts = shift;
  if (scalar(grep { $counts->[$_]->[0] eq 0 } 1 .. 6) == 6) {
    return 0;
  }
  my $ratio = [];
  foreach my $c (0 .. 6) {
    if ($counts->[$c]->[1] > 0) {
      $ratio->[$c] = $counts->[$c]->[0]/$counts->[$c]->[1];
    }
    else {
      $ratio->[$c] = 0;
    }
  }
  my $m1 = &U::max([ map { $counts->[$_]->[0] } (2, 3, 5)]);
  my $m2 = &U::max([ map { $counts->[$_]->[0] } (1, 4, 6)]);
  if ($m1 == $m2) {
    return 0;
  }
  my $m1 = &U::max([ map { $ratio->[$_] } (2, 3, 5)]);
  my $m2 = &U::max([ map { $ratio->[$_] } (1, 4, 6)]);
  if ($m1 == $m2) {
    return 0;
  }
  if ($m1 > $m2) {
    if ($m1 == $ratio->[5] && $m1 > 0.5) {
      return 5;
    }
    if ($m1 == $ratio->[3] && $m1 > 0.5) {
      return 3;
    }
    if ($m1 == $ratio->[2] && $m1 > 0.5) {
      return 2;
    }
  }
  else {
    if ($m2 == $ratio->[6] && $m2 > 0.5) {
      return 6;
    }
    if ($m2 == $ratio->[4] && $m2 > 0.5) {
      return 4;
    }
    if ($m2 == $ratio->[1] && $m2 > 0.5) {
      return 1;
    }
  }
  return 0;
}

sub getRelationship2 {
  my ($counts, $i) = shift;
  if (scalar(grep { $counts->[$_]->[$i] eq 0 } 1 .. 6) == 6) {
    return 0;
  }
  my $ratio = [];
  foreach my $c (0 .. 6) {
    if ($counts->[$c]->[$i + 1] > 0) {
      $ratio->[$c] = $counts->[$c]->[$i]/$counts->[$c]->[$i + 1];
    }
    else {
      $ratio->[$c] = 0;
    }
  }
  my $m1 = &U::max([ map { $counts->[$_]->[$i] } (2, 3, 5)]);
  my $m2 = &U::max([ map { $counts->[$_]->[$i] } (1, 4, 6)]);
  if ($m1 == $m2) {
    return 0;
  }
  my $m1 = &U::max([ map { $ratio->[$_] } (2, 3, 5)]);
  my $m2 = &U::max([ map { $ratio->[$_] } (1, 4, 6)]);
  if ($m1 == $m2) {
    return 0;
  }
  if ($m1 > $m2) {
    if ($m1 == $ratio->[5] && $m1 > 0.5) {
      return 5;
    }
    if ($m1 == $ratio->[3] && $m1 > 0.5) {
      return 3;
    }
    if ($m1 == $ratio->[2] && $m1 > 0.5) {
      return 2;
    }
  }
  else {
    if ($m2 == $ratio->[6] && $m2 > 0.5) {
      return 6;
    }
    if ($m2 == $ratio->[4] && $m2 > 0.5) {
      return 4;
    }
    if ($m2 == $ratio->[1] && $m2 > 0.5) {
      return 1;
    }
  }
  return 0;
}

sub getRelationship3 {
  my $counts = shift;
  my $c1 = &getRelationship2($counts, 0);
  my $c2 = &getRelationship2($counts, 2);
  if ($c1 == $c2) {
    return $c1;
  }
  my $l = &U::intersection([$c1, $c2], [2, 3, 5]);
  if (scalar(@$l) == 2) { return 5; }
  my $l = &U::intersection([$c1, $c2], [1, 4, 6]);
  if (scalar(@$l) == 2) { return 6; }
  return 0;
}

sub getRel2 {
  my ($n, $net, $lu, $lv, $r1, $r2, $balanced, $balancedhash) = @_;
  my $nb = $n->{'numBits'};
  my $vcode = Bit::Vector->new($nb*2*2);
  foreach my $lid (@$lu) {
    my $li = $balancedhash->{$lid};
    foreach my $ljd (@$lv) {
      my $lj = $balancedhash->{$ljd};
      $n->readNCode($net, $li, $lj, $vcode);
      my $code = &convertCodes($vcode);
      if ($code > 0) {
        if (!defined $r1->{$code}) {
          $r1->{$code} = {};
        }
        if (!defined $r2->{$code}) {
          $r2->{$code} = {};
        }
        $r1->{$code}->{$balanced->[$li]} = 1;
        $r2->{$code}->{$balanced->[$lj]} = 1;
      }
    }
  }
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

sub getArijsNet {
  my $pre = "data/arijs-2018-uc";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
    survival => "$pre\-survival.txt");
  my $start = $h->getStart();
  my $end = $h->getEnd();
  my $all = [ $start .. $end ];
  use Network;
  my $file = "results/ibd-network-2.rl";
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
  my $hash = &U::getHash("data/peters-arjis-mapping.txt");
  my $res = { h => $h, n => $n, net => $net, hash => $hash,
    all => $all, balanced => $balanced, balancedhash => $balancedhash};
  return $res;
}

sub getArijsCode {
  my ($obj, $id1, $id2) = @_;
  my $all = $obj->{'all'};
  my $hash = $obj->{'hash'};
  my $h = $obj->{'h'};
  my $n = $obj->{'n'};
  my $net = $obj->{'net'};
  my $balanced = $obj->{'balanced'};
  my $balancedhash = $obj->{'balancedhash'};
  return [0, 0] if (!defined $hash->{$id1});
  return [0, 0] if (!defined $hash->{$id2});
  my $u = $hash->{$id1}->[1];
  my $v = $hash->{$id2}->[1];
  my $ex_u = $h->getExprData($u);
  my @exp_u = map { $ex_u->[$_] } @$all;
  my $ex_v = $h->getExprData($v);
  my @exp_v = map { $ex_v->[$_] } @$all;
  my $corr = &U::correlation(\@exp_u, \@exp_v);
  my $code = &getNCodesPair2($n, $net, $u, $v, $balanced, $balancedhash);
  return [$code, $corr]
}

sub isCodeCompatible {
  my ($code, $res) =  @_;
  if ($code > 0) {
    return 1 if ($code eq $res->[0]);
    if ($res->[1] > 0.5) {
      return 1 if ($code eq 2 || $code eq 3 || $code eq 5);
    }
    if ($res->[1] < -0.25) {
      return 1 if ($code eq 1 || $code eq 4 || $code eq 6);
    }
  }
  return 0;
}

sub downloadData {
  my ($dbid, $pre) = @_;
  my $browser = LWP::UserAgent->new;
  my $url = "http://hegemon.ucsd.edu/Tools/explore.php";
  my $form = {'go' => 'dataDownload', 'id' => $dbid, 'genes' => '', 'groups' => '',
    'key' => 'polyps', 'param' => 'type:expr'};
  my $res = $browser->post($url, $form);
  my $ef = "$pre\-expr.txt";
  open(my $fh, ">$ef");
  print $fh $res->content;
  close($fh);
  my $form = {'go' => 'dataDownload', 'id' => $dbid, 'genes' => '', 'groups' => '',
    'key' => 'polyps', 'param' => 'type:indexHeader'};
  my $res = $browser->post($url, $form);
  my $ef = "$pre\-ih.txt";
  open(my $fh, ">$ef");
  print $fh $res->content;
  close($fh);
  $form->{'param'} = 'type:survival';
  my $res = $browser->post($url, $form);
  my $ef = "$pre\-survival.txt";
  open(my $fh, ">$ef");
  print $fh $res->content;
  close($fh);
}

sub analyzeIBD {
  my @params = @_;
  if (defined $params[0] && $params[0] eq "eq") {
    use Graph;
    my $g = Graph->new();
    my $file = "results/ibd-network-res.txt";
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
  if (defined $params[0] && $params[0] eq "eq-corr") {
    my $pre = "data/arijs-2018-uc";
    my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => "$pre\-survival.txt");
    my $start = $h->getStart();
    my $end = $h->getEnd();
    my $all = [ $start .. $end ];
    my $hash = &U::getHash("data/peters-arjis-mapping.txt");
    my $file = "data/ibd-network-g-eq-4.txt";
    my $index = 0;
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id1 = $list[0];
      my $id2 = $list[2];
      next if (!defined $hash->{$id1});
      next if (!defined $hash->{$id2});
      my $u = $hash->{$id1}->[1];
      my $v = $hash->{$id2}->[1];
      my $ex_u = $h->getExprData($u);
      my @exp_u = map { $ex_u->[$_] } @$all;
      my $ex_v = $h->getExprData($v);
      my @exp_v = map { $ex_v->[$_] } @$all;
      my $corr = &U::correlation(\@exp_u, \@exp_v);
      print join("\t", $id1, $id2, $list[1], $corr), "\n";
      $index++;
    }
    close($fh);
  }
  if (defined $params[0] && $params[0] eq "eq-corr-thr") {
    my $file = "results/ibd-network-g-eq-4-corr.txt";
    my $index = 0;
    open(my $fh, "<$file");
    my $data = [];
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      if ($list[3] > 0) {
        push @$data, $list[3];
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
    my $hash = &U::getHash("data/peters-arjis-mapping.txt");
    use Network;
    my $file = "results/ibd-network-2.rl";
    my $n = Network->new(file => $file, mode => "r");
    $n->readNetworkFile();
    my $balanced = $n->getBalanced();
    my $balancedhash = {};
    for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
      $balancedhash->{$balanced->[$i]} = $i;
      my $id = $balanced->[$i];
      if (defined $hash->{$id}) {
        $uf->makeset($id);
        $nodes->{$id} = {} if (!defined $nodes->{$id});
        $nodes->{$id}->{$id} = 1;
      }
    }
    my $num = $n->getNum();
    my $rhash = {};
    my $net = $n->readAll();
    my $file = "data/ibd-network-g-eq-4.txt";
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id1 = $list[0];
      my $id2 = $list[2];
      if (defined $hash->{$id1}) {
        $uf->makeset($id1);
        $nodes->{$id1} = {} if (!defined $nodes->{$id1});
        $nodes->{$id1}->{$id1} = 1;
      }
      if (defined $hash->{$id2}) {
        $uf->makeset($id2);
        $nodes->{$id2} = {} if (!defined $nodes->{$id2});
        $nodes->{$id2}->{$id2} = 1;
      }
    }
    close($fh);
    my $file = "results/ibd-network-g-eq-4-corr.txt";
    my $index = 0;
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id1 = $list[0];
      my $id2 = $list[1];
      next if (!defined $hash->{$id1});
      next if (!defined $hash->{$id2});
      $uf->makeset($id1);
      $uf->makeset($id2);
      $nodes->{$id1} = {} if (!defined $nodes->{$id1});
      $nodes->{$id2} = {} if (!defined $nodes->{$id2});
      $nodes->{$id1}->{$id1} = 1;
      $nodes->{$id2}->{$id2} = 1;
      my $u = $hash->{$id1}->[1];
      my $v = $hash->{$id2}->[1];
      if ($list[3] > 0.5) {
        push @$data, $list[3];
        $uf->union($id1, $id2);
        $nodes->{$id1}->{$id2} = $list[1];
        $nodes->{$id2}->{$id1} = $list[1];
      }
      elsif (defined $balancedhash->{$u} && defined $balancedhash->{$v}) {
        my $code = &getNCodesPair2($n, $net, $u, $v, $balanced, $balancedhash);
        next if ($code ne 5 && $code ne 3 && $code ne 2);
        push @$data, $list[3];
        $uf->union($id1, $id2);
        $nodes->{$id1}->{$id2} = $list[3];
        $nodes->{$id2}->{$id1} = $list[3];
      }
      $index++;
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
    my $file = "results/ibd-network-g-eq-cls-4.txt";
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
    my $file = "results/ibd-network-1.rl";
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
    my $arijs = &getArijsNet();
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
      my $ru = &getNCodesArijs($arijs, $n, $net, $l, $balanced, $balancedhash);
      if (defined $params[2]) {
        my $v = $params[2];
        print STDERR $v, "\n";
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
  if (defined $params[0] && $params[0] eq "g-filter") {
    my $arijs = &getArijsNet();
    my $file = "results/ibd-network-g-eq-g-4.txt";
    my $index = 0;
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      my $id1 = $list[0];
      my $id2 = $list[2];
      my $res = &getArijsCode($arijs, $id1, $id2);
      if (&isCodeCompatible($list[1], $res)) {
        print join("\t", @list, $res->[1]), "\n";
      }
      $index++;
      if ( ($index % 10000) == 0 ) {
        print STDERR "$index\n";
      }
    }
    close($fh);
  }
  if (defined $params[0] && $params[0] eq "genes") {
    use Graph;
    my $g = Graph->new();
    my $file = "data/ibd-network-g-eq-g-4.txt";
    my $nodelist = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      $g->addEdge($list[0], $list[2], $list[1]);
    }
    close($fh);
    my $file = "data/ibd-network-g-eq-cls-4.txt";
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
    my ($nodes, $edges) = $g->getNodesEdges("5");
    if (defined $params[1]) {
      foreach my $id (@params) {
        next if (!defined $nhash->{$id});
        my $i = $g->getIndex($nhash->{$id});
        $ehash->{$g->getNode($i)}->{$id} = 1;
        if ($rank->{$i} < 2000) {
          $rank->{$i} += 5000;
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
    my $file = "results/ibd-network-g-eq-cls-4.txt";
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
    my $file = "results/ibd-network-1.rl";
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
    my @rels = (1 .. 6);
    if (defined $params[3]) {
      @rels = grep { $_ =~ /[0-9]/ } split("", $params[3]);
    }
    foreach my $c (@rels) {
      foreach my $v (keys %{$ru->{$c}}) {
        next if (!defined $vhash->{$v});
        print join("\t", $u, $c, $v, $ru->{$c}->{$v}), "\n";
      }
    }
  }
  if (defined $params[0] && $params[0] eq "path") {
    use Graph;
    my $g = Graph->new();
    my $file = "results/ibd-network-g-eq-g-4.txt";
    my $nodelist = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      $g->addEdge($list[0], $list[2], $list[1]);
    }
    close($fh);
    my $file = "results/ibd-network-g-eq-cls-4.txt";
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
    my $res = $g->Path("5", $u, $v, 0, 100);
    print join("\t", map { $g->getNode($res->[$_])." ".
     $g->getType($res->[$_], $res->[$_ + 1]) } 0 .. $#{$res}), "\n";
  }
  if (defined $params[0] && $params[0] eq "rank") {
    my $file = "results/ibd-network-g-eq-cls-4.txt";
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
    my $u = $nhash->{$params[1]};
    my $n1 = scalar(@{$nodelist->{$u}});
    print "$u ($n1)\n";
    my $l = $nodelist->{$u};
    my $v = $nhash->{$params[2]};
    my $n2 = scalar(@{$nodelist->{$v}});
    print "$v ($n2)\n";
    my @rels = (1 .. 6);
    if (defined $params[3]) {
      @rels = grep { $_ =~ /[0-9]/ } split("", $params[3]);
    }
    my $vhash = {};
    if (defined $params[4]) {
      $vhash = &U::getHash($params[4], 1);
    }
    else {
      foreach my $v1 (@{$nodelist->{$v}}) {
        $vhash->{$v1} = 1;
      }
    }
    use Network;
    my $file = "results/ibd-network-1.rl";
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
    my $ru = &getNCodes($n, $net, $l, $balanced, $balancedhash);
    foreach my $c (@rels) {
      foreach my $v (keys %{$ru->{$c}}) {
        next if (!defined $vhash->{$v});
        print join("\t", $u, $c, $v, $ru->{$c}->{$v}), "\n";
      }
    }
  }
}

