#!/usr/bin/perl -I ..
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
if ($cmd eq "network") {
  &analyzeNetwork(@ARGV);
}
if ($cmd eq "n-cls") {
  &analyzeNameCLS(@ARGV);
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
    my $file = "nb-net-res.txt";
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
    my $file = "nb-net-eq.txt";
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
    my $file = "nb-net.rl";
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
    my $file = "nb-net-eq.txt";
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
    my $file = "nb-net-cls.txt";
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
    my $file = "nb-net.rl";
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
    my $file = "nb-net-eq-g.txt";
    my $nodelist = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      $g->addEdge($list[0], $list[2], $list[1]);
    }
    close($fh);
    my $file = "nb-net-cls.txt";
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
    my $pre = "zhang-2014-rnaseq-nb";
    my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => "$pre\-survival.txt");
    my $ehash = {};
    my $ehash1 = {};
    my ($nodes, $edges) = $g->getNodesEdges("5");
    if (defined $params[1]) {
      my $ids = [];
      foreach my $id (@params) {
        my $l = $h->getIDs($id);
        next if ! defined $l;
        push @$ids, @$l;
      }
      foreach my $id (@$ids) {
        next if (!defined $nhash->{$id});
        my $i = $g->getIndex($nhash->{$id});
        if (defined $ehash1->{$g->getNode($i)}) {
          $rank->{$i} += 100;
        }
        $ehash1->{$g->getNode($i)}->{$id} = 1;
        $ehash->{$g->getNode($i)}->{$id} = 1;
        if ($rank->{$i} < 2500) {
          $rank->{$i} += 4000;
        }
        foreach my $k (keys %{$edges->{$i}}) {
          if ($rank->{$k} > 5 && $rank->{$k} < 2500) {
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

    foreach my $k (@$res) {
      foreach my $g (@{$k->[2]}) {
        $g->[1] .= " ".$h->getSimpleName($g->[0]);
      }
      foreach my $i (0 .. $#{$k->[3]}) {
        $k->[3]->[$i] .= " (".$h->getSimpleName($k->[3]->[$i]).")";
      }
    }

    my $json = encode_json $res;
    print $json, "\n";
    #foreach (@$res) {
    #  print join("\t", $_->[0], $_->[1], map { $_->[0]." ".$_->[1]} @{$_->[2]}), "\n";
    #}
  }
  if (defined $params[0] && $params[0] eq "query") {
    my $file = "nb-net-cls.txt";
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
    my $file = "nb-net.rl";
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
    my $file = "nb-net-eq-g.txt";
    my $nodelist = {};
    my $ids = [];
    open(my $fh, "<$file");
    while (<$fh>) {
      chomp;
      my @list = split("\t");
      $g->addEdge($list[0], $list[2], $list[1]);
    }
    close($fh);
    my $file = "nb-net-cls.txt";
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
    my $pre = "zhang-2014-rnaseq-nb";
    my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => "$pre\-survival.txt");
    my $u = $nhash->{$params[1]};
    my $n1 = scalar(@{$nodelist->{$u}});
    print "$u ". $h->getSimpleName($u). " ($n1)\n";
    my $l = $nodelist->{$u};
    my $v = $nhash->{$params[2]};
    my $n2 = scalar(@{$nodelist->{$v}});
    print "$v ". $h->getSimpleName($v). " ($n2)\n";
    my $res = $g->Path("All", $u, $v, 0, 100);
    print join("\t", map { $g->getNode($res->[$_])." ".
        $h->getSimpleName($g->getNode($res->[$_])). " ".
     $g->getType($res->[$_], $res->[$_ + 1]) } 0 .. $#{$res}), "\n";
  }
  if (defined $params[0] && $params[0] eq "state") {
    my $file = "nb-net-cls.txt";
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
    my $pre = "zhang-2014-rnaseq-nb";
    my $h = Hegemon->new(expr => "$pre\-expr.txt",
      idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
      survival => "$pre\-survival.txt");
    my $start = $h->getStart();
    my $end = $h->getEnd();
    my $mtype = $h->getSurvivalName("c mycn status");
    my $m0 = [ grep { $mtype->[$_] eq "0" } $start .. $end ];
    my $m1 = [ grep { $mtype->[$_] eq "1" } $start .. $end ];

    my $pG = [ ["All", "red", [ $start .. $end ] ] ];
    my $pG = [
      ["Base", "#0010EE", $m0],
      ["Amp", "#000000", $m1]];
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

sub analyzeNameCLS {
  my $pre = "zhang-2014-rnaseq-nb";
  my $h = Hegemon->new(expr => "$pre\-expr.txt",
    idx => "$pre\-idx.txt", thr => "$pre\-thr.txt",
    survival => "$pre\-survival.txt");
  my $file = "nb-net-cls.txt";
  open(my $fh, "<$file");
  while (<$fh>) {
    chomp;
    my @list = split("\t");
    my $id = $list[0];
    my $l = [@list[2 .. $#list]];
    my $n = $h->getSimpleName($id);
    print join("\t", "$id ($n)", $list[1], 
      map { $h->getSimpleName($_) } @$l), "\n";
  }
  close($fh);
}

