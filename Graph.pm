# Author: Debashis Sahoo <dsahoo@ucsd.edu>
package Graph;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my %params = @_;
  my $self  = {};
  $self->{'file'}  = undef;
  $self->{'nodes'}  = undef;
  $self->{'edges'}  = undef;
  $self->{'index'}  = 0;

  foreach (keys(%params)) {
    $self->{$_} = $params{$_};
  }

  bless ($self, $class);
  return $self;
}

sub init {
  my $self = shift;
  if (defined $self->{'file'}) {
    $self->createGraph($self->{'file'});
  }
  else {
    $self->clear();
  }
}

sub clear {
  my $self = shift;
  $self->{'nodes'}  = {};
  $self->{'nodeID'} = {};
  $self->{'nodeName'} = {};
  $self->{'edges'}  = {};
  $self->{'index'}  = 0;
}

sub createGraph {
  my ($self, $file) = @_;
  if (!defined $file) {
    $self->clear();
    return;
  }
  if ($file =~ /.sif$/i) {
    $self->createGraphSIF($file);
  }
}

sub addNode {
  my ($self, $node, $name) = @_;
  my $index = $self->{'index'};
  if (!defined $self->{'nodeID'}->{$node}) {
    $self->{'nodes'}->{$index} = $node;
    $self->{'nodeID'}->{$node} = $index;
    $self->{'nodeName'}->{$index} = $name;
    $self->{'index'} = $index + 1;
  }
}

sub updateNode {
  my ($self, $node, $name) = @_;
  $self->addNode($node, $name);
  my $index = $self->getIndex($node);
  $self->{'nodeName'}->{$index} = $name;
}

sub deleteNode {
  my ($self, $node) = @_;
  my $index = $self->getIndex($node);
  if (defined $index) {
    delete $self->{'edges'}->{$index};
    foreach my $type (keys %{$self->{'etype'}}) {
      delete $self->{'etype'}->{$type}->{$index};
    }
    foreach my $type (keys %{$self->{'etype'}}) {
      foreach my $i (keys %{$self->{'etype'}->{$type}}) {
        delete $self->{'etype'}->{$type}->{$i}->{$index};
      }
    }
    foreach my $i (keys %{$self->{'edges'}}) {
      delete $self->{'edges'}->{$i}->{$index};
    }
    delete $self->{'nodes'}->{$index};
    delete $self->{'nodeID'}->{$node};
    delete $self->{'nodeName'}->{$index};
  }
}

sub isNode {
  my ($self, $node) = @_;
  return defined $self->{'nodeID'}->{$node};
}
sub getNode {
  my ($self, $index) = @_;
  return $self->{'nodes'}->{$index};
}
sub getName {
  my ($self, $index) = @_;
  return $self->{'nodeName'}->{$index};
}
sub getIndex {
  my ($self, $node) = @_;
  return $self->{'nodeID'}->{$node};
}
sub getType {
  my ($self, $i, $j) = @_;
  return $self->{'edges'}->{$i}->{$j};
}

sub addEdge {
  my ($self, $node1, $node2, $type) = @_;
  if (!defined $type) {
    $type = "e";
  }
  if (!$self->isNode($node1)) {
    $self->addNode($node1);
  }
  if (!$self->isNode($node2)) {
    $self->addNode($node2);
  }
  my $index1 = $self->getIndex($node1);
  my $index2 = $self->getIndex($node2);
  if (!defined $self->{'edges'}->{$index1}) {
    $self->{'edges'}->{$index1} = {};
  }
  $self->{'edges'}->{$index1}->{$index2} = $type;
  if (!defined $self->{'etype'}->{$type}) {
    $self->{'etype'}->{$type} = {};
  }
  if (!defined $self->{'etype'}->{$type}->{$index1}) {
    $self->{'etype'}->{$type}->{$index1} = {};
  }
  $self->{'etype'}->{$type}->{$index1}->{$index2} = 1;
}

sub deleteEdge {
  my ($self, $node1, $node2, $type) = @_;
  #print "deleting $node1 $node2 $type\n";
  my $index1 = $self->getIndex($node1);
  my $index2 = $self->getIndex($node2);
  delete $self->{'edges'}->{$index1}->{$index2};
  if (defined $type) {
    if (defined $self->{'etype'}->{$type}->{$index1}) {
      if (defined $self->{'etype'}->{$type}->{$index1}->{$index2}) {
        delete $self->{'etype'}->{$type}->{$index1}->{$index2};
      }
    }
  }
  else {
    foreach my $type (keys %{$self->{'etype'}}) {
      if (defined $self->{'etype'}->{$type}->{$index1}) {
        if (defined $self->{'etype'}->{$type}->{$index1}->{$index2}) {
          delete $self->{'etype'}->{$type}->{$index1}->{$index2};
        }
      }
    }
  }
}

sub deleteNodesWithEdges {
  my ($self, $type, $hash1, $hash2) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  foreach my $i (keys(%{$edges})) {
    my $ni = $self->getNode($i);
    foreach my $j (keys(%{$edges->{$i}})) {
      my $nj = $self->getNode($j);
      if (defined $hash1->{$ni} && defined $hash2->{$nj}) {
        $self->deleteNode($ni);
        $self->deleteNode($nj);
      }
      if (defined $hash2->{$ni} && defined $hash1->{$nj}) {
        $self->deleteNode($ni);
        $self->deleteNode($nj);
      }
    }
  }
}

sub createGraphSIF {
  my ($self, $file) = @_;
  if (!defined $file) {
    $self->clear();
    return;
  }
  my $fh;
  open($fh, "<$file") || die "Can't open $file\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split(/\s+/);
    if (scalar(@list) <= 0) {
      next;
    }
    $self->addNode($list[0]);
    if (scalar(@list) >= 3) {
      foreach (@list[2 .. $#list]) {
        $self->addNode($_);
        $self->addEdge($list[0], $_, $list[1]);
      }
    }
  }
  close($fh);
}

sub write {
  my ($self, $file) = @_;
  if (!defined $file) {
    die "Undefined file for writing\n";
  }
  if ($file =~ /.sif$/i) {
    $self->writeSIF($file);
  }
  if ($file =~ /.adj$/i) {
    $self->writeSADJ($file);
  }
}

sub writeSIFHandle {
  my ($self, $fh, $edges) = @_;
  foreach my $i (keys %{$edges}) {
    my $ni = $self->getNode($i);
    foreach my $j (keys %{$edges->{$i}}) {
      my $nj = $self->getNode($j);
      my $type = $edges->{$i}->{$j};
      print $fh "$ni $type $nj\n";
    }
  }
}

sub writeSADJHandle {
  my ($self, $fh, $edges) = @_;
  print $fh "0 0\n";
  foreach my $i (keys %{$edges}) {
    my $ni = $self->getNode($i);
    print $fh "$ni " . scalar(keys %{$edges->{$i}});
    foreach my $j (keys %{$edges->{$i}}) {
      my $nj = $self->getNode($j);
      my $type = $edges->{$i}->{$j};
      print $fh " $nj";
    }
    print $fh "\n";
  }
}

sub writeSIF {
  my ($self, $file) = @_;
  if (!defined $file) {
    die "Undefined file for writing\n";
  }
  my $fh;
  open($fh, ">$file") || die "Can't open $file\n";
  $self->writeSIFHandle($fh, $self->{'edges'});
  close($fh);
}

sub writeSADJ {
  my ($self, $file) = @_;
  if (!defined $file) {
    die "Undefined file for writing\n";
  }
  my $fh;
  open($fh, ">$file") || die "Can't open $file\n";
  $self->writeSADJHandle($fh, $self->{'edges'});
  close($fh);
}

sub writeGMLHandle {
  my ($self, $fh, $nodes, $edges, $placed) = @_;
  print $fh "graph [\n";
  print $fh "directed 1\n";
  print $fh "id 1\n";
  foreach my $n (keys(%{$nodes})) {
    my $l = $self->getName($n) . " " . $nodes->{$n};
    my $x = 0;
    my $y = 0;
    my $w = 30.0;
    my $h = 30.0;
    if (defined $placed->{$n}) { 
      $x = $placed->{$n}->[0];
      $y = $placed->{$n}->[1];
    }
    print $fh "node [\n";
    print $fh "\tid $n\n\tlabel \"$l\"\n";
    #print $fh "\tid ".$nodes->{$n}."\n\tlabel \"$l\"\n";
    my $nc = "#3cbbc6";
    if (defined $placed->{'upList'} && defined $placed->{'upList'}->{$n}) {
      $nc = "#3c00c6";
    }
    if (defined $placed->{'upBorder'} && defined $placed->{'upBorder'}->{$n}) {
      $nc = "#ff0000";
    }
    if (defined $placed->{'downBorder'} && defined $placed->{'downBorder'}->{$n}) {
      $nc = "#00ff00";
    }
    print $fh <<END;
\tgraphics [
\t\tx       $x
\t\ty       $y
\t\tw       $w
\t\th       $h
\t\ttype    "ellipse"
\t\twidth   1.00000
\t\tfill    "$nc"
\t\toutline "#000000"
\t]
END
    print $fh "]\n";
  }
  foreach my $n (keys(%{$edges})) {
    foreach my $j (keys(%{$edges->{$n}})) {
      my $nc = "#3b3b3b";
      if (defined $edges->{$j} && defined $edges->{$j}->{$n}) {
        $nc = "#bb0000";
      }
      print $fh "edge [\n";
      print $fh "\tsource $n\n\ttarget $j\n";
      #print $fh "\tsource ".$nodes->{$n}."\n\ttarget ".$nodes->{$j}."\n";
      print $fh <<END;
\tgraphics [
\t\twidth   2
\t\ttype    "line"
\t\tfill    "$nc"
\t\tarrow   "last"
\t]
END
      print $fh "]\n";
    }
  }
  print $fh "]\n";
}

sub writeTIKZHandle {
  my ($self, $fh, $nodes, $edges, $placed) = @_;
  print $fh <<END
\\begin{tikzpicture}[yscale=0.7,xscale=0.7]
\\tikzstyle{VertexStyle}=[shape = circle,%
shading = ball,%
ball color = white,%
very thin,
inner sep=2,%
draw]
\\SetVertexNoLabel
\\SetUpEdge[style={->,>=angle 45,bend right=10},color=gray!40]

END
;
  foreach my $n (keys(%{$nodes})) {
    my $l = $self->getName($n) . " " . $nodes->{$n};
    my $x = 0;
    my $y = 0;
    my $w = 30.0;
    my $h = 30.0;
    if (defined $placed->{$n}) { 
      $x = $placed->{$n}->[0];
      $y = $placed->{$n}->[1];
    }
    $x = $x / 100.0;
    $y = $y / 100.0;
    my $nc = "#3cbbc6";
    if (defined $placed->{'upList'} && defined $placed->{'upList'}->{$n}) {
      $nc = "#3c00c6";
    }
    if (defined $placed->{'upBorder'} && defined $placed->{'upBorder'}->{$n}) {
      $nc = "#ff0000";
    }
    if (defined $placed->{'downBorder'} && defined $placed->{'downBorder'}->{$n}) {
      $nc = "#00ff00";
    }
    $l =~ s/[0-9]*_[a-z_]*at//g;
    $l =~ s/[-,]/ /g;
    print $fh "\\Vertex[x=$x,y=$y,NoLabel=true] {s$n}\n";
    print $fh "\\node[right=0cm of s$n, text width=1cm, text badly centered]{\\tiny $l};\n";

  }
  foreach my $n (keys(%{$edges})) {
    foreach my $j (keys(%{$edges->{$n}})) {
      if (defined $placed && defined $placed->{'crossEdges'}->{$n}->{$j}) {
        print $fh "\\Edge[style=dotted,color=gray!20](s$n)(s$j)\n";
      }
      else {
        print $fh "\\Edge(s$n)(s$j)\n";
      }
    }
  }
  print $fh "\\end{tikzpicture}\n";

}

sub writeGML {
  my ($self, $file, $placed) = @_;
  if (!defined $file) {
    die "Undefined file for writing\n";
  }
  my $fh;
  if ($file ne "-") {
    open($fh, ">$file") || die "Can't open $file\n";
  }
  else {
    $fh = \*STDOUT;
  }
  $self->writeGMLHandle($fh, $self->{'nodes'}, $self->{'edges'}, $placed);
  if ($file ne "-") {
    close($fh);
  }
}

sub writeTIKZ {
  my ($self, $file, $placed) = @_;
  if (!defined $file) {
    die "Undefined file for writing\n";
  }
  my $fh;
  if ($file ne "-") {
    open($fh, ">$file") || die "Can't open $file\n";
  }
  else {
    $fh = \*STDOUT;
  }
  $self->writeTIKZHandle($fh, $self->{'nodes'}, $self->{'edges'}, $placed);
  if ($file ne "-") {
    close($fh);
  }
}

sub writeNodeNames {
  my ($self, $file) = @_;
  if (!defined $file) {
    die "Undefined file for writing\n";
  }
  my $fh;
  open($fh, ">$file") || die "Can't open $file\n";
  foreach my $i (keys %{$self->{'nodes'}}) {
    my $node = $self->{'nodes'}->{$i};
    my $name = $self->{'nodeName'}->{$i};
    print $fh "$node\t$name\n";
  }
  close($fh);
}

use UnionFind;

sub mergeEdge {
  my ($self, $type) = @_;
  # Find equivalence classes
  my $uf = UnionFind->new();
  foreach my $i (keys %{$self->{'edges'}}) {
    my $ni = $self->getNode($i);
    foreach my $j (keys %{$self->{'edges'}->{$i}}) {
      my $nj = $self->getNode($j);
      my $t= $self->getType($i, $j);
      if ($t eq $type) {
        $uf->makeset($ni);
        $uf->makeset($nj);
        $uf->union($ni, $nj);
      }
    }
  }
  return $self->mergeEdgeUsingUF($uf);
}

sub mergeEdgeUsingUF {
  my ($self, $uf) = @_;
  my $nodes = {};
  foreach my $n (keys %{$self->{'nodeID'}}) {
    my $ns = $uf->findset($n);
    $ns = $n if (!defined $ns);
    if (!defined $nodes->{$ns}) {
      if ($n eq $ns) {
        $nodes->{$ns} = $self->getName($self->getIndex($n));
      }
      else {
        $nodes->{$ns} = $n.",".$self->getName($self->getIndex($n));
      }
    }
    else {
      if ($n eq $ns) {
        $nodes->{$ns} = $nodes->{$ns}.",".
          $self->getName($self->getIndex($n));
      }
      else {
        $nodes->{$ns} = $nodes->{$ns}.",".$n.",".
          $self->getName($self->getIndex($n));
      }
    }
  }
  my $res = Graph->new();
  $res->init();
  foreach (keys %{$nodes}) {
    $res->addNode($_, $nodes->{$_});
  }
  foreach my $i (keys %{$self->{'edges'}}) {
    my $ni = $self->getNode($i);
    my $nis = $uf->findset($ni);
    $nis = $ni if (!defined $nis);
    foreach my $j (keys %{$self->{'edges'}->{$i}}) {
      my $nj = $self->getNode($j);
      my $njs = $uf->findset($nj);
      $njs = $nj if (!defined $njs);
      my $t= $self->getType($i, $j);
      if ($nis ne $njs) {
        $res->addEdge($nis, $njs, $t);
      }
    }
  }
  return $res;
}

sub subGraph1 {
  my ($nodes, $edges, $hash) = @_;
  my $nodes_2 = {};
  my $edges_2 = {};
  foreach my $ni (keys %{$hash}) {
    next if (!defined $nodes->{$ni});
    $nodes_2->{$ni} = 1;
  }
  foreach my $i (keys %{$edges}) {
    next if (!defined $nodes_2->{$i});
    foreach my $j (keys %{$edges->{$i}}) {
      next if (!defined $nodes_2->{$j});
      $edges_2->{$i}->{$j} = $edges->{$i}->{$j};
    }
  }
  return ($nodes_2, $edges_2);
}

sub subGraph {
  my ($self, $hash) = @_;
  my $res = Graph->new();
  $res->init();
  my $nodes = {};
  foreach my $ni (keys %{$hash}) {
    my $i = $self->getIndex($ni);
    $nodes->{$i} = 1;
    $res->addNode($ni, $self->getName($i));
  }
  foreach my $i (keys %{$self->{'edges'}}) {
    next if (!defined $nodes->{$i});
    my $ni = $self->getNode($i);
    foreach my $j (keys %{$self->{'edges'}->{$i}}) {
      next if (!defined $nodes->{$j});
      my $nj = $self->getNode($j);
      my $t= $self->getType($i, $j);
      $res->addEdge($ni, $nj, $t);
    }
  }
  return $res;
}

sub getRevConnectedSubgraph {
  my ($self, $hash, $type, $depth) = @_;
  my $nodes = $self->getRevConnectedNodes($hash, $type, $depth);
  print STDERR "#", scalar(keys %{$nodes}), "\n";
  return $self->subGraph($nodes);
}

sub doDFS_rec {
  my ($nodes, $edges, $p, $start, $visited, $startTime, $finishTime, $cluster, $time, $rank) = @_;
  if ($visited->{$start} != 1) {
    #print "Visit:", $start, "\n";
    $visited->{$start} = 1;
    $time->[0]++;
    $startTime->{$start} = $time->[0];
    $cluster->{$start} = $p;
    my @list = keys(%{$edges->{$start}});
    if (defined $rank) {
      @list = sort { $rank->{$b} <=> $rank->{$a} } keys(%{$edges->{$start}});
    }
    foreach my $j (@list) {
      &doDFS_rec($nodes, $edges, $p, $j, $visited, $startTime, $finishTime, $cluster, $time, $rank);
    }
    $time->[0]++;
    $finishTime->{$start} = $time->[0];
  }
}

sub doDFS {
  my ($nodes, $edges, $fTime, $start, $rank) = @_;
  my $visited = {};
  my $startTime = {};
  my $finishTime = {};
  my $cluster = {};
  my $time = [];
  $time->[0] = 0;
  my @list = keys(%{$nodes});
  if (defined $fTime) {
    @list = sort { $fTime->{$b} <=> $fTime->{$a} } keys(%{$nodes});
  }
  elsif (defined $rank) {
    @list = sort { $rank->{$b} <=> $rank->{$a} } keys(%{$nodes});
  }
  if (defined $start) {
    @list = ($start);
  }
  foreach my $n (@list) {
    #print "Start: ", $n, "\n";
    &doDFS_rec($nodes, $edges, $n, $n, $visited, $startTime, $finishTime, $cluster, $time, $rank);
  }
  my $revCluster = {};
  foreach my $c (keys(%{$cluster})) {
    push @{$revCluster->{$cluster->{$c}}}, $c;
  } 
  return ($revCluster, $startTime, $finishTime);
}

sub findSCC {
  my ($nodes, $edges) = @_;
  my ($cls, $sTime, $fTime) = &doDFS($nodes, $edges);
  #print "StartTime: ";
  #foreach my $c (keys(%{$sTime})) {
  #  print "(", $c, "->", $sTime->{$c}, ") ";
  #} 
  #print "\n";
  #print "FinishTime: ";
  #foreach my $c (keys(%{$fTime})) {
  #  print "(", $c, "->", $fTime->{$c}, ") ";
  #} 
  #print "\n";
  my $revedges = &getRevEdges($edges);
  my ($cluster, $startTime, $finishTime) = &doDFS($nodes, $revedges, $fTime);
  #foreach my $i (keys %{$cluster}) {
  #  print $i, "\t", scalar(@{$cluster->{$i}}), "\n";
  #}
  return $cluster;
}

sub doBFSnode {
  my ($nodes, $edges, $start, $lengthMin, $lengthMax) = @_;
  my $visited = {};
  my $depth = {};
  my @queue;
  $visited->{$start} = 1;
  push @queue, $start;
  $depth->{$start} = 0;
  my $res = {};
  while (scalar(@queue) != 0) {
    my $i = shift(@queue);
    #print join(" ", @queue), "\n";
    #print $i, "\n";
    if ($lengthMax >=0 && $depth->{$i} > $lengthMax) {
        return $res;
    }
    if ($depth->{$i} >= $lengthMin) {
        $res->{$i} = $depth->{$i};
    }
    foreach my $j (keys(%{$edges->{$i}})) {
        if ($visited->{$j} != 1) {
          $visited->{$j} = 1;
          $depth->{$j} = $depth->{$i} + 1;
          push @queue, $j;
          #print $j, " ";
        }
    }
    #print "\n";
  }
  return $res;
}

sub getUndirEdges {
  my $edges = shift;
  my $undiredges = {} ;
  foreach my $n (keys(%{$edges})) {
    foreach my $j (keys(%{$edges->{$n}})) {
      my $w = $edges->{$n}->{$j};
      if (!defined $undiredges->{$j}) {
        $undiredges->{$j} = {};
      }
      $undiredges->{$j}->{$n} = $w;
      if (!defined $undiredges->{$n}) {
        $undiredges->{$n} = {};
      }
      $undiredges->{$n}->{$j} = $w;
    }
  }
  return $undiredges;
}

sub getRevEdges {
  my $edges = shift;
  my $revedges = {} ;
  foreach my $n (keys(%{$edges})) {
    foreach my $j (keys(%{$edges->{$n}})) {
      my $w = $edges->{$n}->{$j};
      if (!defined $revedges->{$j}) {
        $revedges->{$j} = {};
      }
      $revedges->{$j}->{$n} = $w;
    }
  }
  return $revedges;
}

sub doBFS {
  my ($nodes, $edges, $lengthMin, $lengthMax) = @_;
  my $visited = {};
  my $depth = {};
  my $revedges = &getRevEdges($edges);
  my @queue;
  my $res = {};
  foreach my $start (keys(%{$nodes})) {
    next if (defined $revedges->{$start});
    next if (defined $visited->{$start});
    #print STDERR $start, "\n";
    $visited->{$start} = 1;
    push @queue, $start;
    $depth->{$start} = 0;
    $res->{$start} = $depth->{$start};
    while (scalar(@queue) != 0) {
      my $i = shift(@queue);
      #print join(" ", @queue), "\n";
      #print STDERR $i, "\n";
      if ($lengthMax >=0 && $depth->{$i} > $lengthMax) {
          last;
      }
      if ($depth->{$i} >= $lengthMin) {
          $res->{$i} = $depth->{$i};
      }
      foreach my $j (keys(%{$edges->{$i}})) {
          if ($visited->{$j} != 1) {
            $visited->{$j} = 1;
            $depth->{$j} = $depth->{$i} + 1;
            push @queue, $j;
            #print STDERR $j, " ";
          }
      }
      #print "\n";
    }
  }
  return $res;
}

sub findPath {
  my ($nodes, $edges, $start, $end, $lengthMin, $lengthMax) = @_;
  my $visited = {};
  my $parent = {};
  my @queue;
  $visited->{$start} = 1;
  push @queue, $start;
  $parent->{$start} = "-1";
  while (scalar(@queue) != 0) {
    my $i = shift(@queue);
    #print join(" ", @queue), "\n";
    #print $i, "\n";
    if ($i eq $end) {
      my $path = [];
      my $count = 0;
      while ($i ne "-1") {
        $path = [$i, @$path];
        $i = $parent->{$i};
        $count++;
      }
      if ($count > $lengthMin && $count < $lengthMax) {
        return $path;
      }
      else {
        $i = $end;
        $visited->{$i} = 0;
        next;
      }
    }
    foreach my $j (keys(%{$edges->{$i}})) {
        if ($visited->{$j} != 1) {
          $visited->{$j} = 1;
          $parent->{$j} = $i;
          push @queue, $j;
          #print $j, " ";
        }
    }
    #print "\n";
  }
  return [];
}

sub getNodesEdges {
  my ($self, $type) = @_;
  my $nodes = $self->{'nodes'};
  my $edges = $self->{'edges'};
  if (defined $type && $type ne "All") {
    if (defined $self->{'etype'}->{$type}) {
      $edges = $self->{'etype'}->{$type};
    }
    else {
      die "Can't find $type edges\n";
    }
  }
  return ($nodes, $edges);
}

sub DFS {
  my ($self, $type, $ftime, $start, $rank) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  return &doDFS($nodes, $edges, $ftime, $self->getIndex($start), $rank);
}

sub BFSnode {
  my ($self, $type, $start, $lengthMin, $lengthMax) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  return &doBFSnode($nodes, $edges, $self->getIndex($start), 
      $lengthMin, $lengthMax);
}

# Returns hash of depth of nodes visited
sub BFS {
  my ($self, $type, $lengthMin, $lengthMax) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  return &doBFS ($nodes, $edges, $lengthMin, $lengthMax);
}

sub Path {
  my ($self, $type, $start, $end, $lengthMin, $lengthMax) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  return &findPath($nodes, $edges, $self->getIndex($start), 
      $self->getIndex($end), $lengthMin, $lengthMax);
}

sub SCC {
  my ($self, $type) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  return &findSCC($nodes, $edges);
}

sub SCCUndir {
  my ($self, $type) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  $edges = &getUndirEdges($edges);
  return &findSCC($nodes, $edges);
}

sub addTo {
  my ($edges, $edges1) = @_;
  foreach my $i (keys(%{$edges1})) {
    foreach my $j (keys(%{$edges1->{$i}})) {
      if (!defined $edges->{$i}) {
        $edges->{$i} = {};
      }
      if (!defined $edges->{$i}->{$j}) {
        $edges->{$i}->{$j} = $edges1->{$i}->{$j};
      }
    }
  }
  return $edges;
}

sub removeFrom {
  my ($edges, $edges1) = @_;
  foreach my $i (keys(%{$edges1})) {
    next if (!defined $edges->{$i});
    foreach my $j (keys(%{$edges1->{$i}})) {
      if (defined $edges->{$i}->{$j}) {
        delete $edges->{$i}->{$j};
      }
    }
  }
  return $edges;
}

sub doCopy {
  my $edges1 = shift;
  my $edges = {};
  foreach my $i (keys(%{$edges1})) {
    $edges->{$i} = {};
    foreach my $j (keys(%{$edges1->{$i}})) {
      $edges->{$i}->{$j} = $edges1->{$i}->{$j};
    }
  }
  return $edges;
}

sub doCompose {
  my ($edges1, $edges2) = @_;
  my $edges = {};
  foreach my $i (keys(%{$edges1})) {
    foreach my $j (keys(%{$edges1->{$i}})) {
      foreach my $k (keys(%{$edges2->{$j}})) {
        if (!defined $edges->{$i}) {
          $edges->{$i} = {};
        }
        $edges->{$i}->{$k} = "o";
      }
    }
  }
  return $edges;
}

sub doUnion {
  my ($edges1, $edges2) = @_;
  my $edges = &doCopy($edges2);
  foreach my $i (keys(%{$edges1})) {
    foreach my $j (keys(%{$edges1->{$i}})) {
      if (!defined $edges->{$i}) {
        $edges->{$i} = {};
      }
      $edges->{$i}->{$j} = $edges1->{$i}->{$j};
    }
  }
  return $edges;
}

# Based on logarithmic algorithms of Ioannidis, Valduriez and Boral
sub doClosure {
  my $edges = shift;
  my $rplus = &doCopy($edges);
  my $del = &doCopy($edges);
  my $delta = &doCopy($edges);
  my $index = 0;
  while ( scalar(keys %{$del}) != 0) {
    print STDERR "TC $index\n";
    $delta = &doCompose($delta, $delta);
    $del = &doCompose($rplus, $delta);
    $rplus = &doUnion($rplus, $del);
    $rplus = &doUnion($rplus, $delta);
    $index++;
  }
  return $rplus;
}

sub transitiveClosure {
  my ($self, $type) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  my $rplus = &doClosure($edges);
  &addTo($edges, $rplus);
  if (defined $type && $type ne "All") {
    &addTo($self->{'edges'}, $rplus);
  }
}

sub getUFfromCluster {
  my ($self, $cluster) = @_;
  my $uf = UnionFind->new();
  foreach my $i (keys %{$cluster}) {
    my $ni = $self->getNode($i);
    $uf->makeset($ni);
    foreach my $j (@{$cluster->{$i}}) {
      my $nj = $self->getNode($j);
      $uf->makeset($nj);
      $uf->union($ni, $nj);
    }
  }
  return $uf;
}

sub getGraphFromCluster {
  my ($self, $cluster, $type) = @_;
  my $g = Graph->new();
  $g->init();
  foreach my $i (keys %{$cluster}) {
    my $ni = $self->getNode($i);
    $g->addNode($ni, $self->getName($i));
    foreach my $j (@{$cluster->{$i}}) {
      my $nj = $self->getNode($j);
      $g->addNode($nj, $self->getName($j));
      $g->addEdge($ni, $nj, $type);
    }
  }
  return $g;
}

sub compactSCC {
  my ($self, $type) = @_;
  my $components = $self->SCC($type);
  my $uf = $self->getUFfromCluster($components);
  return $self->mergeEdgeUsingUF($uf);
}

sub SCCUndirComponents {
  my ($self, $type) = @_;
  my $components = $self->SCCUndir($type);
  my $g = $self->getGraphFromCluster($components, $type);
  return $g;
}

sub doReduce {
  my $edges = shift;
  my $rplus = &doClosure($edges);
  my $del = &doCompose($edges, $rplus);
  &removeFrom($edges, $del);
  return $del;
}

sub transitiveReduce {
  my ($self, $type) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  my $del = &doReduce($edges);
  if (defined $type && $type ne "All") {
    &removeFrom($self->{'edges'}, $del);
  }
}

sub getRoot {
   my ($revedges, $depth, $n) = @_;
   while ($depth->{$n} != 0) {
      if (defined $revedges->{$n}) {
        my @list = keys %{$revedges->{$n}};
        if (scalar(@list) > 0) {
           $n = $list[0]; 
           next;
        }
      }
      die " Error in getRoot \n";
   }
   return $n;
}

sub getLeaf {
   my ($edges, $depth, $n) = @_;
   my $res = {};
   if (!defined $depth->{$n}) {
     die "Error depth not defined for $n\n";
   }
   if (!defined $edges->{$n}) {
     $res->{$n} = 1;
     return;
   }
   my $d = $depth->{$n};
   my $md = 0;
   foreach my $i (keys(%{$edges->{$n}})) {
     if ($md < $depth->{$i}) {
        $md = $depth->{$i};
     }
   }
   if ($d >= $md) {
     $res->{$n} = 1;
     return;
   }
   foreach my $i (keys(%{$edges->{$n}})) {
     my $list = &getLeaf($edges, $depth, $i);
     foreach $j (keys %{$list}) {
       $res->{$j} = 1;
     }
   }
   return $res;
}

sub getChildrenDeep {
   my ($edges, $n, $depthLimit) = @_;
   my $res = {};
   return $res if (!defined $edges->{$n});
   my $visited = {};
   my $depth = {};
   $depth->{$n} = 0;
   my @queue;
   push @queue, $n;
   while (scalar(@queue) > 0) {
     my $e = shift @queue;
     next if (defined $depthLimit && $depth->{$e} > $depthLimit);
     if (!defined $visited->{$e}) {
       $visited->{$e} = 1;
       foreach my $i (keys(%{$edges->{$e}})) {
         if (!defined $visited->{$i}) {
           $depth->{$i} = $depth->{$e} + 1;
         }
         push @queue, $i;
       }
     }
   }
   return $visited;
}

sub getChildren {
   my ($edges, $list) = @_;
   my $res = {};
   foreach my $i (keys(%{$list})) {
     next if (!defined $edges->{$i});
     foreach my $j (keys(%{$edges->{$i}})) {
       $res->{$j} = 1;
     }
   }
   return $res;
}

sub updatePlacing {
   my ($edges, $r, $d, $rx, $ry, $placed, $wx, $wy, $depth) = @_;
   return if (defined $placed->{$r});
   $placed->{$r} = [($rx + $d * $wx), $ry];
   print STDERR "$r $d $rx $ry $wx $wy\n";
   my $ch = &getChildren($edges, {$r => 1});
   my @list = keys(%{$ch});
   for (my $i =0; $i < scalar(@list); $i++) {
     next if ($depth->{$r} >= $depth->{$list[$i]});
     my $y = $ry + ($i - scalar(@list)/2) * $wy;
     &updatePlacing($edges, $list[$i], $d - 1, $rx, $y, $placed, $wx, $wy, $depth);
   }
}

sub upperRound {
  my $val = shift;
  return (int($val) != $val) ? int($val) + 1 : $val;
}

sub placeSideTree {
  my ($self, $refx, $refy, $widthx, $widthy, $nodeList, $type, $depth) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  my $revedges = &getRevEdges($edges);
  my @list = sort { $depth->{$b} <=> $depth->{$a} } @$nodeList;
  my $placed = {};
  foreach my $n (@list) {
    next if (defined $placed->{$n});
    my $d = $depth->{$n};
    if ($d == 0) {
      $placed->{$n} = [$refx, $refy];
      $refy = $refy + $widthy;
      next;
    }
    my $r = &getRoot($revedges, $depth, $n);
    my $ch = &getLeaf($edges, $depth, $r);
    foreach my $i (keys(%{$ch})) {
      if ($d < $depth->{$i}) {
        $d = $depth->{$i};
      }
    }
    my $num = scalar(keys %{$ch});
    my $rx = $refx;
    my $ry = $refy + &upperRound($num/2.0) * $widthy;
    print STDERR "Root $r $d $rx $ry $widthx $widthy\n";
    &updatePlacing($edges, $r, $d, $rx, $ry, $placed, $widthx, $widthy, $depth);
    $refy = $refy + $num * $widthy;
  }
  return $placed;
}

sub createGraphEdges {
  my ($self, $type) = @_;
  my ($nodes, $edges) = $self->getNodesEdges($type);
  my $res = Graph->new();
  $res->init();
  foreach my $i (keys(%{$edges})) {
    my $ni = $self->getNode($i);
    $res->addNode($ni, $self->getName($i));
    foreach my $j (keys(%{$edges->{$i}})) {
      my $nj = $self->getNode($j);
      $res->addNode($nj, $self->getName($j));
      $res->addEdge($ni, $nj, $type);
    }
  }
  return $res;
}

sub getConnectedNodes {
  my ($self, $hash, $type, $depth) = @_;
  my $res = {};
  foreach (keys %{$hash}) {
    $res->{$self->getIndex($_)} = 1;
  }
  my ($nodes, $edges) = $self->getNodesEdges($type);
  my $fres = {};
  foreach (keys %{$res}) {
    if (!defined $fres->{$_}) {
      my $h = &getChildrenDeep($edges, $_, $depth);
      foreach (keys %{$h}) {
        $fres->{$_} = $h->{$_};
      }
    }
  }
  my $ffres = {};
  foreach (keys %{$fres}) {
    my $n = $self->getNode($_);
    $ffres->{$n} = $fres->{$_};
  }
  return $ffres;
}

sub getRevConnectedNodes {
  my ($self, $hash, $type, $depth) = @_;
  my $res = {};
  foreach (keys %{$hash}) {
    $res->{$self->getIndex($_)} = 1;
  }
  my ($nodes, $edges) = $self->getNodesEdges($type);
  $edges = &getRevEdges($edges);
  my $fres = {};
  foreach (keys %{$res}) {
    if (!defined $fres->{$_}) {
      my $h = &getChildrenDeep($edges, $_, $depth);
      foreach (keys %{$h}) {
        $fres->{$_} = $h->{$_};
      }
    }
  }
  my $ffres = {};
  foreach (keys %{$fres}) {
    my $n = $self->getNode($_);
    $ffres->{$n} = $fres->{$_};
  }
  return $ffres;
}

sub getOtherNodes {
  my ($self, $hash) = @_;
  my $res = {};
  my $nodes = $self->{'nodes'};
  foreach (keys %{$nodes}) {
    my $n = $self->getNode($_);
    if (!defined $hash->{$n}) {
      $res->{$n} = $_;
    }
  }
  return $res;
}

sub intersectHash {
  my ($a, $b) = @_;
  my $hash = {};
  foreach (keys %{$a}) {
    if (defined $b->{$_}) {
      $hash->{$_} = $a->{$_};
    }
  }
  return $hash;
}

sub traverseTwo {
  my ($self, $type1, $type2) = @_;
  my ($nodes, $edges1) = $self->getNodesEdges($type1);
  my ($nodes2, $edges2) = $self->getNodesEdges($type2);
  foreach my $s (keys(%{$nodes2})) {
    $nodes->{$s} = $nodes2->{$s};
  }
  my $visited = {};
  my $depth = {};
  my $parent = {};
  my $otherNodes = {};
  my $revedges = &getRevEdges($edges1);
  my @queue;
  foreach my $start (keys(%{$nodes})) {
    next if (defined $revedges->{$start});
    next if (defined $visited->{$start});
    #print "v ", $start, "\n";
    $visited->{$start} = 1;
    $depth->{$start} = 0;
    next if (!defined $edges2->{$start});
    $otherNodes->{$start} = $edges2->{$start};
    push @queue, $start;
    while (scalar(@queue) != 0) {
      my $i = shift(@queue);
      #print "q ", join(" ", @queue), "\n";
      #print "v ", $i, " ", scalar(keys %{$otherNodes->{$i}}), "\n";
      foreach my $j (keys(%{$edges1->{$i}})) {
        next if (!defined $edges2->{$j});
        my $on = &intersectHash($otherNodes->{$i},$edges2->{$j});
        #print $self->getNode($j), ":", scalar(keys %{$on}), ":", scalar(keys %{$otherNodes->{$j}}), " ";
        if ($visited->{$j} != 1) {
          next if scalar(keys %{$on}) <= 0;
          $otherNodes->{$j} = $on;
          $parent->{$j} = $i;
          $visited->{$j} = 1;
          $depth->{$j} = $depth->{$i} + 1;
          push @queue, $j;
          #print $j, " ";
        }
        elsif (scalar(keys %{$on}) > scalar(keys %{$otherNodes->{$j}})) {
          next if scalar(keys %{$on}) <= 0;
          $otherNodes->{$j} = $on;
          $parent->{$j} = $i;
          $visited->{$j} = 1;
          $depth->{$j} = $depth->{$i} + 1;
          push @queue, $j;
          #print $j, " ";
        }
      }
      #print "\n";
    }
  }
  my $res = Graph->new();
  $res->init();
  #print "Final : \n";
  foreach my $i (keys %{$depth}) {
    next if ($depth->{$i} <= 1);
    $res->addNode($self->getNode($i), $self->getName($i));
    #print $self->getNode($i);
    my $k = $i;
    while (defined $parent->{$i}) {
      $res->addNode($self->getNode($parent->{$i}), $self->getName($parent->{$i}));
      $res->addEdge($self->getNode($parent->{$i}), $self->getNode($i), $type1);
      $i = $parent->{$i};
      #print " ", $self->getNode($i);
    }
    #print " hilo ";
    foreach my $n (keys %{$otherNodes->{$k}}) {
      $res->addNode($self->getNode($n), $self->getName($n));
      $res->addEdge($self->getNode($k), $self->getNode($n), $type2);
      $res->addEdge($self->getNode($n), $self->getNode($k), $type2);
      #print " ", $self->getNode($n);
    }
    #print "\n";
  }
  return $res;
}

sub clusterTwo {
  my ($self, $type1, $type2) = @_;
  my ($nodes, $edges1) = $self->getNodesEdges($type1);
  my ($nodes2, $edges2) = $self->getNodesEdges($type2);
  foreach my $s (keys(%{$nodes2})) {
    $nodes->{$s} = $nodes2->{$s};
  }
  my $visited = {};
  my $depth = {};
  my $parent = {};
  my $set = UnionFind->new();
  my $revedges = &getRevEdges($edges1);
  my @queue;
  foreach my $start (keys(%{$nodes})) {
    next if (defined $visited->{$start});
    #print "v ", $start, "\n";
    $visited->{$start} = 1;
    $depth->{$start} = 0;
    push @queue, $start;
    while (scalar(@queue) != 0) {
      my $i = shift(@queue);
      $set->makeset($i);
      #print "q ", join(" ", @queue), "\n";
      #print "v ", $i, " ", scalar(keys %{$otherNodes->{$i}}), "\n";
      foreach my $j (keys(%{$edges1->{$i}})) {
          if ($visited->{$j} != 1) {
            $set->makeset($j);
            my $ec2 = 0;
            if (defined $edges2->{$j}) {
              my $s = $set->findset($i);
              foreach $k (keys %{$edges2->{$j}}) {
                if ($s eq $set->findset($k)) {
                  $ec2 = 1;
                  last;
                }
              }
            }
            if ($ec2 == 0) {
              $set->union($i, $j);
            }
            $parent->{$j} = $i;
            $visited->{$j} = 1;
            $depth->{$j} = $depth->{$i} + 1;
            push @queue, $j;
            #print $j, " ";
          }
      }
      #print "\n";
    }
  }
  return $set;
}

#Find Topological Order
sub topologicalOrder {
  my ($nodes, $edges) = @_;
  my ($cls, $sTime, $fTime) = &doDFS($nodes, $edges);
  my @list = sort { $fTime->{$b} <=> $fTime->{$a} } keys(%{$nodes});
  return [@list];
}

#Find Longest Path in a DAG
sub findLongestPaths {
  my ($nodes, $edges) = @_;
  my $length_to = {};
  my $predecessor = {};
  foreach (keys %{$nodes}) {
    $length_to->{$_} = 0;
  }
  my $topOrder = &topologicalOrder($nodes, $edges);
  foreach my $v (@$topOrder) {
    foreach my $w (keys(%{$edges->{$v}})) {
      if ($length_to->{$w} <= ($length_to->{$v} + 1)) {
        $length_to->{$w} = $length_to->{$v} + 1;
        $predecessor->{$w} = $v;
      }
    }
  }
  my $paths = [];
  my @keys = @$topOrder;
  while (scalar(@keys) > 0) {
    my $max = -1;
    my $maxIndex = -1;
    foreach my $v (@keys) {
      if ($max < $length_to->{$v}) {
        $max = $length_to->{$v};
        $maxIndex = $v;
      }
    }
    my $res = [];
    my $h = {};
    my $v = $maxIndex;
    while (defined $v) {
      last if (defined $h->{$v});
      $h->{$v} = 1;
      unshift @$res, $v;
      $v = $predecessor->{$v};
    }
    push @$paths, $res;
    my @list;
    foreach my $v (@keys) {
      next if (defined $h->{$v});
      push @list, $v;
    }
    @keys = @list;
  }
  return $paths;
}

sub getDeepPath {
  my ($nodes, $edges, $start, $rank) = @_;
  my $visited = {};
  my @stack = ($start);
  $visited->{$start} = 1;
  my $path = [];
  while (scalar(@stack) != 0) {
    my $last = pop(@stack);
    push @$path, $last;
    my @list = keys(%{$edges->{$last}});
    if (defined $rank) {
      @list = sort { $rank->{$b} <=> $rank->{$a} } @list;
    }
    if (scalar(@list) > 0) {
      foreach my $id (@list) {
        if (!defined $visited->{$id}) {
          $visited->{$id} = 1;
          push @stack, $id;
          last;
        }
      }
    }
  }
  return $path;
}

1;
__END__

=head1 NAME

Graph - A Graph read/write/maniputation interface 

Example:

  use Graph;

  my $n = Graph->new(file => "A.sif");
  $n->init();
  $n->transitiveReduce();
  $n->write("B.sif");

=cut
