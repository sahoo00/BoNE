# Author: Debashis Sahoo <dsahoo@ucsd.edu>
package U;

sub uniq {
  my ($a, $rest) = @_;
  my %seen;
  return [grep { !$seen{$_}++ } @$a];
}

sub a_hash {
  my ($a, $rest) = @_;
  my $hash = {};
  foreach (@$a) {
    $hash->{$_} = 1;
  }
  return $hash;
}

sub intersection {
  my ($a, @list) = @_;
  my $hash = [map { &a_hash($_) } ($a, @list)];
  my $res = [];
  foreach my $i (@$a) {
    my $all = 1;
    foreach my $h (@$hash) {
      if (!defined $h->{$i}) {
        $all = 0;
      }
    }
    if ($all) {
      push @$res, $i;
    }
  }
  return $res;
}

sub union {
  my $hash = {};
  my $res = [];
  foreach my $a (@_) {
    foreach (@$a) {
      if (!defined $hash->{$_}) {
        $hash->{$_} = 1;
        push @$res, $_;
      }
    }
  }
  return $res;
}

sub diff {
  my ($a, $b) = @_;
  my $hash = {};
  foreach (@$b) {
    if (!defined $hash->{$_}) {
        $hash->{$_} = 1;
    }
  }
  my $res = [];
  foreach (@$a) {
    if (!defined $hash->{$_}) {
        $hash->{$_} = 1;
        push @$res, $_;
    }
  }
  return $res;
}

sub venn {
  my ($a, $b) = @_;
  return [&diff($a, $b), &intersection($a, $b), &diff($b, $a)];
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

sub unionHash {
  my $hash = {};
  foreach my $a (@_) {
    foreach (keys %{$a}) {
      $hash->{$_} = $a->{$_};
    }
  }
  return $hash;
}

sub diffHash {
  my ($a, $b) = @_;
  my $hash = {};
  foreach (keys %{$a}) {
    if (!defined $b->{$_}) {
      $hash->{$_} = $a->{$_};
    }
  }
  return $hash;
}

sub covariance {
   my ($array1ref, $array2ref) = @_;
   my ($i, $result);
   for ($i = 0; $i < @$array1ref; $i++) {
       $result += $array1ref->[$i] * $array2ref->[$i];
   }
   $result /= @$array1ref;
   $result -= mean($array1ref) * mean($array2ref);
   return $result;
}

sub correlation {
   my ($array1ref, $array2ref) = @_;
   my ($sum1, $sum2);
   my ($sum1_squared, $sum2_squared); 
   foreach (@$array1ref) { $sum1 += $_;  $sum1_squared += $_**2; }
   foreach (@$array2ref) { $sum2 += $_;  $sum2_squared += $_**2; }
   my $numerator = (@$array1ref**2) * covariance($array1ref, $array2ref);
   my $denominator = sqrt(((@$array1ref * $sum1_squared) - ($sum1**2)) *
                          ((@$array1ref * $sum2_squared) - ($sum2**2)));
   my $r;
   if ($denominator == 0) {
     print STDERR "The denominator is 0.\n";
     return "";
   } else {
      $r = $numerator / $denominator;
   }
    return $r;
}

sub mse {
  my ($arrayref, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$arrayref) - 1;
  }
  my $result = 0.0;
  if ($start > $end) {
    return $result;
  }
  my $m = &mean($arrayref, $start, $end);
  my $result = 0.0;
  for (my $i = $start; $i <= $end; $i++) {
    $result += ($arrayref->[$i] - $m) * ($arrayref->[$i] - $m);
  }
  return $result;
}

sub max {
  my ($arrayref, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$arrayref) - 1;
  }
  my $result = $arrayref->[$start];
  if ($start > $end) {
    return $result;
  }
  for (my $i = $start; $i <= $end; $i++) {
    if ($result < $arrayref->[$i]) {
      $result = $arrayref->[$i];
    }
  }
  return $result;
}

sub min {
  my ($arrayref, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$arrayref) - 1;
  }
  my $result = $arrayref->[$start];
  if ($start > $end) {
    return $result;
  }
  for (my $i = $start; $i <= $end; $i++) {
    if ($result > $arrayref->[$i]) {
      $result = $arrayref->[$i];
    }
  }
  return $result;
}

sub median {
  my ($arrayref, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$arrayref) - 1;
  }
  my $result = 0.0;
  if ($start > $end) {
    return $result;
  }
  $result = $arrayref->[int(($start+$end)/2)];
  return $result;
}

sub mean {
  my ($arrayref, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$arrayref) - 1;
  }
  my $result = 0.0;
  if ($start > $end) {
    return $result;
  }
  for (my $i = $start; $i <= $end; $i++) {
    $result += $arrayref->[$i];
  }
  return $result / ($end - $start + 1);
}

sub stddev {
  my ($arrayref, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$arrayref) - 1;
  }
  my $sum = 0.0;
  if ($start >= $end) {
    return $sum;
  }
  my $sumsq = 0;
  for (my $i = $start; $i <= $end; $i++) {
    $sum += $arrayref->[$i];
    $sumsq += $arrayref->[$i] * $arrayref->[$i];
  }
  my $m = $sum / ($end - $start + 1);
  my $msq = $sumsq / ($end - $start + 1);
  my $factor = ($end - $start + 1) / ($end - $start);
  my $res = $factor * ($msq - $m * $m);
  if ($res < 0) { $res = 0; }
  return sqrt($res);
}

sub variance {
  my ($arrayref, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$arrayref) - 1;
  }
  my $sum = 0.0;
  if ($start >= $end) {
    return $sum;
  }
  my $sumsq = 0;
  for (my $i = $start; $i <= $end; $i++) {
    $sum += $arrayref->[$i];
    $sumsq += $arrayref->[$i] * $arrayref->[$i];
  }
  my $m = $sum / ($end - $start + 1);
  my $msq = $sumsq / ($end - $start + 1);
  my $factor = ($end - $start + 1) / ($end - $start);
  return ($msq - $m * $m);
}

sub sum {
  my ($arrayref, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$arrayref) - 1;
  }
  my $result = 0.0;
  for (my $i = $start; $i <= $end; $i++) {
    $result += $arrayref->[$i];
  }
  return $result;
}

sub info {
  my ($arrayref, $start, $end) = @_;
  my $min = &min($arrayref, $start, $end);
  my $max = &max($arrayref, $start, $end);
  my $m = &mean($arrayref, $start, $end);
  my $s = &stddev($arrayref, $start, $end);
  my $sum = &sum($arrayref, $start, $end);
  return [$min, $m, $max, $stddev, $sum];
}

sub fitStep {
  my ($data, $start, $end, $times) = @_;
  my $count = $end - $start + 1;
  if ($count <= 0) {
    return [0, 0, 0, 0, 0, 0, 0, 0];
  }
  if (!defined $times) {
    $times = [0 .. $end];
  }
  else {
    my @sk = sort { $times->[$a] <=> $times->[$b] } $start .. $end;
    $data = [map { $data->[$_] } @sk];
    $times = [map { $times->[$_] } @sk];
    $start = 0;
    $end = $count - 1;
  }
  my @sseArray = map { 0.0 } 0 .. ($count - 1);
  my $sum = &sum($data, $start, $end);
  my $mean = &mean($data, $start, $end);
  my $sstot = &mse($data, $start, $end);
  my $sum1 = 0.0;
  my $count1 = 0;
  my $m1 = 0.0;
  my $sum2 = $sum;
  my $count2 = $count;
  my $m2 = ($sum/$count);
  my $sum1sq = 0.0;
  my $sum2sq = $sstot;
  my $sse = $sum1sq + $sum2sq;
  for (my $i = 0; $i < $count; $i++) {
    my $entry = $data->[$i + $start];
    if (!defined $entry || $entry eq "") {
      $sseArray[$i] = $sse;
      next;
    }
    $count1 ++;
    $count2 --;
    if ($count2 == 0) {
      $sseArray[$i] = $sstot;
      next;
    }
    my $tmp = ($mean - ($entry + $sum1)/$count1);
    $sum1sq = $sum1sq + ($entry - $mean) * ($entry - $mean) - 
      $tmp * $tmp * $count1 + ($count1 - 1) * ($mean - $m1) * ($mean - $m1);
    $tmp = ($mean - ($sum2 - $entry)/$count2);
    $sum2sq = $sum2sq - ($entry - $mean) * ($entry - $mean) - 
      $tmp * $tmp * $count2 + ($count2 + 1) * ($mean - $m2) * ($mean - $m2);
    $sum1 += $entry;
    $sum2 -= $entry;
    $m1 = $sum1/$count1;
    $m2 = $sum2/$count2;
    $sse = $sum1sq + $sum2sq;
    $sseArray[$i] = $sse;
  }
  #print join(",", map { sprintf("%.2f", $_) } @sseArray), "\n-----------\n";
  my $bestSse;
  my $bestIndex = 0;
  my $dof = 0;
  for (my $i = 0; $i < $count ; $i++) {
    if ($i < ($count - 1)) {
       next if ($times->[$i+$start] eq $times->[$i+$start+1]);
    }
    $dof++;
    my $index = $i + $start;
    if (!defined $bestSse) {
      $bestSse = $sseArray[$i];
      $bestIndex = $index;
    }
    if ($sseArray[$i] < $bestSse) {
      $bestSse = $sseArray[$i];
      $bestIndex = $index;
    }
    #{# Debug code start
    #  $entry = $data->[$index];
    #  $sum1sq = &mse($data, $start, $index);
    #  $sum2sq = &mse($data, $index + 1, $end);
    #  $sse = $sum1sq + $sum2sq;
    #  printf "%.2f\t%.2f\t%.2f\t%.2f\n", $sseArray[$i], $entry, $sse, $sum1sq;
    #}# Debug code end
  }
  #printf "bestSse = %.2f\t bestIndex=%d\n", $bestSse, $bestIndex;
  $m1 = &mean($data, $start, $bestIndex);
  $m2 = &mean($data, $bestIndex + 1, $end);
  my $thr = ($m1 + $m2)/2.0;
  #{# Debug code start
  #  $sum1sq = &mse($data, $start, $bestIndex);
  #  $sum2sq = &mse($data, $bestIndex + 1, $end);
  #  $sse = $sum1sq + $sum2sq;
  #  if ($sse != $bestSse) {
  #    print STDERR "SSE calculation is wrong :", $sse, ":", $bestSse, "\n";
  #  }
  #}# Debug code end
  my $label = 0;
  if ($m1 < $m2) {
    $label = 1;
  }
  else {
    $label = 2;
  }
  my $statistic = 0;
  if ($bestSse > 0) {
    if ($dof > 4) {
      $statistic = ($sstot - $bestSse)/3/($bestSse/($dof - 4));
    }
    else {
      $statistic = ($sstot - $bestSse)/2/$bestSse;
    }
  }
  return [$bestIndex, $bestSse, $sstot, $statistic, $m1, $m2, $thr, $label];
}

sub fitStepOld {
  my ($data, $start, $end) = @_;
  my $count = $end - $start + 1;
  if ($count <= 0) {
    return [0, 0, 0, 0, 0, 0, 0, 0];
  }
  my @sseArray = map { 0.0 } 0 .. ($count - 1);
  my $sum = &sum($data, $start, $end);
  my $mean = &mean($data, $start, $end);
  my $sstot = &mse($data, $start, $end);
  my $sum1 = 0.0;
  my $count1 = 0;
  my $m1 = 0.0;
  my $sum2 = $sum;
  my $count2 = $count;
  my $m2 = ($sum/$count);
  my $sum1sq = 0.0;
  my $sum2sq = $sstot;
  my $sse = $sum1sq + $sum2sq;
  for (my $i = 0; $i < $count; $i++) {
    my $entry = $data->[$i + $start];
    if (!defined $entry || $entry eq "") {
      $sseArray[$i] = $sse;
      next;
    }
    $count1 ++;
    $count2 --;
    if ($count2 == 0) {
      $sseArray[$i] = $sstot;
      next;
    }
    my $tmp = ($mean - ($entry + $sum1)/$count1);
    $sum1sq = $sum1sq + ($entry - $mean) * ($entry - $mean) - 
      $tmp * $tmp * $count1 + ($count1 - 1) * ($mean - $m1) * ($mean - $m1);
    $tmp = ($mean - ($sum2 - $entry)/$count2);
    $sum2sq = $sum2sq - ($entry - $mean) * ($entry - $mean) - 
      $tmp * $tmp * $count2 + ($count2 + 1) * ($mean - $m2) * ($mean - $m2);
    $sum1 += $entry;
    $sum2 -= $entry;
    $m1 = $sum1/$count1;
    $m2 = $sum2/$count2;
    $sse = $sum1sq + $sum2sq;
    $sseArray[$i] = $sse;
  }
  #print join(",", map { sprintf("%.2f", $_) } @sseArray), "\n-----------\n";
  my $bestSse;
  my $bestIndex = 0;
  for (my $i = 0; $i < $count ; $i++) {
    my $index = $i + $start;
    if (!defined $bestSse) {
      $bestSse = $sseArray[$i];
      $bestIndex = $index;
    }
    if ($sseArray[$i] < $bestSse) {
      $bestSse = $sseArray[$i];
      $bestIndex = $index;
    }
    #{# Debug code start
    #  $entry = $data->[$index];
    #  $sum1sq = &mse($data, $start, $index);
    #  $sum2sq = &mse($data, $index + 1, $end);
    #  $sse = $sum1sq + $sum2sq;
    #  printf "%.2f\t%.2f\t%.2f\t%.2f\n", $sseArray[$i], $entry, $sse, $sum1sq;
    #}# Debug code end
  }
  #printf "bestSse = %.2f\t bestIndex=%d\n", $bestSse, $bestIndex;
  $m1 = &mean($data, $start, $bestIndex);
  $m2 = &mean($data, $bestIndex + 1, $end);
  my $thr = ($m1 + $m2)/2.0;
  #{# Debug code start
  #  $sum1sq = &mse($data, $start, $bestIndex);
  #  $sum2sq = &mse($data, $bestIndex + 1, $end);
  #  $sse = $sum1sq + $sum2sq;
  #  if ($sse != $bestSse) {
  #    print STDERR "SSE calculation is wrong :", $sse, ":", $bestSse, "\n";
  #  }
  #}# Debug code end
  my $label = 0;
  if ($m1 < $m2) {
    $label = 1;
  }
  else {
    $label = 2;
  }
  my $statistic = 0;
  if ($bestSse > 0) {
    if ($count > 4) {
      $statistic = ($sstot - $bestSse)/3/($bestSse/($count - 4));
    }
    else {
      $statistic = ($sstot - $bestSse)/2/$bestSse;
    }
  }
  return [$bestIndex, $bestSse, $sstot, $statistic, $m1, $m2, $thr, $label];
}

sub getStepMinerThr {
  my ($data, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$data) - 1;
  }
  my @array;
  for ($i = $start; $i <= $end; $i++) {
    next if (!defined $data->[$i] || $data->[$i] eq "");
    push @array, $data->[$i];
  }
  my @sorted = sort { $a <=> $b } @array;
  return &fitStep(\@sorted, 0, $#sorted);
}

sub printThr {
  my ($file, $start, $end, $gap) = @_;
  my $fh;
  open($fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  if ($end > $#headers) {
    $end = $#headers;
  }
  my $index = 0;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my @data = sort { $a <=> $b } @list[ $start .. $end ];
    my $res = &fitStep(\@data, 0, $#data);
    print join("\t", $list[0], $res->[6], $res->[3], $res->[6]-$gap,
        $res->[6]+$gap), "\n";
    if ( ($index % 1000) == 0 ) {
      print STDERR $index, "\n";
    }
    $index++;
  }
  close($fh);
}

sub convertpcl {
  my ($file, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  my $fh;
  open($fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  if (!defined $end || $end > $#headers) {
    $end = $#headers;
  }
  print join("\t", @headers[0, 1]), "\tGWEIGHT\t", join("\t", @headers[$start .. $end]), "\n";
  print join("\t", "EWEIGHT", ""), "\t1\t", join("\t", map { 1 } @headers[$start .. $end]), "\n";
 
  my $index = 0;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
  print join("\t", @list[0, 1]), "\t1\t", join("\t", @list[$start .. $end]), "\n";
    if ( ($index % 1000) == 0 ) {
      print STDERR $index, "\n";
    }
    $index++;
  }
  close($fh);
}

sub convertexpr {
  my ($file, $start, $end) = @_;
  if (!defined $start) {
    $start = 3;
  }
  my $fh;
  open($fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  if (!defined $end || $end > $#headers) {
    $end = $#headers;
  }
  print join("\t", @headers[0, 1]), "\t", join("\t", @headers[$start .. $end]), "\n";
 
  my $index = 0;
  while (<$fh>) {
    next if (/^EWEIGHT/);
    s/[\r\n]//g;
    my @list = split("\t");
    print join("\t", @list[0, 1]), "\t", join("\t", @list[$start .. $end]), "\n";
    if ( ($index % 1000) == 0 ) {
      print STDERR $index, "\n";
    }
    $index++;
  }
  close($fh);
}

sub convertidx {
  my $file = shift;
  open(FL, "< $file") || die "Can't open $file\n";
  my $index = 0;
  my $ptr = tell(FL);
  while (<FL>) {
    my @list = split("\t", $_);
    my ($name, $desc) = split(":", $list[1], 2);
    if ($name =~ /^\s*$/) {
       $name = "---";
    }
    print join("\t", $list[0], $ptr, $name, $desc), "\n";
    $index++;
    $ptr = tell(FL);
  }
  close(FL);
}

sub convertih {
  my $file = shift;
  open(FL, "< $file") || die "Can't open $file\n";
  my $header=<FL>; 
  $header =~ s/[\r\n]//g;
  my @headers = split(/\t/, $header);
  print "ArrayID\tArrayHeader\tClinicalhHeader\n";
  for (my $i= 2; $i < scalar(@headers); $i++) {
    my $arr = $headers[$i];
    print "$arr\t$arr\t$arr\n";
  }
  close(FL);
}

sub convertsurv {
  my $file = shift;
  open(FL, "< $file") || die "Can't open $file\n";
  my $header=<FL>; 
  $header =~ s/[\r\n]//g;
  my @headers = split(/\t/, $header);
  print "ArrayID\ttime\tstatus\tc Title\n";
  for (my $i= 2; $i < scalar(@headers); $i++) {
    my $arr = $headers[$i];
    print "$arr\t\t\t$arr\n";
  }
  close(FL);
}

sub getHeaderFh {
  my $fh = shift;
  seek($fh, 0, 0);
  my $header=<$fh>; 
  $header =~ s/[\r\n]//g;
  my @headers = split(/\t/, $header);
  return [@headers];
}

sub getHeaders {
    my $pclfile = shift;
    my $fh1;
    open($fh1, "<$pclfile") || die "Can't open $pclfile";
    my $header=<$fh1>; 
    $header =~ s/[\r\n]//g;
    my @headers = split(/\t/, $header);
    close($fh1);
    return [@headers];
}

sub getExprChange {
  my ($file, $group1, $group2) = @_;
  open(my $fh, "<$file") || die "Can't open $file";
  my $header=<$fh>; 
  $header =~ s/[\r\n]//g;
  my @headers = split(/\t/, $header);
  my $res = {};
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split(/\t/);
    my ($name, $desc) = split(/: /, $list[1]);
    my $m1 = &mean([map {$list[$_]} @$group1]);
    my $m2 = &mean([map {$list[$_]} @$group2]);
    $res->{$list[0]} = [$m2-$m1, $name, $desc]
  }
  close($fh);
  return $res;
}

sub getHash {
  my ($file, $index, $spstr, $spindex, $rest) = @_;
  if (!defined $index) {
    $index = 0;
  }
  my $hash = {};
  my $fh;
  if ($file eq "-") {
    $fh = \*STDIN;
  }
  else {
    open($fh, "<$file") || die "Can't open $file\n";
  }
  while (my $in = <$fh>) {
    $in =~ s/[\r\n]//g;
    my @list = split("\t", $in, -1);
    my $id = $list[$index];
    if (defined $spstr) {
      $spindex = 0 if (!defined $spindex);
      my @l = split($spstr, $list[$index]);
      $id = $l[$spindex];
      $id =~ s/\s//g;
    }
    $hash->{$id} = [@list];
  }
  if ($file ne "-") {
    close($fh);
  }
  return $hash;
}

sub getHashAndOrder {
  my ($file, $index, $spstr, $spindex, $rest) = @_;
  if (!defined $index) {
    $index = 0;
  }
  my $fh;
  my $hash = {};
  my $order = [];
  open($fh, "<$file") || die "Can't open $file\n";
  while (my $in = <$fh>) {
    $in =~ s/[\r\n]//g;
    my @list = split("\t", $in);
    my $id = $list[$index];
    if (defined $spstr) {
      $spindex = 0 if (!defined $spindex);
      my @l = split($spstr, $list[$index]);
      $id = $l[$spindex];
      $id =~ s/\s//g;
    }
    $hash->{$id} = [@list];
    push @$order, $id;
  }
  close($fh);
  return [$hash, $order];
}

sub getXYstats {
  my ($x, $y, $ph, %params) = @_;
  my ($phr, $phw) = @$ph;
  if (!defined $ph) {
    use IPC::Open2;
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  my (@vx, @vy);
  for (my $i = 0; $i < scalar(@$x); $i++) {
    next if (!defined $x->[$i] || $x->[$i] eq "");
    next if (!defined $y->[$i] || $y->[$i] eq "");
    push @vx, $x->[$i];
    push @vy, $y->[$i];
  }
  print $phw "x <- c(", join(",", @vx), ")\n";
  print $phw "y <- c(", join(",", @vy), ")\n";
  print $phw <<END
f.lm <- lm(y ~ x)
st <- summary(f.lm)
f <- st\$fstatistic
p <- pf(f[1], f[2], f[3], lower=FALSE)
pf <- format.pval(p)
cf <- coef(f.lm)
r <- cor(x, y)
n <- length(x)
#sr <- r / sqrt(1 - r * r ) * sqrt(max(0, n - 2))
#pr <-  (1 - pt(abs(sr), max(0, n - 2))) * 2
cat("BEGIN\\n")
cat(sprintf("R^2 = \%.3f\\n", st\$r.squared))
cat(sprintf("n = \%d\\n", n))
cat(sprintf("a0 = \%f\\n", cf[1]))
cat(sprintf("a1 = \%f\\n", cf[2]))
cat(sprintf("r = \%f\\n", r))
cat(sprintf("pvalue = \%s\\n", pf))
cat(sprintf("m1 = \%f\\n", mean(x)))
cat(sprintf("sd1 = \%f\\n", sd(x)))
cat(sprintf("m2 = \%f\\n", mean(y)))
cat(sprintf("sd2 = \%f\\n", sd(y)))
#cat(sprintf("sr = \%f\\n", sr))
#cat(sprintf("pr = \%s\\n", format.pval(pr)))
cat("box1 =", quantile(x, probs=c(0, 0.25, 0.5, 0.75, 1)), "\\n")
cat("box2 =", quantile(y, probs=c(0, 0.25, 0.5, 0.75, 1)), "\\n")
cat("END\\n")
END
;
  my $hash;
  while (<$phr>) {
    last if (/^END/);
    $hash = {} if (/^BEGIN/);
    next if (!defined $hash);
    next if (!/=/);
    s/[\r\n]//g;
    my ($k, $v) = split(" = ", $_, 2);
    $hash->{$k} = $v;
  }
  if (!defined $ph) {
    close($phw);
    close($phr);
  }
  return $hash;
}

sub getXstats {
  my ($x, $ph, %params) = @_;
  my ($phr, $phw) = @$ph;
  if (!defined $ph) {
    use IPC::Open2;
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  my @vx;
  for (my $i = 0; $i < scalar(@$x); $i++) {
    next if (!defined $x->[$i] || $x->[$i] eq "");
    push @vx, $x->[$i];
  }
  print $phw "x <- c(", join(",", @vx), ")\n";
  print $phw <<END
cat("BEGIN\\n")
cat(sprintf("n = \%d\\n", length(x)))
cat(sprintf("m = \%f\\n", mean(x)))
cat(sprintf("sd = \%f\\n", sd(x)))
cat("box =", quantile(x, probs=c(0, 0.25, 0.5, 0.75, 1)), "\\n")
cat("END\\n")
END
;
  my $hash;
  while (<$phr>) {
    last if (/^END/);
    $hash = {} if (/^BEGIN/);
    next if (!defined $hash);
    next if (!/=/);
    s/[\r\n]//g;
    my ($k, $v) = split(" = ", $_, 2);
    $hash->{$k} = $v;
  }
  if (!defined $ph) {
    close($phw);
    close($phr);
  }
  return $hash;
}

sub chisq {
  my @list = @_;
  use IPC::Open2;
  my ($phr, $phw);
  my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  print $phw "p <- chisq.test(matrix(c(", join(",", @list), "), ncol=2))\$p.value\n";
  print $phw "cat(p, \"\\n\")\n";
  my $p = <$phr>;
  $p =~ s/[\r\n]//g;
  close($phw);
  close($phr);
  return $p;
}

sub ciw {
  my ($arrayref, $start, $end) = @_;
  if (!defined $start) {
    $start = 0;
  }
  if (!defined $end) {
    $end = scalar(@$arrayref) - 1;
  }
  my $n = $end - $start + 1;
  my $stddev = &stddev($arrayref, $start, $end);
  use IPC::Open2;
  my ($phr, $phw);
  my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  print $phw "p <- qt(0.975, $n) * $stddev / sqrt($n)\n";
  print $phw "cat(p, \"\\n\")\n";
  my $p = <$phr>;
  $p =~ s/[\r\n]//g;
  close($phw);
  close($phr);
  return $p;
}

sub kmeans {
  my ($expr, $n) = @_;
  use IPC::Open2;
  my ($phr, $phw);
  my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  print $phw "d <- c(", join(",", @$expr), ")\n";
  print $phw "cl <- kmeans(d, $n)\n";
  print $phw "cat(";
  print $phw "max(d[cl\$cluster == $_]), " foreach (1 .. $n);
  print $phw "'\\n', sep = '\\t')\n";
  my $p = <$phr>;
  $p =~ s/[\r\n]//g;
  my @res = split("\t", $p);
  close($phw);
  close($phr);
  return [@res];
}

sub getExprSoft {
  my ($file, $gsmlist, $logr, %params) = @_;
  my $hhash = {};
  my $nogsm = 1;
  for (my $i = 0; $i < scalar(@$gsmlist); $i++) {
    $hhash->{$gsmlist->[$i]->[4]} = $i;
    if (defined $gsmlist->[$i]->[4] && $gsmlist->[$i]->[4] ne "") {
      $nogsm = 0;
    }
  }
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $phash = {};
  my $platform;
  my $platformids = [];
  my $symbols = {};
  my $collect_platform = 0;
  my $collect_headers = 0;
  my @allheaders;
  my $hash = {};
  my $sampleid;
  my $sampleplatformid;
  my $list = [];
  my $valuefield = "VALUE";
  if (defined $params{"valuefield"}) {
    $valuefield = $params{"valuefield"};
  }
  while (<$fh>) {
    s/[\r\n]//g;
    if (/^\^SERIES = (.*)$/) {
      $collect_headers = 1;
    }
    if ($collect_headers == 1 && /^!Series_sample_id = (.*)$/) {
      push @allheaders, $1;
    }
    if (/^\^PLATFORM = (.*)/) {
      $platform = $1;
      $collect_platform = 1;
    }
    if ($collect_platform == 1 && /^!platform_table_begin/) {
      $collect_platform = 0;
      my $head = <$fh>;
      $head =~ s/[\r\n]//g;
      my @h = split("\t", $head);
      my $hh = {};
      for (my $i = 0; $i < scalar(@h); $i++) {
        $hh->{$h[$i]} = $i;
      }
      my ($namefield, $descfield);
      if (defined $params{"namefield"}) {
        $namefield = $params{"namefield"};
      }
      if (defined $params{"descfield"}) {
        $descfield = $params{"descfield"};
      }
      while (<$fh>) {
        if (/^!platform_table_end/) {
          last;
        }
        s/[\r\n]//g;
        my @l = split("\t");
        my $id = $l[0];
        $id =~ s/\s//g;
        push @$platformids, $id;
        my ($name, $desc);
        if (defined $namefield && defined $descfield) {
          $name = $l[$hh->{$namefield}];
          $desc = $l[$hh->{$descfield}];
        }
        elsif (defined $hh->{"Gene Symbol"} && defined $hh->{"Gene Title"}) {
          ($name, $desc) = ($l[$hh->{"Gene Symbol"}], $l[$hh->{"Gene Title"}]);
        }
        elsif (defined $hh->{"GENE SYMBOL"} && defined $hh->{"SEQUENCE DESCRIPTION"}) {
          ($name, $desc) = ($l[$hh->{"GENE SYMBOL"}], $l[$hh->{"SEQUENCE DESCRIPTION"}]);
        }
        elsif (defined $hh->{"GENE_SYMBOL"} && defined $hh->{"GENE_DESCRIPTION"}) {
          ($name, $desc) = ($l[$hh->{"GENE_SYMBOL"}], $l[$hh->{"GENE_DESCRIPTION"}]);
        }
        elsif (defined $hh->{"GENE_SYMBOL"} && defined $hh->{"GENE_NAME"}) {
          ($name, $desc) = ($l[$hh->{"GENE_SYMBOL"}], $l[$hh->{"GENE_NAME"}]);
        }
        elsif (defined $hh->{"SYMBOL"} && defined $hh->{"GENE_NAME"}) {
          ($name, $desc) = ($l[$hh->{"SYMBOL"}], $l[$hh->{"GENE_NAME"}]);
        }
        elsif (defined $hh->{"Gene Symbol"} && defined $hh->{"Gene Name"}) {
          ($name, $desc) = ($l[$hh->{"Gene Symbol"}], $l[$hh->{"Gene Name"}]);
        }
        elsif (defined $hh->{"V3.0.3_gene_symbol"} && defined $hh->{"V3.0.3_control_description"}) {
          ($name, $desc) = ($l[$hh->{"V3.0.3_gene_symbol"}], $l[$hh->{"V3.0.3_control_description"}]);
        }
        elsif (defined $hh->{"GeneName"} && defined $hh->{"Description"}) {
          ($name, $desc) = ($l[$hh->{"GeneName"}], $l[$hh->{"Description"}]);
        }
        elsif (defined $hh->{"Common"} && defined $hh->{"Description"}) {
          ($name, $desc) = ($l[$hh->{"Common"}], $l[$hh->{"Description"}]);
        }
        elsif (defined $hh->{"Description"} && defined $hh->{"GB_ACC"}) {
          $desc = $l[$hh->{"Description"}];
          $name = $l[$hh->{"GB_ACC"}];
        }
        elsif (defined $hh->{"gene_assignment"}) {
          my ($refseq, $n, $d, $rest) = split(/ \/\/ /, $l[$hh->{"gene_assignment"}]);
          ($name, $desc) = ($n, $d);
        }
        elsif (defined $hh->{"gene_symbol"} && defined $hh->{"GB_ACC"}) {
          $name = $l[$hh->{"gene_symbol"}];
          $desc = $l[$hh->{"GB_ACC"}];
        }
        elsif (defined $hh->{"GeneSymbol"} && defined $hh->{"GB_ACC"}) {
          $name = $l[$hh->{"GeneSymbol"}];
          $desc = $l[$hh->{"GB_ACC"}];
        }
        elsif (defined $hh->{"OligoSet_geneSymbol"} && defined $hh->{"OligoSet_description"}) {
          $name = $l[$hh->{"OligoSet_geneSymbol"}];
          $desc = $l[$hh->{"OligoSet_description"}];
        }
        elsif (defined $hh->{"ORF"} && defined $hh->{"description"}) {
          $name = $l[$hh->{"ORF"}];
          $desc = $l[$hh->{"description"}];
        }
        elsif (defined $hh->{"ORF"} && defined $hh->{"Gene Name"}) {
          $name = $l[$hh->{"ORF"}];
          $desc = $l[$hh->{"Gene Name"}];
        }
        elsif (defined $hh->{"ORF"} && defined $hh->{"GB_ACC"}) {
          $name = $l[$hh->{"ORF"}];
          $desc = $l[$hh->{"GB_ACC"}];
        }
        elsif (defined $hh->{"Symbol"} && defined $hh->{"Definition"}) {
          $name = $l[$hh->{"Symbol"}];
          $desc = $l[$hh->{"Definition"}];
        }
        elsif (defined $hh->{"SYMBOL"} && defined $hh->{"DEFINITION"}) {
          $name = $l[$hh->{"SYMBOL"}];
          $desc = $l[$hh->{"DEFINITION"}];
        }
        elsif (defined $hh->{"miRNA_ID_LIST"} && defined $hh->{"SPOT_ID"}) {
          $name = $l[$hh->{"miRNA_ID_LIST"}];
          $desc = $l[$hh->{"SPOT_ID"}];
        }
        elsif (defined $hh->{"GB_LIST"} && defined $hh->{"SPOT_ID"}) {
          $name = $l[$hh->{"GB_LIST"}];
          $desc = $l[$hh->{"SPOT_ID"}];
        }
        elsif (defined $hh->{"GENE_SYMBOL"} && defined $hh->{"GB_ACC"}) {
          ($name, $desc) = ($l[$hh->{"GENE_SYMBOL"}], $l[$hh->{"GB_ACC"}]);
        }
        if (!defined $name || $name eq "") {
          $name = "---";
        }
        $symbols->{$id} = [$name, $desc];
      }
      $phash->{$platform} = [$platformids, $symbols];
      $platformids = [];
      $symbols = {};
    }

    if (/^\^SAMPLE = (.*)$/) {
      my $sid = $1;
      $sid =~ s/\s//g;
      if (defined $sampleid && ($nogsm == 1 || defined $hhash->{$sampleid})) {
        $hash->{$sampleplatformid}->{$sampleid} = $list;
        print STDERR "$sampleid\n";
      }
      $sampleid = $sid;
      $list = [];
    }
    if (/^!Sample_platform_id = (.*)/) {
      $sampleplatformid = $1;
    }
    if (/^!Sample_title = (.*)$/) { $list->[0] = $1; }
    if (/^!Sample_cha.*tics_ch. = (.*)$/) { push @$list, $1; }
    if (/^!Sample_biomaterial_provider_ch1 = (.*)$/) { push @$list, "biom1: $1"; }
    if (/^!Sample_source_name_ch1 = (.*)$/) { push @$list, "src1: $1"; }
    if (/^#$valuefield = (.*)$/) { push @$list, "valuefield: $1"; }
    if (/^!Sample_description = (.*)$/) { push @$list, "desc: $1"; }
    if (/^!Sample_series_id = (.*)$/) { push @$list, "Series: $1"; }
    if (/^!Sample_relation = BioSample: (.*)$/) { push @$list, "BioSample: $1"; }
    if (/^!Sample_relation = SRA: (.*)$/) { push @$list, "SRA: $1"; }

    if (/^!sample_table_begin/) {
      my $head = <$fh>;
      $head =~ s/[\r\n]//g;
      my @h = split("\t", $head);
      my $hh = {};
      for (my $i = 0; $i < scalar(@h); $i++) {
        $hh->{$h[$i]} = $i;
      }
      my $data = {};
      while (<$fh>) {
        if (/^!sample_table_end/) {
          last;
        }
        s/[\r\n]//g;
        my @l = split("\t");
        my $id = $l[0];
        $id =~ s/\s//g;
        $value = $l[$hh->{$valuefield}];
        if ($value  eq "null") {
          $value = "";
        }
        if (defined $logr && $logr eq 1) {
          if ($value > 0) {
            $data->{$id} = log($value)/log(2);
          }
          else {
            $data->{$id} = $params{"logdefault"};
          }
        }
        elsif (defined $logr && $logr eq 2) {
          if ($value > 1) {
            $data->{$id} = log($value)/log(2);
          }
          elsif ($value > -1) {
            $data->{$id} = $value - 1;
          }
          else {
            $data->{$id} = -2 - log(-$value)/log(2);
          }
        }
        else {
          $data->{$id} = $value;
        }
      }
      push @$list, $data;
    }
  }
  close($fh);
  if (defined $sampleid && ($nogsm == 1 || defined $hhash->{$sampleid})) {
    $hash->{$sampleplatformid}->{$sampleid} = $list;
    print STDERR "$sampleid\n";
  }
  if ($nogsm == 1) {
    my $ann = $gsmlist->[0];
    $ann = [] if (!defined $ann);
    $gsmlist = [];
    foreach (@allheaders) {
      $ann->[4] = $_;
      push @$gsmlist, [map { $_ } @$ann];
    }
    $hhash = {};
    for (my $i = 0; $i < scalar(@$gsmlist); $i++) {
      $hhash->{$gsmlist->[$i]->[4]} = $i;
    }
  }

  return [$hash, $phash, $gsmlist];
}

sub getInfoSoft {
  my ($file, $rest) = @_;
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $hash = {};
  my $sid;
  while (<$fh>) {
    s/[\r\n]//g;
    if (/^\^SAMPLE = (.*)$/) {
      $sid = $1;
      $hash->{$sid} = {};
      print STDERR "$sid\n";
    }
    if (defined $sid) {
      if (/^(.*) = (.*)/) { push @{$hash->{$sid}->{$1}}, $2; }
      if (/^!sample_table_begin/) {
        $sid = undef;
      }
    }
  }
  close($fh);
  return $hash;
}

sub writeExprSoft {
  my ($prefix, $hash, $phash, $gsmlist) = @_;
  foreach my $p (keys %{$hash}) {
    my @headers;
    my $hhead = {};
    for (my $i = 0; $i < scalar(@$gsmlist); $i++) {
      if (defined $hash->{$p}->{$gsmlist->[$i]->[4]} && !defined $hhead->{$gsmlist->[$i]->[4]}) {
        push @headers, $gsmlist->[$i]->[4];
        $hhead->{$gsmlist->[$i]->[4]} = $i;
      }
    }
# collecting array header and writing clinical headers
    my $ofh;
    open($ofh, ">$prefix\-$p\-ih.txt") || die "Can't open $prefix\-$p\-ih.txt\n";
    print $ofh "ArrayID\tArrayHeader\tClinicalhHeader\n";
    for (my $i= 0; $i < scalar(@headers); $i++) {
      my $arr = $headers[$i];
      my $title = $hash->{$p}->{$arr}->[0];
      print $ofh "$arr\t$arr\t$title\n";
    }
    close($ofh);
# writing annotation data
    open($ofh, ">$prefix\-$p\-ann.txt") || die "Can't open $prefix\-$p\-ann.txt\n";
    print $ofh "ArrayID\tannlist\n";
    foreach my $k (keys(%{$hash->{$p}})) {
      my @list = @{$hash->{$p}->{$k}};
      print $ofh $k, "\t", join("\t", @{$gsmlist->[$hhead->{$k}]}), "\t", join("\t", @list[0 .. ($#list - 1)]), "\n";
    }
    close($ofh);

# writing expr and idx data
    open($ofh, ">$prefix\-$p\-expr.txt") || die "Can't open $prefix\-$p\-expr.txt\n";
    print $ofh "ProbeID\tName";
    for (my $i= 0; $i < scalar(@headers); $i++) {
      my $arr = $headers[$i];
      print $ofh "\t$arr";
    }
    print $ofh "\n";
    my $oifh;
    open($oifh, ">$prefix\-$p\-idx.txt") || die "Can't open $prefix\-$p\-idx.txt\n";
    print $oifh "ProbeID\tPtr\tName\tDescription\n";
    foreach my $id (@{$phash->{$p}->[0]}) {
      my $desc = $phash->{$p}->[1]->{$id}->[1];
      my $name = $phash->{$p}->[1]->{$id}->[0];
      my $ptr = tell($ofh);
      print $oifh "$id\t$ptr\t$name\t$desc\n";
      print $ofh "$id\t$name: $desc";
      for (my $i= 0; $i < scalar(@headers); $i++) {
        my $arr = $headers[$i];
        my $l = $hash->{$p}->{$arr};
        my $ex = $l->[scalar(@$l) - 1]->{$id};
        print $ofh "\t$ex";
      }
      print $ofh "\n";
    }
    close($oifh);
    close($ofh);
  }
}

#ID Name Ptr Desc
sub buildIndex {
  my $pclFile = shift;
  my $hash = {};
  open(FL, "< $pclFile") || die "Can't open $pclFile\n";
  my $index = 0;
  my $ptr = tell(FL);
  while (<FL>) {
    my @list = split("\t", $_);
    my ($name, $desc) = split(":", $list[1], 2);
    if ($name =~ /^\s*$/) {
       $name = "---";
    }
    $hash->{$list[0]} = [$list[0], $name, $ptr, $desc];
    $index++;
    $ptr = tell(FL);
  }
  close(FL);
  return $hash;
}

sub writeIhSurvival {
  my ($obj, $if, $cf, %params) = @_;
  my $hash = $obj->{'hash'};
  my $headers = $obj->{'headers'};
  my $cheaders = $obj->{'cheaders'};
  my $narray = $obj->{'narray'};
  my $carray = $obj->{'carray'};

  my $tindex = 0;
  $tindex = $params{'tindex'} if (defined $params{'tindex'});

# collecting array header and writing clinical headers
  my $ofh;
  open($ofh, ">$if") || die "Can't open $if\n";
  print $ofh "ArrayID\tArrayHeader\tClinicalhHeader\n";
  for (my $i= 0; $i < scalar(@$headers); $i++) {
    my $arr = $headers->[$i];
    next if (!defined $hash->{$arr});
    my $title = $hash->{$arr}->[2]->[$tindex];
    print $ofh "$arr\t$arr\t$title\n";
  }
  close($ofh);

# writing survival data
  open($ofh, ">$cf") || die "Can't open $cf\n";
  print $ofh "ArrayID\ttime\tstatus";
  foreach (@$narray) {
    print $ofh "\t", "n ", $cheaders->[$_];
  }
  foreach (@$carray) {
    print $ofh "\t", "c ", $cheaders->[$_];
  }
  print $ofh "\n";
  foreach my $k (keys(%{$hash})) {
    print $ofh $k ."\t".$hash->{$k}->[0]."\t".$hash->{$k}->[1];
    foreach (@$narray, @$carray) {
      print $ofh "\t", $hash->{$k}->[2]->[$_];
    }
    print $ofh "\n";
  }
  close($ofh);

}

sub writeExprIdx {
  my ($obj, $ef, $idxf, %params) = @_;
  my $hash = $obj->{'hash'};
  my $headers = $obj->{'headers'};
  my $platformids = $obj->{'platformids'};
  my $symbols = $obj->{'symbols'};

# writing expr and idx data
  open($ofh, ">$ef") || die "Can't open $ef\n";
  print $ofh "ProbeID\tName";
  for (my $i= 0; $i < scalar(@$headers); $i++) {
    my $arr = $headers->[$i];
    next if (!defined $hash->{$arr});
    print $ofh "\t$arr";
  }
  print $ofh "\n";
  my $oifh;
  open($oifh, ">$idxf") || die "Can't open $idxf\n";
  print $oifh "ProbeID\tPtr\tName\tDescription\n";
  foreach my $id (@$platformids) {
    my $desc = $symbols->{$id}->[1];
    my $name = $symbols->{$id}->[0];
    my $ptr = tell($ofh);
    print $oifh "$id\t$ptr\t$name\t$desc\n";
    print $ofh "$id\t$name: $desc";
    for (my $i= 0; $i < scalar(@$headers); $i++) {
      my $arr = $headers->[$i];
      next if (!defined $hash->{$arr});
      my $l = $hash->{$arr}->[2];
      my $ex = $l->[scalar(@$l) - 1]->{$id};
      print $ofh "\t$ex";
    }
    print $ofh "\n";
  }
  close($oifh);
  close($ofh);
}

sub writeObject {
  my ($obj, $if, $cf, $ef, $idxf, %params) = @_;
  &writeIhSurvival($obj, $if, $cf, %params);
  &writeExprIdx($obj, $ef, $idxf, %params);
}

sub getRandPerm {
  my $min = shift;
  my $max = shift;
  my @numbers = $min .. $max;
  my @rnumbers = map { int(rand($max-$min+1)) } @numbers;
  my @perm = sort { $rnumbers[$a] <=> $rnumbers[$b] } @numbers;
  #return (0 .. ($min-1) , @perm);
  return @perm;
} 

sub getCorrPlot {
  my ($ofile, $v1, $v2, %params) = @_;
  my $width = 650;
  my $height = 600;
  my $dpi = 100;
  $width = $params{'width'} if (defined $params{'width'});
  $height = $params{'height'} if (defined $params{'height'});
  $dpi = $params{'dpi'} if (defined $params{'dpi'});
  my $wd = $width/$dpi;
  my $ht = $height/$dpi;
  $wd = $params{'wd'} if (defined $params{'wd'});
  $ht = $params{'ht'} if (defined $params{'ht'});
  my $lx = 6;
  my $ly = 12;
  $lx = $params{'lx'} if (defined $params{'lx'});
  $ly = $params{'ly'} if (defined $params{'ly'});
  my ($phr, $phw);
  if (defined $params{'ph'}) {
    ($phr, $phw) = @{$params{'ph'}};
  }
  else {
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  my (@vv1, @vv2);
  for (my $i = 0; $i < scalar(@$v1); $i++) {
    next if (!defined $v1->[$i] || $v1->[$i] eq "");
    next if (!defined $v2->[$i] || $v2->[$i] eq "");
    push @vv1, $v1->[$i];
    push @vv2, $v2->[$i];
  }
  my $ivv1 = &info(\@vv1);
  my $ivv2 = &info(\@vv2);
  if (!defined $params{'lx'}) {
    $lx = $ivv1->[0] + 0.7 * ($ivv1->[2] - $ivv1->[0]);
  }
  if (!defined $params{'ly'}) {
    $ly = $ivv2->[0] + 0.2 * ($ivv2->[2] - $ivv1->[0]);
  }
  print $phw "x <- c(", join(",", @vv1), ");\n";
  print $phw "y <- c(", join(",", @vv2), ");\n";
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
  my ($nrpoints, $xt, $yt) = (500, "X", "Y");
  $xt = $params{'xt'} if (defined $params{'xt'});
  $yt = $params{'yt'} if (defined $params{'yt'});
  $nrpoints = $params{'nrpoints'} if (defined $params{'nrpoints'});
  if (1) {
  print $phw <<END
$devstr
f.lm <- lm(y ~ x)
plot(x, y, pch=3, ylab="$yt", xlab="$xt")
abline(f.lm)
st <- summary(f.lm)
f <- st\$fstatistic
p <- pf(f[1], f[2], f[3], lower=FALSE)
pf <- format.pval(p)
cf <- coef(f.lm)
text($lx, $ly, paste(sprintf("R^2 = \%.3f", st\$r.squared), "\n", 
    sprintf("r = \%f", cor(x, y)), "\n",
    sprintf("n = \%d", length(x)), "\n",
    sprintf("y = \%.3f + \%.3f x", cf[1], cf[2]), "\n",
    "p-value: ", pf, sep=""), pos = 4, col="blue")
dev.off()
END
;
  }
  if (0) {
  print $phw <<END
$devstr
require(geneplotter)
f.lm <- lm(y ~ x)
smoothScatter(x, y, nrpoints = $nrpoints, pch=16, ylab="$yt", xlab="$xt")
abline(f.lm)
st <- summary(f.lm)
f <- st\$fstatistic
p <- pf(f[1], f[2], f[3], lower=FALSE)
pf <- format.pval(p)
cf <- coef(f.lm)
text($lx, $ly, paste(sprintf("R^2 = \%.3f", st\$r.squared), "\n", 
    sprintf("r = \%f", cor(x, y)), "\n",
    sprintf("n = \%d", length(x)), "\n",
    sprintf("y = \%.3f + \%.3f x", cf[1], cf[2]), "\n",
    "p-value: ", pf, sep=""), pos = 4, col="blue")
dev.off()
END
;
  }

}

#Source: http://josiahmanson.com/prose/student_t/
sub student_pdf {
  use POSIX;
  my ($t, $n) = @_;
  return (0.3989422804014327*$n*pow(($n + $t*$t)/$n,-0.5 - 0.5*$n))
    / (0.2551219234157433 + $n);
}
sub student_prob_different {
  use POSIX;
  my ($t, $n) = @_;
  if ($t <= 0) {
    return 0;
  }
  my $steps = int(ceil($t * 100));
  my $dt = $t / $steps;
  my $sum = 0; 
  for (my $x = $dt*.5; $x < $t; $x += $dt) {
    $sum += &student_pdf($x, $n);
  }
  $sum *= 2 * $dt;
  return $sum;
}
sub student_prob_same {
  my ($t, $n) = @_;
  return (1 - &student_prob_different($t, $n));
}
sub ttest {
  my ($a, $b) = @_;
  my $ma = &mean($a);
  my $mb = &mean($b);
  my $sa = &stddev($a);
  my $sb = &stddev($b);
  my $na = scalar(@$a);
  my $nb = scalar(@$b);
  if ($na < 3 || $nb < 3) {
    return 1;
  }
  if ($sa <= 0 && $sb <= 0) {
    if ($ma != $mb) { return 0; }
    else { return 1; }
  }
  my $t = abs($ma - $mb) / sqrt($sa * $sa / $na + $sb*$sb / $nb);
  my $dd = ($sa * $sa * $sa * $sa / $na / $na / ($na - 1) +
        $sb * $sb * $sb * $sb / $nb / $nb / ($nb - 1));
  if ($dd <= 0) {
    return 1;
  }     
  my $df = ($sa * $sa / $na + $sb*$sb / $nb) * ($sa * $sa / $na + $sb*$sb / $nb) / $dd;
  return &student_prob_same($t, $df);
}

sub scatterPlot {
  my ($outfile, $expr1, $expr2, %params) = @_;
  my $mew = 1.1;
  my $ms = 4;
  $mew = $params{"mew"} if (defined $params{"mew"});
  $ms = $params{"ms"} if (defined $params{"ms"});
  my $xn = "X";
  my $yn = "Y";
  $xn = $params{"xn"} if (defined $params{"xn"});
  $yn = $params{"yn"} if (defined $params{"yn"});
  my ($phr, $phw);
  my $pid = open2($phr, $phw, "python");
  print $phw <<END
import matplotlib
matplotlib.use('agg')
import re
from pylab import *
from numpy import *
fig = figure(figsize=(6.4,4.8))
ax = fig.add_axes([70.0/640, 54.0/480, 1-2*70.0/640, 1-2*54.0/480])
END
;
    my ($datax, $datay) = ([], []);
    for(my $i = 0; $i < scalar(@{$expr1}); $i++) {
      next if (!defined $expr1->[$i] || $expr1->[$i] eq "");
      next if (!defined $expr2->[$i] || $expr2->[$i] eq "");
      push @$datax, $expr1->[$i];
      push @$datay, $expr2->[$i];
    }
    print $phw "datax = [", join(",", @$datax), "]\n";
    print $phw "datay = [", join(",", @$datay), "]\n";
    print $phw "c = 'blue'\n";
    print $phw "ax.plot(datax,datay, color=c, ls='None', marker='+', mew=$mew, ms=$ms, mec=c)\n";
    print $phw "ax.axis([min(datax)-0.5, max(datax)+0.5, min(datay)-0.5, max(datay)+0.5])\n";
  print $phw <<END
ax.set_xlabel('$xn', fontsize=10)
ax.set_ylabel('$yn', fontsize=10)
fig.savefig('$outfile', dpi=100)
END
;
  close($phr);
  close($phw);
  sleep(1);
  chmod 0666, $outfile;
}

sub BinarySearch {
   my ($A, $value) = @_;
   my $low = 0;
   my $N = scalar(@$A);
   my $high = $N - 1;
   my $mid = 0;
   while ($low <= $high) {
       $mid = int($low + (($high - $low) / 2));
       if ($A->[$mid] > $value) {
           $high = $mid - 1;
       }
       elsif ($A->[$mid] < $value) {
           $low = $mid + 1;
       }
       else {
           return $mid;
       }
   }
   return $high;
}

sub getRangeInfo {
  my ($expr, $thr) = @_;
  my $sum = 0;
  my $sum2 = 0;
  my $num = 0;
  my $thrNum = 0;
  my $min = $expr->[0]; my $max = $expr->[0];
  foreach (@$expr) {
    next if (/^\s*$/);
    $sum += $_;
    $sum2 += $_*$_;
    if ($min > $_) {
      $min = $_;
    }
    if ($max < $_) {
      $max = $_;
    }
    if ($_ < $thr) {
      $thrNum ++;
    }
    $num ++;
  }
  if ($num > 0) {
    my $mean = $sum / $num;
    my $mean2 = $sum2 / $num;
    my $sd = sqrt(abs($mean2 - $mean*$mean));
    return ($min, $max, $mean, $sd, $num, $thrNum);
  }
  return undef;
}

sub kappa {
  my ($a, $b) = @_;
  my $all = &union($a, $b); 
  my $hash = {};
  my $index = 0;
  foreach (@$all) {
    if (!defined $hash->{$_}) {
      $hash->{$_} = $index;
      $index++;
    }
  }
  my $m = [];
  my $len = &min([scalar(@$a), scalar(@$b)]);
  foreach my $i (0 .. ($index - 1)) {
    foreach my $j (0 .. ($index - 1)) {
      $m->[$i]->[$j] = 0;
    }
  }
  foreach my $i (0 .. ($len - 1)) {
    my $k1 = $hash->{$a->[$i]};
    my $k2 = $hash->{$b->[$i]};
    if (!defined $k1 || !defined $k2) {
      next;
    }
    $m->[$k1]->[$k2] ++;
  }
  return &computeKappa($m);
}

sub computeKappa {
  my ($m, $rest) = @_;
  my $n = 0;
  foreach (@$m) {
    foreach (@$_) {
      $n += $_;
    }
  }
  $k = $#{ $m->[0] };

  my $pa = 0;
  for my $j (0 .. $k){
    $pa += $m->[$j]->[$j]/ $n;
  }
  my $p1 = [];
  for my $j (0 .. $k){
    $p1->[$j] = 0;
    for my $i (0 .. $k){
      $p1->[$j] += $m->[$i]->[$j];
    }
  }
  my $p2 = [];
  for my $i (0 .. $k){
    $p2->[$i] = 0;
    for my $j (0 .. $k){
      $p2->[$i] += $m->[$i]->[$j];
    }
  }
  my $pe = 0;
  for my $j (0 .. $k){
    $pe += $p1->[$j] * $p2->[$j] / $n / $n;
  }

  $kappa = ($pa - $pe)/(1 - $pe);
  return $kappa;
}
 
sub getInfo {
  my ($file, $hdr, $idlist, $rest) = @_;
  open($fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @h = split("\t", $head);
  my $hh = {};
  for (my $i = 0; $i < scalar(@h); $i++) {
    $hh->{$h[$i]} = $i;
  }

  return undef if (!defined $hh->{$hdr});

  my $index = 0;
  my $hash = {};
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    $hash->{$list[0]} = $list[$hh->{$hdr}];
    $index++;
  }
  close($fh);

  return [ map { $hash->{$_} } @$idlist ];
}

sub writeData {
  my ($tmpfile, $expr, $rows, $columns) = @_;
  open($ofh, ">$tmpfile") || die "Can't open $tmpfile\n";
  print $ofh join("\t", @$columns), "\n";
  my $hash = {};
  for (my $i = 0; $i <= $#{$rows}; $i++) {
    my $id = $rows->[$i];
    if (!defined $id || $id eq "") { $id = "---"; }
    if (defined $hash->{$id}) {
      $hash->{$id} ++;
      $id = $id.".".$hash->{$id};
    }
    else {
      $hash->{$id} = 0;
    }
    print $ofh join("\t", $id, map { $expr->[$i]->[$_] } 0 .. $#{$columns}), "\n";
  }
  close($ofh);
}

sub plotHistogram {
  my ($file, $par) = @_;
  my %params = %{$par};
  use Statistics::R;
  use HR;
  my $R = HR->new();
  my $R = Statistics::R->new();

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
  my $ofile = "_hist.png";
  if (defined $params{"out"}) {
    $ofile = $params{"out"};
  }
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
  my ($breaks, $xt, $yt) = (10, "value", "count");
  $xt = $params{'xt'} if (defined $params{'xt'});
  $yt = $params{'yt'} if (defined $params{'yt'});
  $breaks = $params{'breaks'} if (defined $params{'breaks'});
  my $main = "Histogram";
  $main = $params{'main'} if (defined $params{'main'});
  if ($params{'logy'}) {
  $R->run(qq{
$devstr
data <- read.table('$file', sep='\t')
l <- data[,1]
h <- hist(l, breaks=$breaks, freq=TRUE, main="$main", xlab="$xt", ylab="$yt")
n <- length(h\$counts)
plot(h\$breaks[1:n], log10(h\$counts), main="$main", xlab="$xt", ylab="$yt")
dev.off()
});
  }
  else {
  $R->run(qq{
$devstr
data <- read.table('$file', sep='\t')
l <- data[,1]
h <- hist(l, breaks=$breaks, freq=TRUE, main="$main", xlab="$xt", ylab="$yt")
n <- length(h\$counts)
plot(h\$breaks[1:n], h\$counts, main="$main", xlab="$xt", ylab="$yt")
dev.off()
});
  }
  my $br = $R->get("h\$breaks");
  my $counts = $R->get("h\$counts");
  my $min = $R->get("min(l)");
  my $m = $R->get("mean(l)");
  my $max = $R->get("max(l)");
  my $sd = $R->get("sd(l)");
  my $sum = $R->get("sum(l)");
  my $info = [$min, $m, $max, $sd, $sum];
  my @breaks = split(" ", join(" ", @$br));
  my @counts = split(" ", join(" ", @$counts));
  return [\@breaks, \@counts, $info];
}

#Replace a string without using RegExp.
sub str_replace { 
  my ($search, $replace, $subject) = @_; 
  my $pos = index($subject, $search); 
  while ($pos > -1) { 
    substr($subject, $pos, length($search), $replace); 
    $pos = index($subject, $search, $pos + length($replace)); 
  } 
  return $subject; 
}

sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

# 
# Series of substitutions to sanitise text for use in LaTeX.
# 
# http://stackoverflow.com/questions/2627135/how-do-i-sanitize-latex-input
# Target document should \usepackage{textcomp}
# 
sub myescape {
  $text = shift;

  # Prepare backslash/newline handling
  $text = str_replace("\n", "\\\\", $text); # Rescue newlines
  $text =~ s/[\x00-\x1F\x7F-\xFF]//g; # Strip all non-printables
  $text = str_replace("\\\\", "\n", $text); # Re-insert newlines and clear \\
  $text = str_replace("\\", "\\\\", $text); # Use double-backslash to signal a backslash in the input (escaped in the final step).

  # Symbols which are used in LaTeX syntax
  $text = str_replace("{", "\\{", $text);
  $text = str_replace("}", "\\}", $text);
  $text = str_replace('$', '\\$', $text);
  $text = str_replace("&", "\\&", $text);
  $text = str_replace("#", "\\#", $text);
  $text = str_replace("^", "\\textasciicircum{}", $text);
  $text = str_replace("_", "\\_", $text);
  $text = str_replace("~", "\\textasciitilde{}", $text);
  $text = str_replace("%", "\\%", $text);

  # Brackets & pipes
  $text = str_replace("<", "\\textless{}", $text);
  $text = str_replace(">", "\\textgreater{}", $text);
  $text = str_replace("|", "\\textbar{}", $text);

  # Quotes
  $text = str_replace("\"", "\\textquotedbl{}", $text);
  $text = str_replace("'", "\\textquotesingle{}", $text);
  $text = str_replace("`", "\\textasciigrave{}", $text);

  # Clean up backslashes from before
  $text = str_replace("\\\\", "\\textbackslash{}", $text); # Substitute backslashes from first step.
  $text = str_replace("\n", "\\\\", trim($text)); # Replace newlines (trim is in case of leading \\)
  return $text;
}

sub getPScolor {
  my ($n, $c) = @_;
  if ($c =~ /^#(..)(..)(..)/) {
    my ($r, $g, $b) = ($1, $2, $3);
    return "\\definecolor{".$n."}{HTML}{"."$r$g$b"."}";
  }
  return "\\definecolor{".$n."}{named}{".$c."}";
}

sub a2n {
  my $n = 0;
  $n = $n * 26 + $_ for map { ord($_) & 0x1F } split //, shift;
  $n;
}

sub n2a {
  my ($n) = @_;
  my @a;
  while ($n > 0) {
    unshift @a, ($n-1) % 26;
    $n = int(($n-1)/26); #Edited for failing corner case
  }
  join '', map chr(ord('a') + $_), @a;
}

sub box_plot_values {
  my $data = shift;
  my $res = {
    'lower_outlier'  => 0,
    'min'            => 0,
    'q1'             => 0,
    'median'         => 0,
    'q3'             => 0,
    'max'            => 0,
    'higher_outlier' => 0,
    'mean'           => 0,
    'stddev'         => 0,
    'size'           => 0,
  };

  use POSIX;

  my $array_count = scalar(@$data);
  my $array = [sort { $a <=> $b } @$data];

  $res->{'size'}           = $array_count;
  $res->{'min'}            = $array->[0];
  $res->{'lower_outlier'}  = [];
  $res->{'max'}            = $array->[$array_count - 1];
  $res->{'higher_outlier'} = [];
  my $middle_index            = POSIX::floor($array_count / 2);
  $res->{'median'}         = $array->[$middle_index];
  my $lower_values            = [];
  my $higher_values           = [];
  $res->{'mean'}           = &mean($array);
  $res->{'stddev'}         = &stddev($array);

  if ($array_count % 2 == 0) {
    $res->{'median'} = 
    (($res->{'median'} + $array->[$middle_index - 1]) / 2);

    foreach my $idx (0 .. $#{$array}) {
      if ($idx <= ($middle_index - 1)) {
        push @$lower_values, $array->[$idx];
      }
      elsif ($idx >= $middle_index) {
        push @$higher_values, $array->[$idx];
      }
    }
  }
  else {
    foreach my $idx (0 .. $#{$array}) {
      if ($idx < $middle_index) {
        push @$lower_values, $array->[$idx];
      }
      elsif ($idx > $middle_index) {
        push @$higher_values, $array->[$idx];
      }
    }
  }

  my $lower_values_count = scalar(@$lower_values);
  my $lower_middle_index = POSIX::floor($lower_values_count / 2);
  $res->{'q1'} = $lower_values->[$lower_middle_index];
  if ($lower_values_count % 2 == 0) {
    $res->{'q1'} = 
    (($res->{'q1'} + $lower_values->[$lower_middle_index - 1]) / 2);
  }

  my $higher_values_count = scalar(@$higher_values);
  my $higher_middle_index = POSIX::floor($higher_values_count / 2);
  $res->{'q3'} = $higher_values->[$higher_middle_index];
  if ($higher_values_count % 2 == 0) {
    $res->{'q3'} = 
    (($res->{'q3'} + $higher_values->[$higher_middle_index - 1]) / 2);
  }

  $iqr = $res->{'q3'} - $res->{'q1'};
  $res->{'iqr'} = $iqr;
  $res->{'min'} = $res->{'q1'} - 1.5*$iqr;
  foreach my $value (@$lower_values) {
    if( $value < $res->{'min'} ) {
      push @{$res->{'lower_outlier'}}, $value ;
    }
    else {
      $res->{'min'} = $value;
      last;
    }
  }

  $res->{'max'} = $res->{'q3'} + 1.5*$iqr;
  foreach my $i (0 .. $#{$higher_values}) {
    my $value = $higher_values->[$#{$higher_values} - $i];
    if( $value > $res->{'max'} ) {
      push @{$res->{'higher_outlier'}}, $value ;
    }
    else {
      $res->{'max'} = $value;
      last;
    }
  }
  return $res;
}

sub getBoxStats {
  my ($pG, $ph) = @_;
  my ($phr, $phw);
  if (defined $ph) {
    ($phr, $phw) = @$ph;
  }
  else {
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  my $num = scalar(@$pG);
  my $res = "cat(\"BEGIN\\n\")\n";
  $res .= "x <- list()\n";
  my $gxa = [];
  my $groups = [];
  my ($i, $j);
  for($i=1; $i <= $num; $i++){
    my $xa = $pG->[$i-1]->[2];
    if (scalar(@$xa) <= 0) {
    $res .= "
cat(paste('info', $i, 0, 0, 0, 0, 0, '\\n', sep='\\t')) 
";
    next;
    }
    push @$gxa, @$xa;
    my $rr = [map { $i} @$xa];
    push @$groups, @$rr;
    $res .= "x[[$i]] <- c(" . join(",", @$xa) . ");\n";
    $res .= "
ex <- x[[$i]]
if (sd(ex) <= 0) {
cat(paste('info', $i, 1, ex[1], sd(ex), ex[1], ex[1], '\\n', sep='\\t')) 
} else {
st <- t.test(ex)
cat(paste('info', $i, length(ex), st\$estimate, sd(ex), 
st\$conf.int[1], st\$conf.int[2], '\\n', sep='\\t')) 
}
";
  }
  for($i=1; $i <= $num; $i++){
    my $xa = $pG->[$i-1]->[2];
    if (scalar(@$xa) <= 0) {
      $res .= "ex <- c()\n";
    }
    else {
      $res .= "ex <- x[[$i]]\n";
    }
    for($j= ($i + 1); $j <= $num; $j++){
    my $xb = $pG->[$j-1]->[2];
    if (scalar(@$xb) <= 0) {
      $res .= "ey <- c()\n";
    }
    else {
      $res .= "ey <- x[[$j]]\n";
    }
      $res .= "
if (length(ex) <= 0 || length(ey) <= 0) {
cat(paste('pvalue', $i, $j, Inf, 1, '\\n', sep='\\t'))
} else if (sd(ey) <= 0 && sd(ex) <= 0) {
if (ex[1] != ey[1]) {
cat(paste('pvalue', $i, $j, Inf, 0, '\\n', sep='\\t'))
} else {
cat(paste('pvalue', $i, $j, Inf, 1, '\\n', sep='\\t'))
}
} else {
st <- t.test(ex, ey)
cat(paste('pvalue', $i, $j, st\$statistic, st\$p.value, '\\n', sep='\\t'))
}
";
    }
  }
  $res .= "values <- c(" . join(",", @$gxa) . ");\n";
  $res .= "groups <- c(" . join(",", @$groups) . ");\n";
  $res .= "
s <- anova(lm(values ~ groups))
cat(paste('anova', s\$'Pr(>F)'[1], s\$'F value'[1], '\\n', sep='\\t'))
";
  $res .= "cat('END\\n')\n";
  print $phw $res;

  my $ihash = [];
  my $phash = [];
  my $ahash = [];
  while (<$phr>) {
    my $line = $_;
    $line =~ s/[\r\n]//g;
    my @list = split("\t", $line);
    if (scalar(@list) > 0 && $list[0] eq "END") {
      last;
    }
    if (scalar(@list) > 0 && $list[0] eq "info") {
      push @$ihash, [@list];
    }
    if (scalar(@list) > 0 && $list[0] eq "pvalue") {
      push @$phash, [@list];
    }
    if (scalar(@list) > 0 && $list[0] eq "anova") {
      push @$ahash, [@list];
    }
  }
  return [$ihash, $phash, $ahash];
}

sub writeTikzBoxplot {
  my ($file, $x_name, $pG, %params) = @_;
  my $bplots = [];
  foreach my $g (@$pG) {
    my $res = &U::box_plot_values($g->[2]);
    push @$bplots, $res;
  }
  my $binfo = &getBoxStats($pG);
  my $xldata = [ 0 .. $#{$pG} ];
  if (defined $params{"sort"}) {
    $xldata = [ sort { 
      scalar(@{$pG->[$a]->[2]}) <=> scalar(@{$pG->[$b]->[2]}) 
      } 0 .. $#{$pG} ];
  }
  open(my $ofh, ">$file") || die "Can't open $file\n";
  print $ofh "\\begin{tikzpicture}\n";
  my $index = 1;
  foreach my $i (@$xldata) {
    my $clr = &U::getPScolor("clr$i", $pG->[$i]->[1]);
    print $ofh "$clr\n";
    if (scalar(@{$pG->[$i]->[2]}) > 0) {
      my $c = &U::n2a($i);
      print $ofh "\\pgfplotstableread{%\n";
      foreach my $v (@{$pG->[$i]->[2]}) {
        print $ofh "$index $v\n";
      }
      print $ofh "}\\g$c"."data%\n";
    }
    $index++;
  }

  my $yl = &U::myescape("$x_name");
  my $res = "
\\begin{axis}[
name=plot1, at={(0,0)},
boxplot/draw direction=y,
x axis line style={opacity=0},
axis x line*=bottom,
axis y line=left,
enlarge y limits,
ylabel=$yl,
ymajorgrids,
";
  print $ofh $res;
  $res = join(",", 1 .. scalar(@$xldata));
  print $ofh "xtick={".$res."},\n";
  $res = join(",", map { $pG->[$_]->[0] } @$xldata);
  print $ofh "xticklabels={".$res."}]\n";
  $index = 1;
  foreach my $i (@$xldata) {
    if (defined $params{"points"} && scalar(@{$pG->[$i]->[2]}) > 0) {
      my $c = &U::n2a($i);
      my $dname = "\\g".$c."data";
      print $ofh "\\addplot+[line width=1.3pt, color=clr$i, mark=*, only marks,mark options={color=clr$i}]  table [x index=0, y index=1] {$dname};\n";
    }
    my $b = $bplots->[$i];
    my $out = [@{$b->{"lower_outlier"}}, @{$b->{"higher_outlier"}}];
    my ($a1, $a2, $a3, $a4, $a5) = ($b->{"min"}, $b->{"q1"}, 
      $b->{"median"}, $b->{"q3"}, $b->{"max"});
    $res = "
    \\addplot+ [ color=clr$i,fill=clr$i!20, mark=o,
      boxplot prepared={
        lower whisker=$a1, lower quartile=$a2,
        median=$a3,
        upper quartile=$a4, upper whisker=$a5,
      },
    ] ";
    print $ofh $res;
    if (scalar(@$out) > 0) {
      $res = join("\\\\ ", @$out) . "\\\\ ";
      print $ofh "table [row sep=\\\\,y index=0] { $res};\n";
    }
    else {
      print $ofh "coordinates {};\n";
    }
    my ($a1, $a2, $a3) = ($b->{"mean"}, $b->{"stddev"}, $b->{"size"});
    my $f = 1.96 * $a2 / sqrt($a3);
    my ($a3, $a4) = ($a1 - $f, $a1 + $f);
    $res = "
    \\node at ($index,$a1) [color=clr$i!50,circle,draw]{};
    \\draw[color=clr$i!50] [<-] ($index,$a3) -- ($index,$a1);
    \\draw[color=clr$i!50] [->] ($index,$a1) -- ($index,$a4);
";
    print $ofh $res;
    $index++;
  }
  print $ofh "\\end{axis}\n";
  $res = "
\\matrix[anchor=north west, xshift=1cm, name=iH] at (plot1.north east)
[matrix of nodes,row sep=0em, column sep=0em] {
info & Group & n & mean & sd & 95\\% CI & \\\\
";
  print $ofh $res;
  foreach my $a (@{$binfo->[0]}) {
    $a->[1] = $pG->[$a->[1] - 1]->[0];
    foreach my $i (3 .. 6) {
      $a->[$i] = sprintf("%.3g", $a->[$i]);
    }
    $res = join(" & ", @$a)." \\\\\n";
    print $ofh $res;
  }
  print $ofh "};\n";
  $res = "
\\matrix[anchor=north west, name=iP] at (iH.south west)
[matrix of nodes,row sep=0em, column sep=0em] {
pvalue & Group 1 & Group 2 & statistic & pvalue \\\\
";
  print $ofh $res;
  foreach my $a (@{$binfo->[1]}) {
    $a->[1] = $pG->[$a->[1] - 1]->[0];
    $a->[2] = $pG->[$a->[2] - 1]->[0];
    $a->[3] = sprintf("%.3g", $a->[3]);
    $a->[4] = sprintf("%.3g", $a->[4]);
    $res = join(" & ", @$a)." \\\\\n";
    print $ofh $res;
  }
  print $ofh "};\n";
  $res = "
\\matrix[anchor=north west, name=iA] at (iP.south west)
[matrix of nodes,row sep=0em, column sep=0em] {
anova & Pr(\$>\$F) & F Value \\\\
";
  print $ofh $res;
  foreach my $a (@{$binfo->[2]}) {
    $a->[1] = sprintf("%.3g", $a->[1]);
    $a->[2] = sprintf("%.3g", $a->[2]);
    $res = join(" & ", @$a)." \\\\\n";
    print $ofh $res;
  }
  print $ofh "};\n";
  print $ofh "\\end{tikzpicture}\n";

}

sub getVal {
  my ($value, $logr, %params) = @_;
  if ($value  eq "null") {
    $value = "";
  }
  if (!defined $params{"logdefault"}) {
    $params{"logdefault"} = "";
  }
  if (defined $logr && $logr eq 1) {
    if ($value > 0) {
      return log($value)/log(2);
    }
    else {
      return $params{"logdefault"};
    }
  }
  elsif (defined $logr && $logr eq 2) {
    if ($value > 1) {
      return log($value)/log(2);
    }
    elsif ($value > -1) {
      return $value - 1;
    }
    else {
      return -2 - log(-$value)/log(2);
    }
  }
  return $value;
}

sub analyzeTPM {
  my ($f, $of) = @_;
  my $ofh = \*STDOUT;
  if (defined $of) {
    open($ofh, ">$of")
  }
  open(my $fh1, "<$f") || die "Can't open $f\n";
  my $h = <$fh1>; $h =~ s/[\r\n]//g;
  my @hdr = split("\t", $h);
  my $hh1 = {};
  my $expr = [];
  my $trpk = [ map { 0 } 2 .. $#hdr];
  my $total = [ map { 0 } 2 .. $#hdr];
  foreach my $i (0 .. $#hdr) { $hh1->{$hdr[$i]} = $i; }
  while (<$fh1>) {
    s/[\r\n]//g;
    my @ll = split("\t");
    next if ($ll[1] <= 0);
    my $e = [ map { $ll[$_]/$ll[1] } 2 .. $#ll];
    $trpk = [ map { $e->[$_] + $trpk->[$_] } 0 .. $#{$trpk}];
    $total = [ map { $total->[$_] + $ll[$_ + 2] } 0 .. $#{$total} ];
    push @$expr, [$ll[0], @$e];
  }
  close($fh1);
  print $ofh join("\t", $hdr[0], "Name", map { $hdr[$_] } 2 .. $#hdr), "\n";
  foreach my $e (@$expr) {
    print $ofh join("\t", $e->[0], $e->[0], map { $e->[$_ + 1]*1e6/($trpk->[$_]+1) } 0 .. $#{$trpk}), "\n";
    #print join("\t", $e->[0], $e->[0], map { $e->[$_ + 1] } 0 .. $#{$trpk}), "\n";
  }
  print $ofh join("\t", "Total", "Total", map { $total->[$_] / 1e6 } 0 .. $#{$total}), "\n";
  print STDERR join("\t", "Total", "Total", map { $total->[$_] / 1e6 } 0 .. $#{$total}), "\n";

  if (defined $of) {
    close($ofh);
  }
}

sub analyzeCPM {
  my ($f, $of) = @_;
  my $ofh = \*STDOUT;
  if (defined $of) {
    open($ofh, ">$of")
  }
  open(my $fh1, "<$f") || die "Can't open $f\n";
  my $h = <$fh1>; $h =~ s/[\r\n]//g;
  my @hdr = split("\t", $h);
  my $hh1 = {};
  my $total = [ map { 0 } 2 .. $#hdr];
  foreach my $i (0 .. $#hdr) { $hh1->{$hdr[$i]} = $i; }
  while (<$fh1>) {
    s/[\r\n]//g;
    my @ll = split("\t");
    $total = [ map { $total->[$_] + $ll[$_ + 2] } 0 .. $#{$total} ];
  }
  close($fh1);
  print $ofh join("\t", map { $hdr[$_] } 0 .. $#hdr), "\n";
  open(my $fh1, "<$f") || die "Can't open $f\n";
  my $h = <$fh1>; $h =~ s/[\r\n]//g;
  while (<$fh1>) {
    s/[\r\n]//g;
    my @ll = split("\t");
    print $ofh join("\t", $ll[0], $ll[1], map { $ll[$_ + 2]*1e6/($total->[$_]+1) } 0 .. $#{$total}), "\n";
  }
  close($fh1);
  print $ofh join("\t", "Total", "Total", map { $total->[$_] / 1e6 } 0 .. $#{$total}), "\n";
  print STDERR join("\t", "Total", "Total", map { $total->[$_] / 1e6 } 0 .. $#{$total}), "\n";

  if (defined $of) {
    close($ofh);
  }
}

sub countLen {
  my ($start, $end) = @_;
  my @l1 = grep { $_ ne "" } split(",", $start);
  my @l2 = grep { $_ ne "" } split(",", $end);
  my @order = sort { $l1[$a] <=> $l1[$b] } 0 .. $#l1;
  my @exon_start = map { $l1[$_] } @order;
  my @exon_end = map { $l2[$_] } @order;
  my $count = 0;
  my $max = $exon_start[0];
  my @list;
  for (my $i = 0; $i < scalar(@exon_start); $i++) {
    if ($max < $exon_start[$i]) {
      $max = $exon_start[$i];
    }
    if ($max < $exon_end[$i]) {
      $count += $exon_end[$i] - $max + 1;
      push @list, $exon_end[$i] - $max + 1;
      $max = $exon_end[$i];
    }
    else {
      push @list, 0;
    }
  }
  return $count;
}

1;
