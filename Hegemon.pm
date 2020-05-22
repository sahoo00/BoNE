# Author: Debashis Sahoo <dsahoo@ucsd.edu>
package Hegemon;

use IPC::Open2;
use U;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my %params = @_;
  my $self  = {};
  $self->{'expr'}  = undef;
  $self->{'idx'}  = undef;
  $self->{'thr'}  = undef;
  $self->{'survival'}  = undef;
  $self->{'fh'}  = undef;

  foreach (keys(%params)) {
    $self->{$_} = $params{$_};
  }

  &open($self, $self->{'expr'}, $self->{'idx'}, $self->{'thr'},
      $self->{'survival'});

  bless ($self, $class);
  return $self;
}

sub getIndex {
  my ($idxfile, $thrfile, $survivalfile) = @_;
  my $hash = {};
  if (defined $idxfile && $idxfile ne "") {
    open(my $fh, "<$idxfile") || die "Can't open $idxfile\n";
    my $h = <$fh>;
    while (<$fh>) {
      s/[\r\n]//g;
      my ($id, $ptr, $name, $desc) = split("\t");
      $hash->{"idx"}->{$id} = [$ptr, $name, $desc];
      my @names = split(" /// ", $name);
      foreach (@names) {
        if ($_ ne "") {
          push @{$hash->{"name"}->{$_}}, $id;
        }
      }
    }
    close($fh);
  }
  if (defined $thrfile && $thrfile ne "") {
    open(my $fh, "<$thrfile") || die "Can't open $thrfile\n";
    while (<$fh>) {
      s/[\r\n]//g;
      my ($id, $thr1, $stat, $thr0, $thr2) = split("\t");
      $hash->{"thr"}->{$id} = [$thr1, $stat, $thr0, $thr2];
    }
    close($fh);
  }
  if (defined $survivalfile && $survivalfile ne "") {
    open(my $fh, "<$survivalfile") || die "Can't open $survivalfile\n";
    while (<$fh>) {
      s/[\r\n]//g;
      my ($id, @list) = split("\t");
      $hash->{"survival"}->{$id} = [@list];
    }
    close($fh);
  }
  return $hash;
}

sub open {
  my ($self, $exprfile, $idxfile, $thrfile, $survivalfile) = @_;
  my $fh;
  my $hhash = {};
  my @headers;
  if (defined $exprfile) {
    open($fh, "<$exprfile") || die "Can't open $exprfile\n";
    my $head = <$fh>;
    $head =~ s/[\r\n]//g;
    @headers = split("\t", $head);
    for (my $i = 2; $i < scalar(@headers); $i++) {
      $hhash->{$headers[$i]} = $i;
    }
  }
  my $hash = &getIndex($idxfile, $thrfile, $survivalfile);
  $self->{'fh'} = $fh;
  $self->{'idHash'} = $hash->{'idx'};
  $self->{'nameHash'} = $hash->{'name'};
  $self->{'thrHash'} = $hash->{'thr'};
  $self->{'survivalHash'} = $hash->{'survival'};
  if (defined $survivalfile) {
    my $h = &U::getHeaders($survivalfile);
    $self->{'survivalHdrs'} = {};
    for (my $i = 0; $i <= scalar(@$h); $i++) {
      $self->{'survivalHdrs'}->{$h->[$i]} = $i;
    }
  }
  $self->{'headers'} = \@headers;
  $self->{'hhash'} = $hhash;

  my @ids = keys(%{$self->{"idHash"}});
  my @thr = map { $self->{'thrHash'}->{$_}->[1] } @ids;
  my @sthr = sort { $thr[$b] <=> $thr[$a] } 0 .. $#thr;
  my $pval = {};
  foreach my $i (0 .. $#sthr) {
    $pval->{$ids[$sthr[$i]]} = $i/scalar(@sthr);
  }
  $self->{'thrStatPval'} = $pval;
  $self->{'start'} = 2;
}

sub readPlatformFile {
  my ($self, $file) = @_;
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $h = <$fh>;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $id = $list[0];
    foreach my $i (1 .. $#list) {
      my $name = $list[$i];
      my @names = split(" /// ", $name);
      foreach my $n (@names) {
        $n =~ s/^\s*//g;
        $n =~ s/\s*$//g;
        $n = uc($n);
        next if ($n eq "" || $n eq "---");
        if (defined $self->{"nameHash"}->{$n}) {
          $self->{"nameHash"}->{$n} = &U::union([$id], $self->{"nameHash"}->{$n});
        }
        else {
          $self->{"nameHash"}->{$n} = [$id];
        }
      }
    }
  }
  close($fh);
}

sub close {
  my $self = shift;
  close($self->{'fh'});
}

sub getName {
  my ($self, $id) = @_;
  return undef if (!defined $self->{"idHash"});
  return undef if (!defined $self->{"idHash"}->{$id});
  return $self->{"idHash"}->{$id}->[1];
}

sub getDesc {
  my ($self, $id) = @_;
  return undef if (!defined $self->{"idHash"});
  return undef if (!defined $self->{"idHash"}->{$id});
  return $self->{"idHash"}->{$id}->[2];
}

sub getSimpleName {
  my ($self, $id) = @_;
  my $name = $self->getName($id);
  return undef if (!defined $name);
  $name = (split ":", $name)[0];
  $name = (split " /// ", $name)[0];
  return $name;
}

sub getIDs {
  my ($self, $name) = @_;
  if (defined $self->{"idHash"} && defined $self->{"idHash"}->{$name}) {
    return [$name];
  }
  if (defined $self->{"nameHash"} && defined $self->{"nameHash"}->{$name}) {
    return $self->{"nameHash"}->{$name};
  }
  return undef;
}

sub getStart {
  my $self = shift;
  return $self->{'start'};
}

sub getEnd {
  my $self = shift;
  my $end = scalar(@{$self->{'headers'}}) - 1;
  return $end;
}

sub getNum {
  my $self = shift;
  my $num = scalar(@{$self->{'headers'}});
  return $num;
}

sub getExprData {
  my ($self, $id) = @_;
  return undef if (!defined $self->{"idHash"});
  return undef if (!defined $self->{"idHash"}->{$id});
  my $ptr = $self->{"idHash"}->{$id}->[0];
  return &getData($self->{'fh'}, $ptr);  
}

sub getData {
   my ($fh, $ptr) = @_;
   &my_seek($fh, $ptr);
   my $in = <$fh>;
   $in =~ s/[\r\n]//g;
   my @list = split("\t", $in, -1);
   return [@list];
}

sub my_seek {
  my ($fh, $x) = @_;
  seek($fh, $x, 0);
  $pos = tell($fh);
  if ($pos != $x) {
    $off = $x - $pos;
    if ($off > 0) {
      $chunk = 1000000000;
      $num = int($off / $chunk);
      for ($i = 0; $i < $num; $i++) {
        seek($fh, $chunk, 0);
      }
      $off = $off - $num * $chunk;
      seek($fp, $off, 0);
    }
    else {
      print STDERR "Error in seek $x < $pos\n";
      exit(1);
    }
  }
}

# Get Patient IDs that have value $str in the $index column
sub getSurvivalArray {
  my ($self, $index, $str) = @_;
  return undef if (!defined $self->{"survivalHash"});
  $res = [];
  for (my $i = 0; $i < scalar(@{$self->{'headers'}}); $i++) {
    my $arr = $self->{'headers'}->[$i];
    if ($self->{'survivalHash'}->{$arr}->[$index] eq $str) {
      push @$res, $i;
    }
  }
  return $res;
}

# Get values of the $index column
sub getSurvivalData {
  my ($self, $index) = @_;
  return undef if (!defined $self->{'survivalHash'});
  $res = [];
  for (my $i = 0; $i < scalar(@{$self->{'headers'}}); $i++) {
    my $arr = $self->{'headers'}->[$i];
    $res->[$i] = $self->{'survivalHash'}->{$arr}->[$index];
  }
  return $res;
}

sub getSurvivalName {
  my ($self, $name) = @_;
  if (!defined $self->{'survivalHdrs'}) {
    die "Can't find survival headers\n";
  }
  if (!defined $self->{'survivalHdrs'}->{$name}) {
    die "Can't find survival headers $name\n";
  }
  return $self->getSurvivalData($self->{'survivalHdrs'}->{$name} - 1);
}

sub getTimeData {
  my $self = shift;
  return $self->getSurvivalData(0);
}

sub getStatusData {
  my $self = shift;
  return $self->getSurvivalData(1);
}

sub getThrData {
  my ($self, $id) = @_;
  return undef if (!defined $self->{"thrHash"});
  return undef if (!defined $self->{"thrHash"}->{$id});
  return $self->{"thrHash"}->{$id};
}

sub getBinaryData {
  my ($self, $id, $thr) = @_;
  my $expr = $self->getExprData($id);
  my $thr_step = $self->getThrData($id);
  my $t = &getThrCode($thr_step, $thr_step->[0], $thr);
  my $start = $self->getStart();
  my $end = $self->getEnd();
  my $res = [$expr->[0], $expr->[1]];
  for(my $i = $start; $i <= $end; $i++) {
    next if (!defined $expr->[$i] || $expr->[$i] eq "");
    if ($expr->[$i] >= $t) { $res->[$i] = 1; }
    if ($expr->[$i] <  $t) { $res->[$i] = 0; }
  }
  return $res;
}

sub getTertiaryData {
  my ($self, $id, $thr, $gap) = @_;
  $gap = 0.5 if (!defined $gap);
  my $expr = $self->getExprData($id);
  my $thr_step = $self->getThrData($id);
  my $t = &getThrCode($thr_step, $thr_step->[0], $thr);
  my $start = $self->getStart();
  my $end = $self->getEnd();
  my $res = [$expr->[0], $expr->[1]];
  for(my $i = $start; $i <= $end; $i++) {
    next if (!defined $expr->[$i] || $expr->[$i] eq "");
    if ($expr->[$i] >= ($t+$gap)) { $res->[$i] = 2; }
    if ($expr->[$i] <  ($t+$gap) && $expr->[$i] >= ($t-$gap)) { $res->[$i] = 1; }
    if ($expr->[$i] <  ($t-$gap)) { $res->[$i] = 0; }
  }
  return $res;
}

sub getGroupsData {
  my ($self, $id, $groups) = @_;
  my $expr = $self->getExprData($id);
  my $res = [];
  foreach my $g (@$groups) {
    my $datax = [];
    foreach my $i (@{$g->[2]}) {
      next if (!defined $expr->[$i] || $expr->[$i] eq "");
      push @$datax, $expr->[$i];
    }
    push @$res, $datax;
  }
  return $res;
}

sub getCox {
  my ($self, $id) = @_;
  my $expr = $self->getExprData($id);
  my $time = $self->getTimeData();
  my $status = $self->getStatusData();
  my $thr = $self->getThrData($id);
  my $res = [];
  for(my $i = 2; $i < scalar(@{$self->{'headers'}}); $i++) {
    next if (!defined $expr->[$i] || $expr->[$i] eq "");
    next if (!defined $time->[$i] || $time->[$i] eq "");
    next if (!defined $status->[$i] || $status->[$i] eq "");
    push @$res, [$expr->[$i], $time->[$i], $status->[$i]];
  }

  my ($phr, $phw);
  my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  print $phw "library(survival)\n";
  print $phw "rm(list=ls())\n";
  print $phw "expr <- c(", join(",", map {$_->[0]} @$res), ")\n";
  print $phw "time <- c(", join(",", map {$_->[1]} @$res), ")\n";
  print $phw "status <- c(", join(",", map {$_->[2]} @$res), ")\n";
  print $phw "groups <- as.numeric(expr < ", $thr->[0], ") + 2 * as.numeric(expr >= ", $thr->[0], ")\n";
print $phw <<END
idx <- (!is.na(groups) & groups > 0)
time <- time[idx]
status <- status[idx]
groups <- groups[idx]
x <- coxph(Surv(time, status) ~ groups)
s <- summary(x)
codes <- c("0", "***", "**", "*", ".", " ")
p <- s\$coefficients[,5]
pg <- as.numeric(p <= 0) + 2 * as.numeric(p > 0 & p <= 0.001) +
    3 * as.numeric(p > 0.001 & p <= 0.01) +
    4 * as.numeric(p > 0.01 & p <= 0.05) +
    5 * as.numeric(p > 0.05 & p <= 0.1) +
    6 * as.numeric(p > 0.1)
cat("exp(coef)", "lower .95", "upper .95", "p", "c", "\\n", sep="\\t")
cat(s\$conf.int[1], s\$conf.int[3], s\$conf.int[4], p, codes[pg], "\\n", sep="\\t")
END
;
  my $h = <$phr>;
  my $h = <$phr>;
  $h =~ s/[\r\n]//g;
  my @list = split("\t", $h, -1);
  close($phr);
  close($phw);
  return [@list[0 .. 4]];
}

sub getHR {
  my ($self, $id) = @_;
  my $res = $self->getCox($id);
  return $res->[0];
}

sub getHRList {
  my ($self, $file, $index) = @_;
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $hash = {};
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $id = $list[$index];
    my $hr = $self->getHR($id);
    my $n = $self->getName($id);
    #print "$id\t$hr\t$n\n";
    $hash->{$id} = [$hr, @list];
  }
  close($fh);
  return $hash;
}

sub plotPair {
  my ($self, $outfile, $id1, $id2, $groups, %params) = @_;
  my $marker = "+";
  if (defined $params{"marker"}) {
    $marker = $params{"marker"};
  }
  my $dpi = 100;
  if (defined $params{"dpi"}) {
    $dpi = $params{"dpi"};
  }
  my $mew = 1.1;
  my $ms = 4;
  $mew = $params{"mew"} if (defined $params{"mew"});
  $ms = $params{"ms"} if (defined $params{"ms"});
  my $expr1 = $self->getExprData($id1);
  my $expr2 = $self->getExprData($id2);
  my $xn = $self->getName($id1);
  my $yn = $self->getName($id2);
  my ($x_id, $y_id) = ($id1, $id2);
  my ($phr, $phw);
  my $pid = open2($phr, $phw, "python");
  print $phw <<END
import matplotlib
matplotlib.use('agg')
import re
from pylab import *
from numpy import *
END
;
  if ($outfile =~ /.pdf$/) {
  $x_id =~ s/_/\\_/g;
  $y_id =~ s/_/\\_/g;
  print $phw <<END
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
END
;
  }
  if (defined $params{"square"}) {
  print $phw <<END
fig = figure(figsize=(4.8,4.8))
ax = fig.add_axes([54.0/480, 54.0/480, 1-2*54.0/480, 1-2*54.0/480])
END
;
  }
  else {
  print $phw <<END
fig = figure(figsize=(6.4,4.8))
ax = fig.add_axes([70.0/640, 54.0/480, 1-2*70.0/640, 1-2*54.0/480])
END
;
  }
  if (!defined $groups) {
    my ($datax, $datay) = ([], []);
    for(my $i = 2; $i < scalar(@{$self->{'headers'}}); $i++) {
      next if (!defined $expr1->[$i] || $expr1->[$i] eq "");
      next if (!defined $expr2->[$i] || $expr2->[$i] eq "");
      push @$datax, $expr1->[$i];
      push @$datay, $expr2->[$i];
    }
    print $phw "datax = [", join(",", @$datax), "]\n";
    print $phw "datay = [", join(",", @$datay), "]\n";
    print $phw "c = 'blue'\n";
    print $phw "ax.plot(datax,datay, color=c, ls='None', marker='$marker', mew=$mew, ms=$ms, mec=c)\n";
    print $phw "ax.axis([min(datax)-0.5, max(datax)+0.5, min(datay)-0.5, max(datay)+0.5])\n";
  }
  else {
    print $phw "lim = [100, -100, 100, -100]\n";
    foreach my $g (@$groups) {
      my ($datax, $datay) = ([], []);
      foreach my $i (@{$g->[2]}) {
        next if (!defined $expr1->[$i] || $expr1->[$i] eq "");
        next if (!defined $expr2->[$i] || $expr2->[$i] eq "");
        push @$datax, $expr1->[$i];
        push @$datay, $expr2->[$i];
      }
      next if (scalar(@$datax) == 0 || scalar(@$datay) == 0);
      print $phw "datax = [", join(",", @$datax), "]\n";
      print $phw "datay = [", join(",", @$datay), "]\n";
      print $phw "c = '", $g->[1], "'\n";
      print $phw "lim = [min(lim[0], min(datax)), max(lim[1], max(datax)), min(lim[2], min(datay)), max(lim[3], max(datay))]\n";
      print $phw "ax.plot(datax,datay, color=c, ls='None', marker='$marker', mew=$mew, ms=$ms, mec=c)\n";
    }
    print $phw "ax.axis([lim[0]-0.5, lim[1]+0.5, lim[2]-0.5, lim[3]+0.5])\n";
  }
  print $phw <<END
ax.set_xlabel('$x_id:  $xn', fontsize=10)
ax.set_ylabel('$y_id:  $yn', fontsize=10)
fig.savefig('$outfile', dpi=$dpi)
END
;
  close($phr);
  close($phw);
  sleep(1);
  chmod 0666, $outfile;
}

sub plotTikz {
  my ($self, $ofile, $id1, $id2, $groups, %params) = @_;
  open(my $ofh, ">$ofile") || die "Can't write $ofile\n";
  my $marker = "+";
  if (defined $params{"marker"}) {
    $marker = $params{"marker"};
  }
  my $mew = 1.1;
  my $ms = 4;
  $mew = $params{"mew"} if (defined $params{"mew"});
  $ms = $params{"ms"} if (defined $params{"ms"});
  my $expr1 = $self->getExprData($id1);
  my $expr2 = $self->getExprData($id2);
  my $thr1 = $self->getThrData($id1);
  my $thr2 = $self->getThrData($id2);
  $thr1 = $params{"thr1"} if (defined $params{"thr1"});
  $thr2 = $params{"thr2"} if (defined $params{"thr2"});
  my $xn = $self->getName($id1);
  my $yn = $self->getName($id2);
  my ($x_id, $y_id) = ($id1, $id2);
  my ($datax, $datay) = ([], []);
  if (!defined $groups) {
    my ($g, $c) = (0, 'a');
    my $clr = &U::getPScolor("clr$g", 'blue');
    print $ofh "$clr\n";
    print $ofh "\\pgfplotstableread{%\n";
    for(my $i = 2; $i < scalar(@{$self->{'headers'}}); $i++) {
      next if (!defined $expr1->[$i] || $expr1->[$i] eq "");
      next if (!defined $expr2->[$i] || $expr2->[$i] eq "");
      push @$datax, $expr1->[$i];
      push @$datay, $expr2->[$i];
      print $ofh $expr1->[$i]." ".$expr2->[$i]."\n";
    }
    print $ofh "}\\g$c"."data%\n";
  }
  else {
    for (my $g = 0; $g < scalar(@$groups); $g++) {
      my $c = chr(ord('a')+$g);
      my $clr = &U::getPScolor("clr$g", $groups->[$g]->[1]);
      print $ofh "$clr\n";
      print $ofh "\\pgfplotstableread{%\n";
      foreach my $i (@{$groups->[$g]->[2]}) {
	next if (!defined $expr1->[$i] || $expr1->[$i] eq "");
	next if (!defined $expr2->[$i] || $expr2->[$i] eq "");
	push @$datax, $expr1->[$i];
	push @$datay, $expr2->[$i];
	print $ofh $expr1->[$i]." ".$expr2->[$i]."\n";
      }
      print $ofh "}\\g$c"."data%\n";
    }
  }
  my ($thrx0, $thrx1, $thrx2) = ($thr1->[2], $thr1->[0], $thr1->[3]);
  my ($thry0, $thry1, $thry2) = ($thr2->[2], $thr2->[0], $thr2->[3]);
  my $clr = &U::getPScolor("thr0", 'red');
  print $ofh "$clr\n";
  my $clr = &U::getPScolor("thr1", 'cyan');
  print $ofh "$clr\n";
  my $maxx = &U::max($datax) + 0.5;
  my $maxy = &U::max($datay) + 0.5;
  my $minx = &U::min($datax) - 0.5;
  my $miny = &U::min($datay) - 0.5;
  my $lmaxx = ($maxx - $minx);
  my $lmaxy = ($maxy - $miny);
  if ($lmaxx <= 0) { $lmaxx = 1; }
  if ($lmaxy <= 0) { $lmaxy = 1; }
  my $xunit = (640 - 2 * 70)/$lmaxx/100;
  my $yunit = (480 - 2 * 54)/$lmaxy/100;
  my $ox = 70/100/$xunit;
  my $oy = 54/100/$yunit;
  my $tx = (640 - 70)/100/$xunit;
  my $ty = (480 - 54)/100/$yunit;
  my $xl = &U::myescape("$x_id: $xn");
  my $yl = &U::myescape("$y_id: $yn");
  print $ofh <<END
\\def\\ttsize{\\tiny}%
\\def\\ttsizea{\\fontsize{9pt}{9pt}\\selectfont}%
\\def\\ttsizeb{\\fontsize{8pt}{8pt}\\selectfont}%
\\def\\ttsizec{\\fontsize{7pt}{7pt}\\selectfont}%
\\def\\pshlabel#1{\\ttsize #1}%
\\def\\psvlabel#1{\\ttsize #1}%
\\begin{tikzpicture}[x=$xunit in, y=$yunit in]
\\begin{axis}[x=$xunit in, y=$yunit in,
  axis x line = bottom,axis y line = left,
  ticklabel style={font=\\tiny},
  ymin=$miny, ymax=$maxy, xmin=$minx, xmax=$maxx,
  xlabel={\\bf \\ttsizeb $xl},
  ylabel={\\bf \\ttsizea $yl},
  xlabel style={below,anchor=north},
  ylabel style={left,anchor=south}]
\\draw[black,line width=0.75pt] ($minx, $miny) rectangle ($maxx, $maxy);
END
;
  if (!defined $groups) {
    my $dname = "\\gadata";
    my $g = 0;
    print $ofh "\\addplot+[line width=1.3pt, color=clr$g, mark=*, only marks,mark options={color=clr$g}]  table [x index=0, y index=1] {$dname};\n";
  }
  else {
    for (my $g = 0; $g < scalar(@$groups); $g++) {
      my $c = chr(ord('a')+$g);
      my $dname = "\\g".$c."data";
      print $ofh "\\addplot+[line width=1.3pt, color=clr$g, mark=*, only marks,mark options={color=clr$g}]  table [x index=0, y index=1] {$dname};\n";
    }
  }
  print $ofh "\\draw[thr1,line width=0.5pt] ($thrx0, $miny) -- ($thrx0, $maxy);\n";
  print $ofh "\\draw[thr0,line width=0.5pt] ($thrx1, $miny) -- ($thrx1, $maxy);\n";
  print $ofh "\\draw[thr1,line width=0.5pt] ($thrx2, $miny) -- ($thrx2, $maxy);\n";
  print $ofh "\\draw[thr1,line width=0.5pt] ($minx, $thry0) -- ($maxx, $thry0);\n";
  print $ofh "\\draw[thr0,line width=0.5pt] ($minx, $thry1) -- ($maxx, $thry1);\n";
  print $ofh "\\draw[thr1,line width=0.5pt] ($minx, $thry2) -- ($maxx, $thry2);\n";
  print $ofh "\\end{axis}\n";
  print $ofh "\\end{tikzpicture}\n";
  close($ofh);
  chmod 0666, $ofile;
}

sub countQs {
  my ($datax, $datay, $thrx0, $thrx2, $thry0, $thry2) = @_;
  my $num1 = scalar(@$datax);
  my $num2 = scalar(@$datay);
  if ($num1 != $num2) {
    die "Error datax($num1) and datay($num2) sizes don't match\n";
  }
  my $res = [0, 0, 0, 0];
  for (my $i = 0; $i < $num1; $i++) {
    if ($datax->[$i] < $thrx0 && $datay->[$i] < $thry0) { $res->[0] ++; }
    if ($datax->[$i] < $thrx0 && $datay->[$i] >= $thry2) { $res->[1] ++; }
    if ($datax->[$i] >= $thrx2 && $datay->[$i] < $thry0) { $res->[2] ++; }
    if ($datax->[$i] >= $thrx2 && $datay->[$i] >= $thry2) { $res->[3] ++; }
  }
  return $res;
}

sub convertBooleanStats {
   my ($c0, $c1, $c2, $c3)= @_;
   $total = $c0 + $c1 + $c2 + $c3;
   if ($total <= 0) {
     return [[0,0,0,0], [0,0,0,0], [0,0,0,0]];
   }
   $e0 = ($c0 + $c1) * ($c0 + $c2) /$total;
   $e1 = ($c1 + $c0) * ($c1 + $c3) /$total;
   $e2 = ($c2 + $c0) * ($c2 + $c3) /$total;
   $e3 = ($c3 + $c1) * ($c3 + $c2) /$total;
   $s0 = (int($e0) + 1 - $c0)/sqrt(int($e0) + 1);
   $s1 = (int($e1) + 1 - $c1)/sqrt(int($e1) + 1);
   $s2 = (int($e2) + 1 - $c2)/sqrt(int($e2) + 1);
   $s3 = (int($e3) + 1 - $c3)/sqrt(int($e3) + 1);
   $p0 = ($c0/($c0 + $c1 + 1) + $c0/($c0 + $c2 + 1))/2;
   $p1 = ($c1/($c1 + $c0 + 1) + $c1/($c1 + $c3 + 1))/2;
   $p2 = ($c2/($c2 + $c0 + 1) + $c2/($c2 + $c3 + 1))/2;
   $p3 = ($c3/($c3 + $c1 + 1) + $c3/($c3 + $c2 + 1))/2;
   return [[$e0, $e1, $e2, $e3], [$s0, $s1, $s2, $s3], [$p0, $p1, $p2, $p3]];
}

sub getThrCode {
  my ($thr_step, $value, $code) = @_;
  return $value if (!defined $code);
  my $thr = $value;
  if ($code eq "thr1") {
    $thr = $thr_step->[0];
  }
  elsif ($code eq "thr0") {
    $thr = $thr_step->[2];
  }
  elsif ($code eq "thr2") {
    $thr = $thr_step->[3];
  }
  else {
    $thr = $code;
  }
  return $thr;
}

sub getBooleanStats {
  my ($self, $id1, $id2, $groups, $thrx0, $thrx2, $thry0, $thry2) = @_;
  my $expr1 = $self->getExprData($id1);
  my $expr2 = $self->getExprData($id2);
  my $thr1 = $self->getThrData($id1);
  my $thr2 = $self->getThrData($id2);
  $thrx0 = &getThrCode($thr1, $thr1->[2], $thrx0);
  $thrx2 = &getThrCode($thr1, $thr1->[3], $thrx2);
  $thry0 = &getThrCode($thr2, $thr2->[2], $thry0);
  $thry2 = &getThrCode($thr2, $thr2->[3], $thry2);
  my $qs = [0, 0, 0, 0];
  if (!defined $groups) {
    my ($datax, $datay) = ([], []);
    for(my $i = 2; $i < scalar(@{$self->{'headers'}}); $i++) {
      next if (!defined $expr1->[$i] || $expr1->[$i] eq "");
      next if (!defined $expr2->[$i] || $expr2->[$i] eq "");
      push @$datax, $expr1->[$i];
      push @$datay, $expr2->[$i];
    }
    $qs = &countQs($datax, $datay, $thrx0, $thrx2, $thry0, $thry2);
  }
  else {
    my ($dx, $dy) = ([], []);
    foreach my $g (@$groups) {
      foreach my $i (@{$g->[2]}) {
        next if (!defined $expr1->[$i] || $expr1->[$i] eq "");
        next if (!defined $expr2->[$i] || $expr2->[$i] eq "");
        push @$dx, $expr1->[$i];
        push @$dy, $expr2->[$i];
      }
    }
    $qs = &countQs($dx, $dy, $thrx0, $thrx2, $thry0, $thry2);
  }
  my ($c0, $c1, $c2, $c3) = @$qs;
  my $bs = &convertBooleanStats($c0, $c1, $c2, $c3);
  return [$qs, @$bs];
}

sub plotBooleanPair {
  my ($self, $outfile, $id1, $id2, $groups, $thr1, $thr2, $mew, $ms, %params) = @_;
  my $marker = "+";
  if (defined $params{"marker"}) {
    $marker = $params{"marker"};
  }
  my $dpi = 100;
  if (defined $params{"dpi"}) {
    $dpi = $params{"dpi"};
  }
  $mew = 1 if (!defined $mew);
  $ms = 4 if (!defined $ms);
  $thr1 = $self->getThrData($id1) if (!defined $thr1);
  $thr2 = $self->getThrData($id2) if (!defined $thr2);
  my $expr1 = $self->getExprData($id1);
  my $expr2 = $self->getExprData($id2);
  my ($thrx1, $sx, $thrx0, $thrx2) = @$thr1; 
  my ($thry1, $sy, $thry0, $thry2) = @$thr2; 
  my $qs = [0, 0, 0, 0];
  my $xn = $self->getName($id1);
  my $yn = $self->getName($id2);
  my ($x_id, $y_id) = ($id1, $id2);
  my ($phr, $phw);
  if (defined $params{"ph"}) {
    ($phr, $phw) = @{$params{"ph"}};
  }
  else {
    my $pid = open2($phr, $phw, "python");
  }
  my $blw = 1;
  if (defined $params{"blw"}) {
    $blw = $params{"blw"};
  }
  #$phw = \*STDOUT;
  print $phw <<END
import matplotlib
matplotlib.use('agg')
import re
from pylab import *
from numpy import *
rcParams['axes.linewidth'] = $blw
END
;
  if ($outfile =~ /.pdf$/) {
  $x_id =~ s/_/\\_/g;
  $y_id =~ s/_/\\_/g;
  print $phw <<END
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
END
;
  }
  if (defined $params{"square"}) {
  print $phw <<END
fig = figure(figsize=(4.8,4.8))
ax = fig.add_axes([54.0/480, 54.0/480, 1-2*54.0/480, 1-2*54.0/480])
END
;
  }
  else {
  print $phw <<END
fig = figure(figsize=(6.4,4.8))
ax = fig.add_axes([70.0/640, 54.0/480, 1-2*70.0/640, 1-2*54.0/480])
END
;
  }
  if (!defined $groups) {
    my ($datax, $datay) = ([], []);
    for(my $i = 2; $i < scalar(@{$self->{'headers'}}); $i++) {
      next if (!defined $expr1->[$i] || $expr1->[$i] eq "");
      next if (!defined $expr2->[$i] || $expr2->[$i] eq "");
      push @$datax, $expr1->[$i];
      push @$datay, $expr2->[$i];
    }
    $qs = &countQs($datax, $datay, $thrx0, $thrx2, $thry0, $thry2);
    print $phw "datax = [", join(",", @$datax), "]\n";
    print $phw "datay = [", join(",", @$datay), "]\n";
    print $phw "c = 'blue'\n";
    print $phw "ax.plot(datax,datay, color=c, ls='None', marker='$marker', mew=$mew, ms=$ms, mec=c)\n";
    print $phw "lim = [min(datax), max(datax), min(datay), max(datay)]\n";
  }
  else {
    print $phw "lim = [100, -100, 100, -100]\n";
    my ($dx, $dy) = ([], []);
    foreach my $g (@$groups) {
      my ($datax, $datay) = ([], []);
      foreach my $i (@{$g->[2]}) {
        next if (!defined $expr1->[$i] || $expr1->[$i] eq "");
        next if (!defined $expr2->[$i] || $expr2->[$i] eq "");
        push @$datax, $expr1->[$i];
        push @$datay, $expr2->[$i];
        push @$dx, $expr1->[$i];
        push @$dy, $expr2->[$i];
      }
      next if (scalar(@$datax) == 0 || scalar(@$datay) == 0);
      print $phw "datax = [", join(",", @$datax), "]\n";
      print $phw "datay = [", join(",", @$datay), "]\n";
      print $phw "c = '", $g->[1], "'\n";
      print $phw "lim = [min(lim[0], min(datax)), max(lim[1], max(datax)), min(lim[2], min(datay)), max(lim[3], max(datay))]\n";
      print $phw "ax.plot(datax,datay, color=c, ls='None', marker='$marker', mew=$mew, ms=$ms, mec=c)\n";
    }
    $qs = &countQs($dx, $dy, $thrx0, $thrx2, $thry0, $thry2);
  }
  my ($c0, $c1, $c2, $c3) = @$qs;
  my $bs = &convertBooleanStats($c0, $c1, $c2, $c3);
  my ($e0, $e1, $e2, $e3) = map { sprintf("%.3f", $_) } @{$bs->[0]};
  my ($s0, $s1, $s2, $s3) = map { sprintf("%.3f", $_) } @{$bs->[1]};
  my ($p0, $p1, $p2, $p3) = map { sprintf("%.3f", $_) } @{$bs->[2]};
  my @cls = ('#cccccc', '#cccccc', '#cccccc', '#cccccc');
  for (my $i = 0; $i < 4; $i++) {
    if ($bs->[1]->[$i] > 2 && $bs->[2]->[$i] < 0.1) { $cls[$i] = 'r'; }
  }
  my ($cl0, $cl1, $cl2, $cl3) = @cls;
  #print "$xn,$yn,$x_id,$y_id,$s0,$p0,$s1,$p1,$s2,$p2,$s3,$p3\n";
  my $tcl1 = '#dddddd';
  my $tcl2 = '#aaaaaa';
  my $tcl1 = '#777777';
  my $tcl2 = '#777777';
  my $tcl1 = 'c';
  my $tcl2 = 'r';

  if (defined $params{'limit'}) {
    my ($minx, $maxx, $miny, $maxy) = @{$params{'limit'}};
  print $phw <<END
minx,maxx,miny,maxy=$minx, $maxx, $miny, $maxy
END
;
  }
  else {
  print $phw <<END
minx,maxx,miny,maxy=0, 16, 0, 16
minx,maxx,miny,maxy=lim[0]-0.5, lim[1]+0.5, lim[2]-0.5, lim[3]+0.5
END
;
  }
  if (defined $params{"insert"}) {
    print $phw $params{"insert"}, "\n";
  }
  my $lw = 1;
  my $lfs = 10;
  if (defined $params{"lw"}) {
    $lw = $params{"lw"};
  }
  if (defined $params{"lfs"}) {
    $lfs = $params{"lfs"};
  }
  print $phw <<END
ax.plot([minx,maxx], [$thry0, $thry0],  '$tcl1', linewidth=$lw)
ax.plot([$thrx0, $thrx0], [miny, maxy], '$tcl1', linewidth=$lw)
ax.plot([minx,maxx], [$thry1, $thry1],  '$tcl2', linewidth=$lw)
ax.plot([$thrx1, $thrx1], [miny, maxy], '$tcl2', linewidth=$lw)
ax.plot([minx,maxx], [$thry2, $thry2],  '$tcl1', linewidth=$lw)
ax.plot([$thrx2, $thrx2], [miny, maxy], '$tcl1', linewidth=$lw)
#text(0.02, 0.02,'$c0,$e0,$s0,$p0', color='$cl0', ha='left', va='bottom', transform = ax.transAxes)
#text(0.02, 0.98,'$c1,$e1,$s1,$p1', color='$cl1', ha='left', va='top', transform = ax.transAxes)
#text(0.98, 0.02,'$c2,$e2,$s2,$p2', color='$cl2', ha='right', va='bottom', transform = ax.transAxes)
#text(0.98, 0.98,'$c3,$e3,$s3,$p3', color='$cl3', ha='right', va='top', transform = ax.transAxes)
ax.axis([minx,maxx,miny,maxy])
ax.set_xlabel('$x_id:  $xn', fontsize=$lfs)
ax.set_ylabel('$y_id:  $yn', fontsize=$lfs)
fig.savefig('$outfile', dpi=$dpi)
END
;
  close($phr);
  close($phw);
  sleep(1);
  chmod 0666, $outfile;
}

sub writeSelectedExpr {
  my ($self, $outexprfile, $outidxfile, $order) = @_;
  my $fh = $self->{'fh'};
  seek($fh, 0, 0);
  open(my $oexprfh, ">$outexprfile") || die "Can't open $outexprfile\n";
  open(my $oidxfh, ">$outidxfile") || die "Can't open $outidxfile\n";
  print $oidxfh "ProbeID\tPtr\tName\tDescription\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @header = split("\t", $head);
  $head = join("\t", @header[0,1], map {$header[$_]} @$order)."\n";
  print $oexprfh $head;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $ptr = tell($oexprfh);
    my @info = @{$self->{"idHash"}->{$list[0]}};
    print $oidxfh join("\t", $list[0], $ptr, @info[1,2]), "\n";
    print $oexprfh join("\t", @list[0,1], map {$list[$_]} @$order), "\n";
  }
  close($oexprfh);
  close($oidxfh);
}

sub writeSelectedSurvival {
  my ($self, $outsurvivalfile, $survivalfile, $order) = @_;
  my $hash = {};
  foreach my $i (@$order) {
    my $arr = $self->{'headers'}->[$i];
    $hash->{$arr} = 1;
  }
  open(my $fh, "<$survivalfile") || die "Can't open $survivalfile\n";
  open(my $ofh, ">$outsurvivalfile") || die "Can't open $outsurvivalfile\n";
  my $head = <$fh>;
  print $ofh $head;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    if (defined $hash->{$list[0]}) {
      print $ofh join("\t", @list), "\n";
    }
  }
  close($ofh);
  close($fh);
}

sub writeSelectedIh {
  my ($self, $outihfile, $order) = @_;
  my $if = $outihfile;
  my $headers = $self->{'headers'};
  my $ofh;
  open($ofh, ">$if") || die "Can't open $if\n";
  print $ofh "ArrayID\tArrayHeader\tClinicalhHeader\n";
  foreach my $i (@$order) {
    my $arr = $headers->[$i];
    print $ofh "$arr\t$arr\t$arr\n";
  }
  close($ofh);
}

sub writeSelectedPCL {
  my ($ofile, $file, $order) = @_;
  open(my $ofh, ">$ofile") || die "Can't write $ofile\n";
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  if (!defined $order) {
    print $ofh join("\t", @headers[0, 1]), "\tGWEIGHT\t", join("\t", @headers[2 .. $#headers]), "\n";
    print $ofh join("\t", "EWEIGHT", ""), "\t1\t", join("\t", map { 1 } @headers[2 .. $#headers]), "\n";
  }
  else {
    print $ofh join("\t", @headers[0, 1]), "\tGWEIGHT\t", 
          join("\t", map {$headers[$_]} @$order), "\n";
    print $ofh join("\t", "EWEIGHT", ""), "\t1\t", 
          join("\t", map { 1 } @$order), "\n";
  }
 
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    if (!defined $order) {
      print $ofh join("\t", @list[0, 1]), "\t1\t", join("\t", @list[2 .. $#headers]), "\n";
    }
    else {
      print $ofh join("\t", @list[0, 1]), "\t1\t", join("\t", map {$list[$_]} @$order), "\n";
    }
  }
  close($fh);
  close($ofh);
}

sub writeSelectedThr {
  my ($ofile, $file, $order, $gap) = @_;
  open(my $ofh, ">$ofile") || die "Can't write $ofile\n";
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my @data = sort { $a <=> $b } map {$list[$_]} @$order;
    my $res = &U::fitStep(\@data, 0, $#data);
    print $ofh join("\t", $list[0], $res->[6], $res->[3], $res->[6]-$gap,
        $res->[6]+$gap), "\n";
  }
  close($fh);
  close($ofh);
}

sub getArraysRank {
  my ($self, $id, $rank1, $rank2, $type) = @_;
  $type = "abs" if (!defined $type);
  my $res = [];
  my $start = $self->getStart();
  my $end = $self->getEnd();
  my $len = $end - $start + 1;
  my $expr = $self->getExprData($id);
  my @sorted = sort { $expr->[$a] <=> $expr->[$b] } $start .. $end;
  if ($type eq "abs") {
    $rank1 = ($len + $rank1) % $len;
    $rank2 = ($len + $rank2) % $len;
    if ($rank1 > $rank2) {
      my $tmp = $rank1; $rank1 = $rank2; $rank2 = $tmp;
    }
    $res = [ map { $sorted[$_] } $rank1 .. $rank2 ];
  }
  if ($type eq "perc") {
    $rank1 = int( ($len - 1) * $rank1) % $len;
    $rank2 = int( ($len - 1) * $rank2) % $len;
    if ($rank1 > $rank2) {
      my $tmp = $rank1; $rank1 = $rank2; $rank2 = $tmp;
    }
    $res = [ map { $sorted[$_] } $rank1 .. $rank2 ];
  }
  return $res;
}

sub getArraysThr {
  my ($self, $id, $thr, $type) = @_;
  my $res = [];
  my $start = 2;
  my $end = scalar(@{$self->{'headers'}}) - 1;
  my $expr = $self->getExprData($id);
  my $thr_step = $self->getThrData($id);
  $thr = &getThrCode($thr_step, $thr_step->[0], $thr);
  for (my $i = $start; $i <= $end ; $i++) {
    if (!defined $thr) { push @$res, $i; }
    elsif (!defined $expr->[$i] || $expr->[$i] eq "") { next; }
    elsif ($type eq "hi" && $expr->[$i] >= $thr) { push @$res, $i; }
    elsif ($type eq "lo" && $expr->[$i] <  $thr) { push @$res, $i; }
    elsif (defined $type && $type ne "lo" && $type ne "hi" && 
        $expr->[$i] >= $thr && $expr->[$i] <= $type) { 
      push @$res, $i;
    }
  }
  return $res;
}

sub getArraysAll {
  my ($self, @data) = @_;
  my $res = $self->getArraysThr();
  for (my $i = 0; $i < scalar(@data); $i+=3) {
    my $r = $self->getArraysThr($data[$i], $data[$i+1], $data[$i+2]);
    $res = &U::intersection($res, $r);
  }
  return $res;
}

sub filterSurvivalData {
  my ($self, $groups) = @_;
  my $time =   $self->getTimeData();
  my $status = $self->getStatusData();
  my $res = [];
  for (my $i = 0; $i < scalar(@$groups); $i++) {
    my $g = $groups->[$i];
    my $count = 0;
    my $r = [];
    foreach my $j (@{$g->[2]}) {
      next if (!defined $time->[$j] || $time->[$j] eq "");
      next if (!defined $status->[$j] || $status->[$j] eq "");
      push @$r, $j;
      $count++;
    }
    push @$res, [$g->[0], $g->[1], $r];
  }
  return $res;
}

sub plotSurvival {
  my ($self, $outfile, $groups, $ct, $maxt, $ph, $times, %params) = @_;
  my $time =   $self->getTimeData();
  my $status = $self->getStatusData();
  return &plotSurvivalR($outfile, $time, $status, $groups, $ct, $maxt, $ph, $times, %params);
}

sub plotSurvivalR {
  my ($outfile, $time, $status, $groups, $ct, $maxt, $ph, $times, %params) = @_;
  my $res = [];
  my $numgroups = 0;
  for (my $i = 0; $i < scalar(@$groups); $i++) {
    my $g = $groups->[$i];
    my $count = 0;
    foreach my $j (@{$g->[2]}) {
      next if (!defined $time->[$j] || $time->[$j] eq "");
      next if (!defined $status->[$j] || $status->[$j] eq "");
      $time->[$j] =~ s/\s//g;
      $status->[$j] =~ s/\s//g;
      $time->[$j] =~ s/NA//ig;
      $status->[$j] =~ s/NA//ig;
      next if (!defined $time->[$j] || $time->[$j] eq "");
      next if (!defined $status->[$j] || $status->[$j] eq "");
      push @$res, [$i+1, $time->[$j], $status->[$j]];
      $count++;
    }
    if ($count > 0) {
      $numgroups++;
    }
  }
  if ($numgroups <= 1) {
    return undef;
  }
  my ($phr, $phw);
  if (defined $ph) {
    ($phr, $phw) = @$ph;
  }
  else {
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  print $phw "library(survival)\n";
  print $phw "rm(list=ls())\n";
  print $phw "groups <- c(", join(",", map {$_->[0]} @$res), ")\n";
  print $phw "time <- c(", join(",", map {$_->[1]} @$res), ")\n";
  print $phw "status <- c(", join(",", map {$_->[2]} @$res), ")\n";
  print $phw "nm <- c(", join(",", map {"\"".$_->[0]."\""} @$groups), ")\n";
  print $phw "clr <- c(", join(",", map {"\"".$_->[1]."\""} @$groups), ")\n";

  my $devstr = "";
  if (defined $outfile) {
    my @l = split("\\.", $outfile);
    if ($l[$#l] eq "png") {
      $devstr = "png(filename=\"$outfile\", width=640, height=480, pointsize=15)";
    }
    if ($l[$#l] eq "ps") {
      $devstr = "postscript(file=\"$outfile\", width=6.4, height=4.8, pointsize=11)";
    }
    if ($l[$#l] eq "pdf") {
      $devstr = "pdf(file=\"$outfile\", width=6.4, height=4.8, pointsize=11)";
    }
  }

  print $phw "idx <- (!is.na(groups) & groups > 0)\n";

  if (defined $ct && $ct ne "") {
print $phw <<END
ct <- $ct
status[time > ct] <- 0;
time[time > ct] <- ct;
END
;
  }
  if (defined $maxt && $maxt ne "") {
    print $phw "maxt <- $maxt\n";
  }
  else {
    print $phw "maxt <- max(time)\n";
  }

print $phw <<END
time <- time[idx]
status <- status[idx]
groups <- groups[idx]

$devstr

s <- survfit(Surv(time, status) ~ groups)
l <- names(s\$strata)
x <- (max(s\$time) - min(s\$time)) * 0.7
par(mar=c(3, 3, 3, 3))
st <- survdiff(Surv(time, status) ~ groups)
g <- as.numeric(gsub("groups=", "", names(st\$n)))
nm <- nm[g]
clr <- clr[g]
p <- 1 - pchisq(st\$chisq, length(l)-1)
nm <- paste(nm, "(", st\$obs, "/", st\$n, ")",
    sprintf(rep("%.2f%%", length(st\$n)), s\$surv[cumsum(s\$strata)]*100))
#nm <- paste(nm, sprintf(rep("%.2f%%", length(st\$n)), st\$obs/st\$n*100), 
#    "(", st\$obs, "/", st\$n, ")")
t.ticks <- c(0)
s.ticks <- c(0)
END
;
  if ($devstr ne "") {
  print $phw <<END
plot(s, lty = 2:(length(l)+1), col=clr, lwd=3, xlim=c(0,maxt), xaxt="n")
t.ticks <- axTicks(1)
s.ticks <- axTicks(2)
END
;
  if (defined $times) {
    print $phw "t.ticks <- c(", join(",", @$times), ")\n";
  }
  if (defined $params{'ticks'}) {
    print $phw "t.ticks.l <- c(", join(",", @{$params{'ticks'}}), ")\n";
  }
  else {
    print $phw "t.ticks.l <- t.ticks\n";
  }
  print $phw <<END
axis(1, t.ticks, labels=t.ticks.l)
legend("topright", nm, lty = 2:(length(l)+1), col=clr, lwd=3, inset=0.02)
legend("topright", sprintf("pvalue = %.4f", p), inset=c(0.02, length(l)*0.042+0.12))
#text(90,0.8, labels=sprintf("p = %.4f", p), pos=2)
dev.off()
END
;
  }
  print $phw <<END
cat("pvalue", p, length(l), "\\n", sep="\\t")
cat("groups", nm, "\\n", sep="\\t")
cat("colors", clr, "\\n", sep="\\t")
cat("strata", s\$strata, "\\n", sep="\\t")
cat("time", s\$time, "\\n", sep="\\t")
cat("n.risk", s\$n.risk, "\\n", sep="\\t")
cat("n.event", s\$n.event, "\\n", sep="\\t")
cat("n.censor", s\$n.censor, "\\n", sep="\\t")
cat("surv", s\$surv, "\\n", sep="\\t")
cat("lower", s\$lower, "\\n", sep="\\t")
cat("upper", s\$upper, "\\n", sep="\\t")
cat("std.err", s\$std.err, "\\n", sep="\\t")
cat("t.ticks", t.ticks, "\\n", sep="\\t")
cat("s.ticks", s.ticks, "\\n", sep="\\t")

x <- coxph(Surv(time, status) ~ groups)
s <- summary(x)
codes <- c("0", "***", "**", "*", ".", " ")
p <- s\$coefficients[,5]
pg <- as.numeric(p <= 0) + 2 * as.numeric(p > 0 & p <= 0.001) +
    3 * as.numeric(p > 0.001 & p <= 0.01) +
    4 * as.numeric(p > 0.01 & p <= 0.05) +
    5 * as.numeric(p > 0.05 & p <= 0.1) +
    6 * as.numeric(p > 0.1)
cat("names", rownames(s\$conf.int), "\\n", sep="\\t")
cat("exp(coef)", s\$conf.int[1], "\\n", sep="\\t")
cat("lower .95", s\$conf.int[3], "\\n", sep="\\t")
cat("upper .95", s\$conf.int[4], "\\n", sep="\\t")
cat("cp", p, "\\n", sep="\\t")
cat("c", codes[pg], "\\n", sep="\\t")
cat("end", "\\n",  sep="")
END
;
  if ($devstr ne "") {
    my $h = <$phr>;
    my $h = <$phr>;
  }
  my $hash = {};
  while (my $h = <$phr>) {
    $h =~ s/[\r\n]//g;
    if ($h eq "end") {
      last;
    }
    my ($k, $val) = split("\t", $h, 2);
    push @{$hash->{$k}}, $val;
  }
  if (!defined $ph) {
    close($phr);
    close($phw);
  }
  return $hash;
}

sub printSurvStats {
  my ($res, $times) = @_;
  &printUMStats($res);
  print "pvalue ", join("\t", @{$res->{"pvalue"}}), "\n";
  print "groups ", join("\t", @{$res->{"groups"}}), "\n";
  #print "strata ", join("\t", @{$res->{"strata"}}), "\n";
  #print "time ", join("\t", @{$res->{"time"}}), "\n";
  #print "risk ", join("\t", @{$res->{"n.risk"}}), "\n";
  #print "surv ", join("\t", @{$res->{"surv"}}), "\n";
  if (!defined $times) {
    $times = [split("\t", $res->{"t.ticks"}->[0])];
  }
  my @groups = split("\t", $res->{"groups"}->[0]);
  my @colors = split("\t", $res->{"colors"}->[0]);
  my @strata = split("\t", $res->{"strata"}->[0]);
  my @time = split("\t", $res->{"time"}->[0]);
  my @risk = split("\t", $res->{"n.risk"}->[0]);
  #my @event = split("\t", $res->{"n.event"}->[0]);
  #my @censor = split("\t", $res->{"n.censor"}->[0]);
  my @surv = split("\t", $res->{"surv"}->[0]);
  my @upper = split("\t", $res->{"upper"}->[0]);
  my @lower = split("\t", $res->{"lower"}->[0]);
  my @stderr = split("\t", $res->{"std.err"}->[0]);
  #print "Time\tRisk\tEvent\tCensor\tSurv\n";
  #for (my $i = 0; $i <= scalar(@risk); $i++) {
  #  printf "\%d\t\%d\t\%d\t\%d\t%.2f\n", $time[$i], 
  #         $risk[$i], $event[$i], $censor[$i], $surv[$i];
  #}
  my @stimes = sort { $a <=> $b } @$times;
  print "No. at Risk\t", join("\t", @stimes), "\n";
  my $index = 0;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $start = $index;
    my $limit = $index + $strata[$g];
    my $gr = $groups[$g];
    $gr =~ s/ \(.*//g;
    print $gr;
    foreach my $t (@stimes) {
      while ($index < $limit && $time[$index] < $t) { $index++ }
      if ($index >= $limit) {
        print "\t", 0;
      }
      else {
        print "\t", $risk[$index];
      }
    }
    print "\n";
    $index = $limit;
  }
  print "Surv prob\t", join("\t", @stimes), "\n";
  my $index = 0;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $start = $index;
    my $limit = $index + $strata[$g];
    my $gr = $groups[$g];
    $gr =~ s/ \(.*//g;
    print $gr;
    foreach my $t (@stimes) {
      while ($index < $limit && $time[$index] < $t) { $index++ }
      if ($index >= $limit) {
        #printf "\t%.2f", $surv[$limit - 1];
        #printf "\t%.2f (%.2f - %.2f)", $surv[$limit - 1],
        #$lower[$limit - 1], $upper[$limit - 1];
        printf "\t%.1f(%.1f-%.1f)%.1f", $surv[$limit - 1] * 100,
        $lower[$limit - 1] * 100,
        $upper[$limit - 1] * 100,
        $stderr[$limit - 1] * $surv[$limit - 1] * 100;
      }
      else {
        #printf "\t%.2f", $surv[$index];
        #printf "\t%.2f (%.2f - %.2f)", $surv[$index],
        #$lower[$index], $upper[$index];
        printf "\t%.1f(%.1f-%.1f)%.1f", $surv[$index] * 100,
        $lower[$index] * 100,
        $upper[$index] * 100,
        $stderr[$index] * $surv[$index] * 100;
      }
    }
    print "\n";
    $index = $limit;
  }
}

sub getPScolor {
  my ($n, $c) = @_;
  if ($c =~ /#/) {
    $c =~ s/^(\#|Ox)//;
    $_ = $c;
    my ($r, $g, $b) = m/(\w{2})(\w{2})(\w{2})/;
    my ($rf, $gf, $bf) = map { hex($_) } ($r, $g, $b);
    return "\\definecolor{$n}{RGB}{$rf,$gf,$bf}";
  }
  return "\\definecolor{$n}{named}{$c}";
}

sub printSurvTikz {
  my ($ofile, $res, $times, %params) = @_;
  open(my $ofh, ">$ofile") || die "Can't write $ofile\n";
  if (!defined $times) {
    $times = [split("\t", $res->{"t.ticks"}->[0])];
  }
  #print "groups ", join("\t", @{$res->{"groups"}}), "\n";
  #print "strata ", join("\t", @{$res->{"strata"}}), "\n";
  my @groups = split("\t", $res->{"groups"}->[0]);
  my @colors = split("\t", $res->{"colors"}->[0]);
  my @strata = split("\t", $res->{"strata"}->[0]);
  my @time = split("\t", $res->{"time"}->[0]);
  my @risk = split("\t", $res->{"n.risk"}->[0]);
  my @event = split("\t", $res->{"n.event"}->[0]);
  my @censor = split("\t", $res->{"n.censor"}->[0]);
  my @surv = split("\t", $res->{"surv"}->[0]);
  #print "Time\tRisk\tEvent\tCensor\tSurv\n";
  #for (my $i = 0; $i <= scalar(@risk); $i++) {
  #  printf "$i\t\%d\t\%d\t\%d\t\%d\t%.2f\n", $time[$i], 
  #         $risk[$i], $event[$i], $censor[$i], $surv[$i];
  #}
  my $index = 0;
  my $cnum = [];
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $start = $index;
    my $limit = $index + $strata[$g];
    my $c = chr(ord('a')+$g);
    print $ofh "\\pgfplotstableread\{%\n";
    my $s = 1.0;
    my $t = 0;
    print $ofh join(" ", $t, $s), "\n";
    for (my $i = $start; $i < $limit; $i++) {
      print $ofh join(" ", $time[$i], $surv[$i]), "\n";
      ($t, $s) = ($time[$i], $surv[$i]);
    }
    print $ofh "}\\g$c\data%\n";
    my $cname = "\\g$c"."censor";
    my $str = "\\pgfplotstableread\{%\n";
    my $n = 0;
    for (my $i = $start; $i < $limit; $i++) {
      if ($censor[$i] > 0) {
        $str = $str . join(" ", $time[$i], $surv[$i]). "\n";
        $n++;
      }
    }
    $cnum->[$g] = $n;
    if ($n > 0) {
      print $ofh "$str}$cname%\n";
    }
    $index = $limit;
  }
  my @stimes = sort { $a <=> $b } @$times;
  my $maxtime = $stimes[$#stimes];
  my $xmax = $maxtime * 1.1;
  my $xunit = 1.5/$maxtime;
  my $yunit = 1.5;
  my $middle = $maxtime/2;
  my $leftx = 0.3/$xunit;
  my $bottom = 0.4/$yunit;
  my $l1 = 0.75/$xunit;
  my $dx = $stimes[1];
  my $rs = 0.08;
  my $boty = $bottom + (scalar(@strata)+1) * $rs;
  my $pvalue = $res->{"pvalue"}->[0];
  $pvalue =~ s/\s.*//g;
  $pvalue = sprintf("%.2g", $pvalue);
  $params{'ylabel'} = 'Overall Survival' if (!defined $params{'ylabel'});
  $params{'xlabel'} = 'Months' if (!defined $params{'xlabel'});
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $clr = &getPScolor("clr$g", $colors[$g]);
    print $ofh "$clr\n";
  }
  print $ofh <<END
\\def\\ttsize{\\tiny}%
\\def\\ttsizea{\\fontsize{9pt}{9pt}\\selectfont}%
\\def\\ttsizeb{\\fontsize{8pt}{8pt}\\selectfont}%
\\def\\ttsizec{\\fontsize{7pt}{7pt}\\selectfont}%
\\def\\pshlabel#1{\\ttsize #1}%
\\def\\psvlabel#1{\\ttsize #1}%
\\begin{tikzpicture}[x=$xunit\in, y=$yunit\in]
\\node[anchor=west]  at ($leftx,0){\\bf \\ttsizea p = $pvalue};
\\node[anchor=west] at (-$l1,-$bottom){\\bf \\ttsizec No at risk};
\\draw[black,line width=0.75pt] (0, -0.1) rectangle ($xmax, 1.1);
\\begin{axis}[x=$xunit\in, y=$yunit\in,
  anchor=origin, axis x line = bottom,axis y line = left,
  ticklabel style={font=\\tiny},
  ymin=-0.1, ymax=1.1, xmin=0, xmax=$xmax,
  xlabel={\\bf \\ttsizeb $params{'xlabel'}},
  ylabel={\\bf \\ttsizea $params{'ylabel'}},
  xlabel style={below,anchor=north},
  ylabel style={left,anchor=south},
  xtick={0,$dx,...,$maxtime},
  ytick={0,0.2,...,1.0}]
END
;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $c = chr(ord('a')+$g);
    my $cname = "\\g$c"."censor";
    print $ofh "\\addplot+[line width=1.3pt, color=clr$g, const plot, no marks]  table [x index=0, y index=1] {\\g$c\data};\n";
    if ($cnum->[$g] > 0) {
      print $ofh "\\addplot+[line width=0.75pt, color=clr$g, mark=+, only marks]  table [x index=0, y index=1] {$cname};\n";
    }
  }
  print $ofh "\\end{axis}\n";
  my $nrisk = [["No. at Risk", [@stimes], []]];
  my $index = 0;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $start = $index;
    my $limit = $index + $strata[$g];
    my $gr = $groups[$g];
    $gr =~ s/ \(.*//g;
    my @numbers;
    my @snumbers;
    foreach my $t (@stimes) {
      while ($index < $limit && $time[$index] < $t) { $index++ }
      if ($index >= $limit) {
        push @numbers, 0;
        push @snumbers, $surv[$limit - 1];
      }
      else {
        push @numbers, $risk[$index];
        push @snumbers, $surv[$index];
      }
    }
    push @$nrisk, [$gr, [@numbers], [@snumbers]];
    $index = $limit;
  }
  for (my $i = 1; $i < scalar(@$nrisk); $i++) {
    my $g = $i - 1;
    my $y = $bottom + $i * $rs;
    my $ly = 1.0 - $g * $rs;
    my $gr = $nrisk->[$i]->[0];
    my $pr = sprintf("\%.1f\\%%", $nrisk->[$i]->[2]->[$#stimes] * 100);
    print $ofh "\\node[anchor=east,color=clr$g] at ($maxtime,$ly){\\bf \\ttsizeb $gr ($pr)};\n";
    print $ofh "\\node[anchor=west] at (-$l1,-$y){\\bf \\ttsizec $gr};\n";
  }
  for (my $j = 0; $j < scalar(@{$nrisk->[0]->[1]}); $j++) {
    for (my $i = 1; $i < scalar(@$nrisk); $i++) {
      my $y = $bottom + $i * $rs;
      my $t = $nrisk->[0]->[1]->[$j];
      my $n = $nrisk->[$i]->[1]->[$j];
      print $ofh "\\node[anchor=west,inner sep=0pt] at ($t,-$y) {\\bf \\ttsizec $n};\n";
    } 
  }
  print $ofh "\\end{tikzpicture}\n";
  close($ofh);
}

sub printSurvPstricks {
  my ($ofile, $res, $times, %params) = @_;
  open(my $ofh, ">$ofile") || die "Can't write $ofile\n";
  if (!defined $times) {
    $times = [split("\t", $res->{"t.ticks"}->[0])];
  }
  #print "groups ", join("\t", @{$res->{"groups"}}), "\n";
  #print "strata ", join("\t", @{$res->{"strata"}}), "\n";
  my @groups = split("\t", $res->{"groups"}->[0]);
  my @colors = split("\t", $res->{"colors"}->[0]);
  my @strata = split("\t", $res->{"strata"}->[0]);
  my @time = split("\t", $res->{"time"}->[0]);
  my @risk = split("\t", $res->{"n.risk"}->[0]);
  my @event = split("\t", $res->{"n.event"}->[0]);
  my @censor = split("\t", $res->{"n.censor"}->[0]);
  my @surv = split("\t", $res->{"surv"}->[0]);
  #print "Time\tRisk\tEvent\tCensor\tSurv\n";
  #for (my $i = 0; $i <= scalar(@risk); $i++) {
  #  printf "$i\t\%d\t\%d\t\%d\t\%d\t%.2f\n", $time[$i], 
  #         $risk[$i], $event[$i], $censor[$i], $surv[$i];
  #}
  my $index = 0;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $start = $index;
    my $limit = $index + $strata[$g];
    my $c = chr(ord('a')+$g);
    print $ofh "\\def\\g$c\data{%\n";
    my $s = 1.0;
    my $t = 0;
    print $ofh join(" ", $t, $s), " ";
    for (my $i = $start; $i < $limit; $i++) {
      print $ofh join(" ", $time[$i], $s, $time[$i], $surv[$i]), " ";
      ($t, $s) = ($time[$i], $surv[$i]);
    }
    print $ofh "}%\n";
    my $cname = "\\g$c"."censor";
    print $ofh "\\def$cname\{%\n";
    for (my $i = $start; $i < $limit; $i++) {
      if ($censor[$i] > 0) {
        print $ofh join(" ", $time[$i], $surv[$i]), " ";
      }
    }
    print $ofh "}%\n";
    $index = $limit;
  }
  my @stimes = sort { $a <=> $b } @$times;
  my $maxtime = $stimes[$#stimes];
  my $xmax = $maxtime * 1.1;
  my $xunit = 1.5/$maxtime;
  my $yunit = 1.5;
  my $middle = $maxtime/2;
  my $leftx = 0.3/$xunit;
  my $bottom = 0.4/$yunit;
  my $l1 = 0.5/$xunit;
  my $dx = $stimes[1];
  my $rs = 0.08;
  my $boty = $bottom + (scalar(@strata)+1) * $rs;
  my $pvalue = $res->{"pvalue"}->[0];
  $pvalue =~ s/\s.*//g;
  $pvalue = sprintf("%.2g", $pvalue);
  $params{'ylabel'} = 'Overall Survival' if (!defined $params{'ylabel'});
  $params{'xlabel'} = 'Months' if (!defined $params{'xlabel'});
  print $ofh <<END
\\def\\ttsize{\\tiny}%
\\def\\ttsizea{\\fontsize{9pt}{9pt}\\selectfont}%
\\def\\ttsizeb{\\fontsize{8pt}{8pt}\\selectfont}%
\\def\\ttsizec{\\fontsize{7pt}{7pt}\\selectfont}%
\\def\\pshlabel#1{\\ttsize #1}%
\\def\\psvlabel#1{\\ttsize #1}%
\\psset{xunit=$xunit\in,yunit=$yunit\in}%
\\begin{pspicture}(-$l1,-$boty)($xmax,1.1)
\\psframe[dimen=outer, linewidth=0.75pt, linecolor=black](0, -0.1)($xmax, 1.1)
\\psaxes[xAxis=false,subticks=5,ticksize=-2pt 0,linewidth=0.5pt,linecolor=black,Dy=0.2]{-}($maxtime,1)
\\psaxes[yAxis=false,subticks=5,ticksize=-2pt 0,linewidth=0.5pt,linecolor=black,Dx=$dx]{-}(0,-0.1)(0,-0.1)($maxtime,1)
\%\\psaxes[yAxis=false,subticks=5,ticksize=-2pt 0,linewidth=0.5pt,linecolor=black,Ox=0,Dx=1,dx=$dx]{-}(0,-0.1)(0,-0.1)($maxtime,1)
\\rput[c]($leftx,0){\\bf \\ttsizea p = $pvalue}
\\rput[c]{90}(-$leftx,0.5){\\bf \\ttsizea $params{'ylabel'}}
\\rput[B]($middle,-$bottom){\\bf \\ttsizeb $params{'xlabel'}}
\\rput[Bl](-$l1,-$bottom){\\bf \\ttsizec No at risk}
\\makeatletter
\\\@namedef{psds\@cd}{% 
/DS DS 1.253 mul def
\\pst\@gdot{DS 0 moveto DS neg 0 L stroke 0 DS moveto 0 DS neg L stroke}}
\\def\\psdots\@iii{%
 \\psk\@dotsize
 \\\@nameuse{psds\@\\psk\@dotstyle}
 newpath
 n { transform floor .5 add exch floor
     .5 add exch itransform Dot stroke} repeat }
\\makeatother
END
;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $clr = &getPScolor("clr$g", $colors[$g]);
    my $y = $bottom + ($g + 1) * $rs;
    my $gr = $groups[$g];
    $gr =~ s/ \(.*//g;
    my $c = chr(ord('a')+$g);
    my $cname = "\\g$c"."censor";
    print $ofh <<END
$clr
\\listplot[linecolor=clr$g,linewidth=1.3pt]{\\g$c\data}
\\listplot[linecolor=clr$g,linewidth=0.75pt,plotstyle=dots,dotstyle=cd,dotsize=3pt]{$cname}
\\rput[Bl](-$l1,-$y){\\bf \\ttsizec $gr}
END
;
  }
  my $nrisk = [["No. at Risk", [@stimes], []]];
  my $index = 0;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $start = $index;
    my $limit = $index + $strata[$g];
    my $gr = $groups[$g];
    $gr =~ s/ \(.*//g;
    my @numbers;
    my @snumbers;
    foreach my $t (@stimes) {
      while ($index < $limit && $time[$index] < $t) { $index++ }
      if ($index >= $limit) {
        push @numbers, 0;
        push @snumbers, $surv[$limit - 1];
      }
      else {
        push @numbers, $risk[$index];
        push @snumbers, $surv[$index];
      }
    }
    push @$nrisk, [$gr, [@numbers], [@snumbers]];
    $index = $limit;
  }
  for (my $i = 1; $i < scalar(@$nrisk); $i++) {
    my $g = $i - 1;
    my $ly = 1.0 - $g * $rs;
    my $gr = $nrisk->[$i]->[0];
    my $pr = sprintf("\%.1f\\%%", $nrisk->[$i]->[2]->[$#stimes] * 100);
    print $ofh "\\rput[tr]($maxtime,$ly){\\bf \\ttsizeb \\color{clr$g}$gr ($pr)}\n";
  }
  for (my $j = 0; $j < scalar(@{$nrisk->[0]->[1]}); $j++) {
    for (my $i = 1; $i < scalar(@$nrisk); $i++) {
      my $y = $bottom + $i * $rs;
      my $t = $nrisk->[0]->[1]->[$j];
      my $n = $nrisk->[$i]->[1]->[$j];
      print $ofh "\\rput[B]($t,-$y){\\bf \\ttsizec $n}\n";
    } 
  }
  print $ofh "\\end{pspicture}\n";
  close($ofh);
}

sub printSurvPstricksOld {
  my ($ofile, $res, $times, %params) = @_;
  open(my $ofh, ">$ofile") || die "Can't write $ofile\n";
  if (!defined $times) {
    $times = [split("\t", $res->{"t.ticks"}->[0])];
  }
  #print "groups ", join("\t", @{$res->{"groups"}}), "\n";
  #print "strata ", join("\t", @{$res->{"strata"}}), "\n";
  my @groups = split("\t", $res->{"groups"}->[0]);
  my @colors = split("\t", $res->{"colors"}->[0]);
  my @strata = split("\t", $res->{"strata"}->[0]);
  my @time = split("\t", $res->{"time"}->[0]);
  my @risk = split("\t", $res->{"n.risk"}->[0]);
  my @event = split("\t", $res->{"n.event"}->[0]);
  my @censor = split("\t", $res->{"n.censor"}->[0]);
  my @surv = split("\t", $res->{"surv"}->[0]);
  #print "Time\tRisk\tEvent\tCensor\tSurv\n";
  #for (my $i = 0; $i <= scalar(@risk); $i++) {
  #  printf "$i\t\%d\t\%d\t\%d\t\%d\t%.2f\n", $time[$i], 
  #         $risk[$i], $event[$i], $censor[$i], $surv[$i];
  #}
  my $index = 0;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $start = $index;
    my $limit = $index + $strata[$g];
    my $c = chr(ord('a')+$g);
    print $ofh "\\def\\g$c\data{%\n";
    my $s = 1.0;
    my $t = 0;
    print $ofh join(" ", $t, $s), " ";
    for (my $i = $start; $i < $limit; $i++) {
      print $ofh join(" ", $time[$i], $s, $time[$i], $surv[$i]), " ";
      ($t, $s) = ($time[$i], $surv[$i]);
    }
    print $ofh "}%\n";
    my $cname = "\\g$c"."censor";
    print $ofh "\\def$cname\{%\n";
    for (my $i = $start; $i < $limit; $i++) {
      if ($censor[$i] > 0) {
        print $ofh join(" ", $time[$i], $surv[$i]), " ";
      }
    }
    print $ofh "}%\n";
    $index = $limit;
  }
  my @stimes = sort { $a <=> $b } @$times;
  my $maxtime = $stimes[$#stimes];
  my $xmax = $maxtime * 1.1;
  my $xunit = 1.5/$maxtime;
  my $yunit = 1.5;
  my $middle = $maxtime/2;
  my $leftx = 0.3/$xunit;
  my $bottom = 0.4/$yunit;
  my $l1 = 0.5/$xunit;
  my $dx = $stimes[1];
  my $rs = 0.08;
  my $boty = $bottom + (scalar(@strata)+1) * $rs;
  my $pvalue = $res->{"pvalue"}->[0];
  $pvalue =~ s/\s.*//g;
  $pvalue = sprintf("%.2g", $pvalue);
  $params{'ylabel'} = 'Overall Survival' if (!defined $params{'ylabel'});
  $params{'xlabel'} = 'Months' if (!defined $params{'xlabel'});
  print $ofh <<END
\\def\\ttsize{\\tiny}%
\\def\\pshlabel#1{\\ttsize #1}%
\\def\\psvlabel#1{\\ttsize #1}%
\\psset{xunit=$xunit\in,yunit=$yunit\in}%
\\begin{pspicture}(-$l1,-$boty)($xmax,1.1)
\\psframe[dimen=outer, linewidth=0.5pt, linecolor=black](0, -0.1)($xmax, 1.1)
\\psaxes[xAxis=false,subticks=5,ticksize=-2pt 0,linewidth=0.5pt,linecolor=black,Dy=0.2]{-}($maxtime,1)
\\psaxes[yAxis=false,subticks=5,ticksize=-2pt 0,linewidth=0.5pt,linecolor=black,Dx=$dx]{-}(0,-0.1)(0,-0.1)($maxtime,1)
\%\\psaxes[yAxis=false,subticks=5,ticksize=-2pt 0,linewidth=0.5pt,linecolor=black,Ox=0,Dx=1,dx=$dx]{-}(0,-0.1)(0,-0.1)($maxtime,1)
\\rput[c]($leftx,0){\\bf \\ttsize p = $pvalue}
\\rput[c]{90}(-$leftx,0.5){\\bf \\ttsize $params{'ylabel'}}
\\rput[B]($middle,-$bottom){\\bf \\ttsize $params{'xlabel'}}
\\rput[Bl](-$l1,-$bottom){\\bf \\ttsize No at risk}
\\makeatletter
\\\@namedef{psds\@cd}{% 
/DS DS 1.253 mul def
\\pst\@gdot{DS 0 moveto DS neg 0 L stroke 0 DS moveto 0 DS neg L stroke}}
\\def\\psdots\@iii{%
 \\psk\@dotsize
 \\\@nameuse{psds\@\\psk\@dotstyle}
 newpath
 n { transform floor .5 add exch floor
     .5 add exch itransform Dot stroke} repeat }
\\makeatother
END
;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $clr = &getPScolor("clr$g", $colors[$g]);
    my $y = $bottom + ($g + 1) * $rs;
    my $gr = $groups[$g];
    $gr =~ s/ \(.*//g;
    my $c = chr(ord('a')+$g);
    my $cname = "\\g$c"."censor";
    print $ofh <<END
$clr
\\listplot[linecolor=clr$g,linewidth=1pt]{\\g$c\data}
\\listplot[linecolor=clr$g,linewidth=0.5pt,plotstyle=dots,dotstyle=cd,dotsize=3pt]{$cname}
\\rput[Bl](-$l1,-$y){\\bf \\ttsize $gr}
END
;
  }
  my $nrisk = [["No. at Risk", [@stimes], []]];
  my $index = 0;
  for (my $g = 0; $g < scalar(@strata); $g++) {
    my $start = $index;
    my $limit = $index + $strata[$g];
    my $gr = $groups[$g];
    $gr =~ s/ \(.*//g;
    my @numbers;
    my @snumbers;
    foreach my $t (@stimes) {
      while ($index < $limit && $time[$index] < $t) { $index++ }
      if ($index >= $limit) {
        push @numbers, 0;
        push @snumbers, $surv[$limit - 1];
      }
      else {
        push @numbers, $risk[$index];
        push @snumbers, $surv[$index];
      }
    }
    push @$nrisk, [$gr, [@numbers], [@snumbers]];
    $index = $limit;
  }
  for (my $i = 1; $i < scalar(@$nrisk); $i++) {
    my $g = $i - 1;
    my $ly = 1.0 - $g * $rs;
    my $gr = $nrisk->[$i]->[0];
    my $pr = sprintf("\%.1f\\%%", $nrisk->[$i]->[2]->[$#stimes] * 100);
    print $ofh "\\rput[tr]($maxtime,$ly){\\bf \\ttsize \\color{clr$g}$gr ($pr)}\n";
  }
  for (my $j = 0; $j < scalar(@{$nrisk->[0]->[1]}); $j++) {
    for (my $i = 1; $i < scalar(@$nrisk); $i++) {
      my $y = $bottom + $i * $rs;
      my $t = $nrisk->[0]->[1]->[$j];
      my $n = $nrisk->[$i]->[1]->[$j];
      print $ofh "\\rput[B]($t,-$y){\\bf \\ttsize $n}\n";
    } 
  }
  print $ofh "\\end{pspicture}\n";
  close($ofh);
}

sub plotSurvivalThr {
  my ($self, $outfile, $id, $thr, $ct, $maxt, $ph) = @_;
  my $hi = $self->getArraysThr($id, $thr, "hi");
  my $lo = $self->getArraysThr($id, $thr, "lo");
  if (!defined $hi || !defined $lo || scalar(@$hi) == 0 || scalar(@$lo) == 0) {
    return undef;
  }
  my $name = $self->getName($id);
  my $thr_step = $self->getThrData($id);
  $thr = &getThrCode($thr_step, $thr_step->[0], $thr);
  my $groups = [ ["$name low ", "green", $lo],
                 ["$name high", "red", $hi] ];
  my $res = $self->plotSurvival($outfile, $groups, $ct, $maxt, $ph);
  $res->{'id'} = [$id];
  $res->{'thr'} = [$thr];
  return $res;
}

sub writeIndexExpr {
  my ($ofile, $file) = @_;
  open(my $ofh, ">$ofile") || die "Can't write $ofile\n";
  open(my $fh, "<$file") || die "Can't open $file\n";
  my $head = <$fh>;
  print $ofh "ProbeID\tPtr\tName\tDescription\n";
  my $ptr = tell($fh);
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my ($name, $desc) = split(/: /, $list[1]);
    print $ofh "$list[0]\t$ptr\t$name\t$desc\n";
    $ptr = tell($fh);
  }
  close($fh);
  close($ofh);
}

sub correlation {
  my ($self, $id, $arrays) = @_;
  my $expr = $self->getExprData($id);
  my @expr1 = map { $expr->[$_] } @$arrays;
  my $fh = $self->{'fh'};
  seek($fh, 0, 0);
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @header = split("\t", $head);
  my $res = [];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my @expr2 = map { $list[$_] } @$arrays;
    my $corr = &U::correlation(\@expr1, \@expr2);
    push @$res, [$corr, $list[0], $list[1]];
  }
  return $res;
}

sub convertValues {
  my ($arr, %params) = @_;
  my $res = [];
  for(my $i = 0; $i < scalar(@$arr); $i++) {
    if (defined $params{$arr->[$i]}) {
      $res->[$i] = $params{$arr->[$i]};
    }
    else {
      $res->[$i] = "";
    }
  }
  return $res;
}

sub convertValuesOnly {
  my ($arr, %params) = @_;
  my $res = [];
  for(my $i = 0; $i < scalar(@$arr); $i++) {
    if (defined $params{$arr->[$i]}) {
      $res->[$i] = $params{$arr->[$i]};
    }
    else {
      $res->[$i] = $arr->[$i];
    }
  }
  return $res;
}

sub countValues {
  my ($self, $arr, $rest) = @_;
  my $res = {};
  for(my $i = $self->getStart(); $i <= $self->getEnd(); $i++) {
    my $k = "Other";
    if (defined $arr->[$i]) {
      $k = $arr->[$i];
    }
    if (!defined $res->{$k}) {
      $res->{$k} = 0;
    }
    $res->{$k}++;
  }
  return $res;
}

sub countTable {
  my ($self, $arr1, $arr2) = @_;
  my $res = {};
  for(my $i = $self->getStart(); $i <= $self->getEnd(); $i++) {
    my $k1 = "Other";
    if (defined $arr1->[$i]) {
      $k1 = $arr1->[$i];
    }
    my $k2 = "Other";
    if (defined $arr2->[$i]) {
      $k2 = $arr2->[$i];
    }
    if (!defined $res->{$k1}) {
      $res->{$k1} = {};
    }
    if (!defined $res->{$k1}->{$k2}) {
      $res->{$k1}->{$k2} = 0;
    }
    $res->{$k1}->{$k2}++;
  }
  return $res;
}

sub meanTable {
  my ($self, $arr1, $arr2) = @_;
  my $res = {};
  for(my $i = $self->getStart(); $i <= $self->getEnd(); $i++) {
    my $k1 = "Other";
    if (defined $arr1->[$i]) {
      $k1 = $arr1->[$i];
    }
    my $k2 = "";
    if (defined $arr2->[$i]) {
      $k2 = $arr2->[$i];
    }
    if (!defined $res->{$k1}) {
      $res->{$k1} = [];
    }
    if ($k2 ne "") {
      push @{$res->{$k1}}, $k2;
    }
  }
  my $m = {};
  my $c = {};
  foreach my $k (keys %{$res}) {
    $m->{$k} = &U::mean($res->{$k});
    $c->{$k} = scalar(@{$res->{$k}});
  }
  return [$m, $c];
}

sub countGTable {
  my ($self, $arr1, $groups) = @_;
  my $start = $self->getStart();
  my $end = $self->getEnd();
  my $hash = {};
  my $gr = [];
  foreach my $g (@$groups) {
    $hash->{$g->[0]} = {};
    foreach my $i (@{$g->[2]}) {
      $hash->{$g->[0]}->{$i} = 1;
      $gr->[$i] = 1;
    }
  }
  $hash->{"Other"} = {};
  for(my $i = $start; $i <= $end; $i++) {
    next if (defined $gr->[$i]);
    $hash->{"Other"}->{$i} = 1;
  }
  my $res = {};
  for(my $i = $start; $i <= $end; $i++) {
    my $k1 = "Other";
    if (defined $arr1->[$i]) {
      $k1 = $arr1->[$i];
    }
    if (!defined $res->{$k1}) {
      $res->{$k1} = {};
    }
    foreach my $k2 (%{$hash}) {
      next if (!defined $hash->{$k2}->{$i});
      if (!defined $res->{$k1}->{$k2}) {
        $res->{$k1}->{$k2} = 0;
      }
      $res->{$k1}->{$k2}++;
    }
  }
  return $res;
}

sub printGTable {
  my ($self, $arr1, $groups) = @_;
  my $count = $self->countGTable($arr1, $groups);
  my @keys = ((map { $_->[0] } @$groups), "Other");
  print join("\t", "Keys", @keys, "Total"), "\n";
  foreach my $k (keys %{$count}) {
    my $a = [map { $count->{$k}->{$_} } @keys];
    print "$k\t", join("\t", @$a, &U::sum($a)), "\n";
  }
}

sub multivariate {
  my ($self, $ct, %params) = @_;
  my $time = $self->getTimeData();
  my $status = $self->getStatusData();
  return &multivariateR($time, $status, $ct, %params);
}

sub multivariateR {
  my ($time, $status, $ct, %params) = @_;
  my $res = [];
  for(my $i = 2; $i < scalar(@$time); $i++) {
    next if (!defined $time->[$i] || $time->[$i] eq "");
    next if (!defined $status->[$i] || $status->[$i] eq "");
    my $empty = 0;
    foreach my $k (keys(%params)) {
      if (!defined $params{$k}->[$i] || $params{$k}->[$i] eq "") {
        $empty = 1;
        last;
      }
    }
    next if ($empty == 1);
    push @$res, $i;
  }
  print "Total Arrays ", scalar(@$res), "\n";

  my ($phr, $phw);
  my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  #$phw = \*STDOUT;
  print $phw "library(survival)\n";
  print $phw "rm(list=ls())\n";
  print $phw "time <- c(", join(",", map {$time->[$_]} @$res), ")\n";
  print $phw "status <- c(", join(",", map {$status->[$_]} @$res), ")\n";
  foreach my $k (keys(%params)) {
    print $phw "$k <- c(", join(",", map {$params{$k}->[$_]} @$res), ")\n";
  }
  if (defined $ct && $ct ne "") {
print $phw <<END
ct <- $ct
status[time > ct] <- 0;
time[time > ct] <- ct;
END
;
  }
  my $str = "";
  foreach my $k (keys(%params)) {
    if ($str eq "") { $str = $k; }
    else { $str = "$str+$k"; }
  }
print $phw <<END
x <- coxph(Surv(time, status) ~ $str)
s <- summary(x)
codes <- c("0", "***", "**", "*", ".", " ")
p <- s\$coefficients[,5]
pg <- as.numeric(p <= 0) + 2 * as.numeric(p > 0 & p <= 0.001) +
    3 * as.numeric(p > 0.001 & p <= 0.01) +
    4 * as.numeric(p > 0.01 & p <= 0.05) +
    5 * as.numeric(p > 0.05 & p <= 0.1) +
    6 * as.numeric(p > 0.1)
END
;
  if (scalar(keys(%params)) == 1) {
print $phw <<END
cat("names", rownames(s\$conf.int), "\\n", sep="\\t")
cat("exp(coef)", s\$conf.int[1], "\\n", sep="\\t")
cat("lower .95", s\$conf.int[3], "\\n", sep="\\t")
cat("upper .95", s\$conf.int[4], "\\n", sep="\\t")
cat("cp", p, "\\n", sep="\\t")
cat("c", codes[pg], "\\n", sep="\\t")
END
;
  }
  else {
print $phw <<END
cat("names", rownames(s\$conf.int), "\\n", sep="\\t")
cat("exp(coef)", s\$conf.int[,1], "\\n", sep="\\t")
cat("lower .95", s\$conf.int[,3], "\\n", sep="\\t")
cat("upper .95", s\$conf.int[,4], "\\n", sep="\\t")
cat("cp", p, "\\n", sep="\\t")
cat("c", codes[pg], "\\n", sep="\\t")
END
;
  }
  my $hash = {};
  for (my $i = 0; $i < 6; $i++) {
    my $h = <$phr>;
    $h =~ s/[\r\n]//g;
    my ($k, @list) = split("\t", $h);
    $hash->{$k} = [@list];
  }
  close($phr);
  close($phw);
  return $hash;
}

sub univariate {
  my ($self, $ct, %params) = @_;
  my $time = $self->getTimeData();
  my $status = $self->getStatusData();
  return &univariateR($time, $status, $ct, %params);
}

sub univariateR {
  my ($time, $status, $ct, %params) = @_;
  my $res = [];
  for(my $i = 2; $i < scalar(@$time); $i++) {
    next if (!defined $time->[$i] || $time->[$i] eq "");
    next if (!defined $status->[$i] || $status->[$i] eq "");
    my $empty = 0;
    foreach my $k (keys(%params)) {
      if (!defined $params{$k}->[$i] || $params{$k}->[$i] eq "") {
        $empty = 1;
        last;
      }
    }
    next if ($empty == 1);
    push @$res, $i;
  }
  print "Total Arrays ", scalar(@$res), "\n";

  my ($phr, $phw);
  my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  #$phw = \*STDOUT;
  print $phw "library(survival)\n";
  print $phw "rm(list=ls())\n";
  print $phw "time <- c(", join(",", map {$time->[$_]} @$res), ")\n";
  print $phw "status <- c(", join(",", map {$status->[$_]} @$res), ")\n";
  foreach my $k (keys(%params)) {
    print $phw "$k <- c(", join(",", map {$params{$k}->[$_]} @$res), ")\n";
  }
  if (defined $ct && $ct ne "") {
print $phw <<END
ct <- $ct
status[time > ct] <- 0;
time[time > ct] <- ct;
END
;
  }
  foreach my $str (keys(%params)) {
print $phw <<END
x <- coxph(Surv(time, status) ~ $str)
s <- summary(x)
codes <- c("0", "***", "**", "*", ".", " ")
p <- s\$coefficients[,5]
pg <- as.numeric(p <= 0) + 2 * as.numeric(p > 0 & p <= 0.001) +
    3 * as.numeric(p > 0.001 & p <= 0.01) +
    4 * as.numeric(p > 0.01 & p <= 0.05) +
    5 * as.numeric(p > 0.05 & p <= 0.1) +
    6 * as.numeric(p > 0.1)
cat("names", rownames(s\$conf.int), "\\n", sep="\\t")
cat("exp(coef)", s\$conf.int[1], "\\n", sep="\\t")
cat("lower .95", s\$conf.int[3], "\\n", sep="\\t")
cat("upper .95", s\$conf.int[4], "\\n", sep="\\t")
cat("cp", p, "\\n", sep="\\t")
cat("c", codes[pg], "\\n", sep="\\t")
END
;
  }
  my $hash = {};
  foreach my $str (keys(%params)) {
    for (my $i = 0; $i < 6; $i++) {
      my $h = <$phr>;
      $h =~ s/[\r\n]//g;
      my ($k, $val) = split("\t", $h);
      push @{$hash->{$k}}, $val;
    }
  }
  close($phr);
  close($phw);
  return $hash;
}

sub printUMStats {
  my $res = shift;
  my $len = 0;
  my @names = @{$res->{"names"}};
  foreach (@names) {
    if ($len < length($_)) {
      $len = length($_);
    }
  }
  $len += 1;
  printf "\%$len\s\tHR(lo-hi)\tp\tc\n", "";
  for (my $i = 0; $i < scalar(@names); $i++) {
    printf "\%$len\s\t%.2f(%.2f-%.2f)\t%.2g\t%3s\n", $res->{"names"}->[$i],
        $res->{"exp(coef)"}->[$i], $res->{"lower .95"}->[$i],
        $res->{"upper .95"}->[$i], $res->{"cp"}->[$i], $res->{"c"}->[$i];
  }
}

sub printUMStatsAll {
  my ($res1, $res2, @list) = @_;
  my $len = 0;
  my @names = @{$res1->{"names"}};
  my $hash1 = {};
  for (my $i = 0; $i < scalar(@names); $i++) {
    $hash1->{$names[$i]} = $i;
    if ($len < length($names[$i])) {
      $len = length($names[$i]);
    }
  }
  $len += 1;
  my @names = @{$res2->{"names"}};
  my $hash2 = {};
  for (my $i = 0; $i < scalar(@names); $i++) {
    $hash2->{$names[$i]} = $i;
  }
  if (0) {
    printf "\%$len\s\tHR(lo-hi)\tp\tc\tHR(lo-hi)\tp\tc\n", "";
    my $nd = 2;
    foreach (@list) {
      my $i1 = $hash1->{$_};
      my $i2 = $hash2->{$_};
      printf "\%$len\s\t\%.". $nd ."f(\%." . $nd . "f-\%.". $nd ."f)\t\%.". $nd ."g\t%3s\t", $res1->{"names"}->[$i1],
             $res1->{"exp(coef)"}->[$i1], $res1->{"lower .95"}->[$i1],
             $res1->{"upper .95"}->[$i1], $res1->{"cp"}->[$i1], $res1->{"c"}->[$i1];
      printf "%.". $nd ."f(%.". $nd ."f-%.". $nd ."f)\t%.". $nd ."g\t%3s\n",
             $res2->{"exp(coef)"}->[$i2], $res2->{"lower .95"}->[$i2],
             $res2->{"upper .95"}->[$i2], $res2->{"cp"}->[$i2], $res2->{"c"}->[$i2];
    }
  }
  if (1) {
    printf "\%$len\s\tHR\t95% CI\tp\tc\tHR\t95% CI\tp\tc\n", "";
    my $nd = 2;
    foreach (@list) {
      my $i1 = $hash1->{$_};
      my $i2 = $hash2->{$_};
      printf "\%$len\s\t\%.". $nd ."f\t\%." . $nd . "f - \%.". $nd ."f\t\%.". $nd ."g\t%3s\t", $res1->{"names"}->[$i1],
             $res1->{"exp(coef)"}->[$i1], $res1->{"lower .95"}->[$i1],
             $res1->{"upper .95"}->[$i1], $res1->{"cp"}->[$i1], $res1->{"c"}->[$i1];
      printf "%.". $nd ."f\t%.". $nd ."f - %.". $nd ."f\t%.". $nd ."g\t%3s\n",
             $res2->{"exp(coef)"}->[$i2], $res2->{"lower .95"}->[$i2],
             $res2->{"upper .95"}->[$i2], $res2->{"cp"}->[$i2], $res2->{"c"}->[$i2];
    }
  }
}

sub printUMStatsTexAll {
  my ($ofile, $res1, $res2, @list) = @_;
  my $len = 0;
  my @names = @{$res1->{"names"}};
  my $hash1 = {};
  for (my $i = 0; $i < scalar(@names); $i++) {
    $hash1->{$names[$i]} = $i;
    if ($len < length($names[$i])) {
      $len = length($names[$i]);
    }
  }
  $len += 1;
  my @names = @{$res2->{"names"}};
  my $hash2 = {};
  for (my $i = 0; $i < scalar(@names); $i++) {
    $hash2->{$names[$i]} = $i;
  }
  open(my $ofh, ">$ofile") || die "Can't write $ofile\n";
  print $ofh <<END
\\begin{tabular}{|l|*{8}{c|}}
\\hline
\\multirow{2}{*}{\\bf Vars} & \\multicolumn{4}{|c|}{\\bf Univariate} & \\multicolumn{4}{|c|}{\\bf Multivariate} \\\\
\\cline{2-9}
& \\bf HR & \\bf 95\\% CI & \\bf p & \\bf c & \\bf HR & \\bf 95\\% CI & \\bf p & \\bf c \\\\
\\hline
END
;
  my $nd = 2;
  foreach (@list) {
    my $i1 = $hash1->{$_};
    my $i2 = $hash2->{$_};
    printf $ofh "\%$len\s & \%.". $nd ."f & \%." . $nd . "f - \%.". $nd ."f & \%.". $nd ."g & %3s &", $res1->{"names"}->[$i1],
           $res1->{"exp(coef)"}->[$i1], $res1->{"lower .95"}->[$i1],
           $res1->{"upper .95"}->[$i1], $res1->{"cp"}->[$i1], $res1->{"c"}->[$i1];
    printf $ofh "%.". $nd ."f & %.". $nd ."f - %.". $nd ."f & %.". $nd ."g & %3s \\\\\n",
           $res2->{"exp(coef)"}->[$i2], $res2->{"lower .95"}->[$i2],
           $res2->{"upper .95"}->[$i2], $res2->{"cp"}->[$i2], $res2->{"c"}->[$i2];
    print $ofh "\\hline\n";
  }
  print $ofh <<END
\\end{tabular}
END
;
  close($ofh);
}

sub analyzeAllIds {
  my ($self, $thr, $ct, $maxt) = @_;
  my @headers = ("id", "pvalue", "exp(coef)", "lower .95", "upper .95",
      "cp", "c", "groups");
  my @ids = keys(%{$self->{"idHash"}});
  print join("\t", "id", "name", "pvalue", "numGroups", @headers[2 .. $#headers]), "\n";
  my ($phr, $phw);
  my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  #$phw = \*STDOUT;
  my $ph = [$phr, $phw];
  my $index = 0;
  foreach my $id (@ids) {
    $index++;
    my $res = $self->plotSurvivalThr("null", $id, $thr, $ct, $maxt, $ph);
    my $n = $self->getName($id);
    if (!defined $res) {
      print join("\t", $id, $n), "\n";
    }
    else {
      print join("\t", $id, $n, map {$res->{$headers[$_]}->[0]} 1 .. $#headers), "\n";
    }
  }
  close($phr);
  close($phw);
}

sub getMeans {
  my ($self, $id, @stages) = @_;
  my $expr = $self->getExprData($id);
  my $thr = $self->getThrData($id);

  my $m = [];
  foreach my $s (@stages) {
    my $m1 = &U::mean([map { $expr->[$_] } @$s]);
    push @$m, $m1;
  }
  my $sorted = 1;
  for (my $i = 1; $i < scalar(@$m); $i++) {
    if ($m->[$i - 1] < $m->[$i]) {
      $sorted = 0;
      last;
    }
  }
  my $down = 1;
  for (my $i = 1; $i < scalar(@$m); $i++) {
    if ($m->[0] < $m->[$i]) {
      $down = 0;
      last;
    }
  }
  my $start = $m->[0];
  my $end = $m->[scalar(@$m) - 1];
  my $diff = sprintf("%.2f", $start - $end);
  return [$id, $self->getName($id), $sorted, $down, $diff, $thr->[0], @$m];
}

sub getMarkers {
  my ($self, $cmd, @stages) = @_;
  my $fh = $self->{'fh'};
  seek($fh, 0, 0);
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @header = split("\t", $head);
  my $res = [];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $m = [];
    foreach my $s (@stages) {
      my $m1 = &U::mean([map { $list[$_] } @$s]);
      push @$m, $m1;
    }
    if ($cmd eq "high") {
      my $thr = $self->getThrData($list[0]);
      my $allHigh = 1;
      for (my $i = 0; $i < scalar(@$m); $i++) {
        if ($m->[$i] <= $thr->[0]) {
           $allHigh = 0;
        }
      }
      my $diff = &U::min($m) - $thr->[0];
      if ($allHigh) {
        push @$res, [$diff, $list[0], $list[1]];
      }
      next;
    }
    if ($cmd eq "low") {
      my $thr = $self->getThrData($list[0]);
      my $allLow = 1;
      for (my $i = 0; $i < scalar(@$m); $i++) {
        if ($m->[$i] > $thr->[0]) {
           $allLow = 0;
        }
      }
      my $diff = &U::min($m) - $thr->[0];
      if ($allLow) {
        push @$res, [$diff, $list[0], $list[1]];
      }
      next;
    }
    my $sorted = 1;
    for (my $i = 1; $i < scalar(@$m); $i++) {
      if ($m->[$i - 1] < $m->[$i]) {
        $sorted = 0;
        last;
      }
    }
    my $down = 1;
    for (my $i = 1; $i < scalar(@$m); $i++) {
      if ($m->[0] < $m->[$i]) {
        $down = 0;
        last;
      }
    }
    my $start = $m->[0];
    my $end = $m->[scalar(@$m) - 1];
    my $diff = sprintf("%.2f", $start - $end);
    if ($cmd eq "All") {
      push @$res, [$diff, $list[0], $list[1]];
    }
    if ($cmd eq "sAll" && $sorted == 1) {
      push @$res, [$diff, $list[0], $list[1]];
    }
    if ($cmd eq "shilo" && $sorted == 1) {
      my $thr = $self->getThrData($list[0]);
      if ($start > $thr->[0] && $end < $thr->[0]) {
        push @$res, [$diff, $list[0], $list[1]];
      }
    }
    if ($cmd eq "shihi" && $sorted == 1) {
      my $thr = $self->getThrData($list[0]);
      if ($start > $thr->[0] && $end > $thr->[0]) {
        push @$res, [$diff, $list[0], $list[1]];
      }
    }
    if ($cmd eq "dhilo" && $down == 1) {
      my $thr = $self->getThrData($list[0]);
      if ($start > $thr->[0] && $end < $thr->[0]) {
        push @$res, [$diff, $list[0], $list[1]];
      }
    }
    if ($cmd eq "dhihi" && $down == 1) {
      my $thr = $self->getThrData($list[0]);
      if ($start > $thr->[0] && $end > $thr->[0]) {
        push @$res, [$diff, $list[0], $list[1]];
      }
    }
    if ($cmd eq "hilo" && $sorted == 1) {
      push @$res, [$diff, $list[0], $list[1]];
    }
  }
  return $res;
}

sub getMarkersPerc {
  my ($self, $cmd, $stage, $thrcode) = @_;
  my $fh = $self->{'fh'};
  seek($fh, 0, 0);
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @header = split("\t", $head);
  my $res = [];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $thr1 = $self->getThrData($list[0]);
    my $thr = &getThrCode($thr1, $thr1->[0], $thrcode);
    my $count = 0;
    for (my $i = 0; $i < scalar(@$stage); $i++) {
      if ($cmd eq "high" && $list[$stage->[$i]] > $thr) {
        $count++;
      }
      if ($cmd eq "low" && $list[$stage->[$i]] <= $thr) {
        $count++;
      }
    }
    my $diff = sprintf("%.2f", $count/scalar(@$stage));
    push @$res, [$diff, $list[0], $list[1]];
  }
  return $res;
}

sub readThr {
  my ($self, $thrfile) = @_;
  open(my $fh, "<$thrfile") || die "Can't open $thrfile\n";
  my $hash = {};
  while (<$fh>) {
    s/[\r\n]//g;
    my ($id, $thr1, $stat, $thr0, $thr2) = split("\t");
    $hash->{$id} = [$thr1, $stat, $thr0, $thr2];
  }
  close($fh);
  $self->{'thrHash'} = $hash;
}

sub readThrOld {
  my ($self, $idxfile, $thrfile) = @_;
  open(my $fh, "<$idxfile") || die "Can't open $idxfile\n";
  open(my $fh1, "<$thrfile") || die "Can't open $thrfile\n";
  my $head = <$fh>;
  my $probe = <$fh>;
  my $hash = {};
  while (<$fh>) {
    s/[\r\n]//g;
    my ($id, @list) = split("\t");
    my $line = <$fh1>;
    $line =~ s/[\r\n]//g;
    my ($n, @thrlist) = split("\t", $line);
    $hash->{$id} = [@thrlist];
  }
  close($fh);
  close($fh1);
  $self->{'thrHash'} = $hash;
}

sub subplot {
  my ($self, $cmd, $list) = @_;
  my $idlist = [];
  for (my $i = 0; $i < scalar(@$list); $i++) {
    my $l = $self->getIDs($list->[$i]);
    push @$idlist, @$l;
  }
  print <<END
<html>
<head> <title> subplot </title> </head>
<body>
<table border=0> <tr>
END
;
  if (!defined $cmd || $cmd eq "" || $cmd eq "linear") {
    my $id0 = $idlist->[0];
    my $expr = $self->{'expr'};
    for (my $i = 1; $i < scalar(@$idlist); $i++) {
      my $id1 = $idlist->[$i];
      next if (!defined $id1);
      my $xn = $self->getName($id0);
      my $yn = $self->getName($id1);
      $xn =~ s/ .*//g;
      $yn =~ s/ .*//g;
      my $x = $self->{"idHash"}->{$id0}->[0];
      my $y = $self->{"idHash"}->{$id1}->[0];
      my $url = "http://genepyramid.stanford.edu/microarray/Explore/explore.php?go=plot&file=$expr&xn=$xn&yn=$yn&x=$x&y=$y";
      print "<td> $xn and $yn ($id0, $id1) <br/>";
      print "<img height=240 width=320 src=$url/>";
      print "</td>";
      if ( ($i % 2) == 0 ) {
        print "</tr>\n<tr>";
      }
    }
  }
  if ($cmd eq "pairs") {
    my $expr = $self->{'expr'};
    for (my $i = 0; $i < scalar(@$idlist); $i+=2) {
      my $id0 = $idlist->[$i];
      my $id1 = $idlist->[$i+1];
      next if (!defined $id1);
      next if (!defined $id0);
      my $xn = $self->getName($id0);
      my $yn = $self->getName($id1);
      $xn =~ s/ .*//g;
      $yn =~ s/ .*//g;
      my $x = $self->{"idHash"}->{$id0}->[0];
      my $y = $self->{"idHash"}->{$id1}->[0];
      my $url = "http://genepyramid.stanford.edu/microarray/Explore/explore.php?go=plot&file=$expr&xn=$xn&yn=$yn&x=$x&y=$y";
      print "<td> $xn and $yn ($id0, $id1) <br/>";
      print "<img height=240 width=320 src=$url/>";
      print "</td>";
      if ( ($i % 2) == 0 ) {
        print "</tr>\n<tr>";
      }
    }
  }
  print <<END
</tr> </table>
</body>
</html>
END
;
}

sub boxplot {
  my ($self, $outfile, $id, $groups, $ci, $ph) = @_;
  $ci = 0 if (!defined $ci);
  my $expr = $self->getExprData($id);
  my $thr = $self->getThrData($id);
  my $name = $self->getName($id);
  my $res = [];
  my $ncdata = [];
  my $numgroups = 0;
  for (my $i = 0; $i < scalar(@$groups); $i++) {
    my $g = $groups->[$i];
    my $count = 0;
    foreach my $j (@{$g->[2]}) {
      next if (!defined $expr->[$j] || $expr->[$j] eq "");
      push @$res, [$i+1, $expr->[$j], $j];
      $count++;
    }
    if ($count > 0) {
      push @$ncdata, $g;
      $numgroups++;
    }
  }
  if ($numgroups <= 1) {
    return undef;
  }
  my ($phr, $phw);
  if (defined $ph) {
    ($phr, $phw) = @$ph;
  }
  else {
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  #$phw = \*STDOUT;
  print $phw "groups <- c(", join(",", map {$_->[0]} @$res), ")\n";
  print $phw "expr <- c(", join(",", map {$_->[1]} @$res), ")\n";
  print $phw "indices <- c(", join(",", map {$_->[2]} @$res), ")\n";
  print $phw "nm <- c(", join(",", map {"\"".$_->[0]."\""} @$ncdata), ")\n";
  print $phw "clr <- c(", join(",", map {"\"".$_->[1]."\""} @$ncdata), ")\n";

  my $devstr = "";
  if (defined $outfile) {
    my @l = split("\\.", $outfile);
    if ($l[$#l] eq "png") {
      $devstr = "png(filename=\"$outfile\", width=640, height=480, pointsize=15)";
    }
    if ($l[$#l] eq "ps") {
      $devstr = "postscript(file=\"$outfile\", width=6.4, height=4.8, pointsize=11)";
    }
    if ($l[$#l] eq "pdf") {
      $devstr = "pdf(file=\"$outfile\", width=6.4, height=4.8, pointsize=11)";
    }
  }


print $phw <<END
$devstr
END
;
  if ($devstr ne "") {
  if ($ci == 0) {
  print $phw <<END
par(font.lab=2)
boxplot(expr ~ groups, names=nm, col=clr, notch=FALSE, ylab=\"$id : $name Gene Expression\")
dev.off()
END
;
  }
  else {
  print $phw <<END
par(font.lab=2)
require(gplots)
tmp   <- split(expr, groups)
means <- sapply(tmp, mean)
stdev <- sqrt(sapply(tmp, var))
n     <- sapply(tmp,length)
ciw   <- qt(0.975, n) * stdev / sqrt(n)
num   <- length(nm)
yl <- range($thr->[0], means, means+ciw, means-ciw, na.rm = TRUE)
pr1 <- plotCI(x=means, uiw=ciw, col=clr, xlab=\"\", ylab=\"$id : $name Gene Expression\", barcol=clr, lwd=3, xaxt=\"n\", xlim=c(0, num+1), ylim=yl, yaxs="r")
axis(side=1, at=1:num, labels=nm, cex=0.7)
abline(h=$thr->[0], col='red')
abline(h=$thr->[2], col='cyan')
abline(h=$thr->[3], col='cyan')
d <- means+ciw
d[is.na(d)] <- means[is.na(d)]
text(x=1:num, y=d, labels=paste("n=", n, sep=""), pos=3, offset=0.5)
dev.off()
END
;
  }
  }
print $phw <<END
f <- factor(groups)
l <- levels(f)
for (x in l) {
  if (sum(groups==x) <= 1) { next; }
  ex <- expr[groups==x]
  if (sd(ex) <= 0) {
    cat("mean", x, "\t", ex[1], "\\n")
    cat("conf.int", x, "\t", ex[1], ex[1], "\\n")
  }
  else {
    st <- t.test(ex)
    cat("mean", x, "\t", st\$estimate, "\\n")
    cat("conf.int", x, "\t", st\$conf.int, "\\n")
  }
  for (y in l) {
    if (x >= y) { next; }
    if (sum(groups==y) <= 1) { next; }
    ey <- expr[groups==y]
    if (sd(ey) <= 0 && sd(ex) <= 0) {
      cat("statistic", x, y, "\t", Inf, "\\n")
      if (ex[1] != ey[1]) {
        cat("pvalue", x, y, "\t", 0, "\\n")
      }
      else {
        cat("pvalue", x, y, "\t", 1, "\\n")
      }
    }
    else {
      st <- t.test(ex, ey)
      cat("statistic", x, y, "\t", st\$statistic, "\\n")
      cat("pvalue", x, y, "\t", st\$p.value, "\\n")
    }
  }
}
cat("END\\n")
END
;
  if ($devstr ne "") {
    my $h = <$phr>;
    my $h = <$phr>;
  }
  my $hash = {};
  my $num = scalar(@$groups);
  for (my $i = 0; $i < (2*$num*($num - 1)/2+2*$num); $i++) {
    my $h = <$phr>;
    $h =~ s/[\r\n]//g;
    if ($h eq "END") { last; }
    my ($k, $val) = split("\t", $h, 2);
    push @{$hash->{$k}}, $val;
  }
  if (!defined $ph) {
    close($phr);
    close($phw);
  }
  return $hash;
}

sub getDynamicRangeInfo {
#AffyID name thr mean mean-thr perc min max sd #thrNum #hi #int #lo
  my ($self, $id, $groups) = @_;
  my $expr = $self->getExprData($id);
  my $thr = $self->getThrData($id);
  my $name = $self->getName($id);
  my $res = [];
  my $ng = 1;
  if (defined $groups) {
    $ng = scalar(@$groups);
    for (my $i = 0; $i < scalar(@$groups); $i++) {
      my $g = $groups->[$i];
      foreach my $j (@{$g->[2]}) {
        next if (!defined $expr->[$j] || $expr->[$j] eq "");
        push @$res, [$i+1, $expr->[$j], $j];
      }
    }
  }
  else {
    for (my $i = $self->getStart(); $i < $self->getNum(); $i++) {
      next if (!defined $expr->[$i] || $expr->[$i] eq "");
      push @$res, [1, $expr->[$i], $i];
    }
  }
  my $res1 = [];
  for (my $i = 0; $i < $ng; $i++) {
    my @keys = map { $_->[2] } grep { $_->[0] == ($i+1) } @$res;
    my $sum = 0;
    my $sum2 = 0;
    my $num = 0;
    my $min = 1000; my $max = 0;
    foreach (@keys) {
      $sum += $expr->[$_];
      $sum2 += $expr->[$_]*$expr->[$_];
      if ($min > $expr->[$_]) {
        $min = $expr->[$_];
      }
      if ($max < $expr->[$_]) {
        $max = $expr->[$_];
      }
      $num ++;
    }
    if ($num <= 0) {
      $res1->[$i] = [$id, $name, $thr->[0], 0, - $thr->[0], -1, 0, 0, 0, 0, 0, 0];
      next;
    }
    my $mean = $sum / $num;
    my $mean2 = $sum2 / $num;
    my $sd = 0;
    if ($mean2 > ($mean*$mean)) {
      $sd = sqrt($mean2 - $mean*$mean);
    }
    my $num = 0;
    my $thrNum = 0;
    my $total = scalar(@keys);
    my $high = 0;
    my $int = 0;
    my $low = 0;
    foreach my $i (@keys) {
      if ($expr->[$i] < $mean) { $num ++; }
      if ($expr->[$i] < $thr->[0]) { $thrNum ++; }
      if ($expr->[$i] < $thr->[2]) { $low ++; }
      if ($expr->[$i] >= $thr->[3]) { $high ++; }
      if ($expr->[$i] >= $thr->[2] && $expr->[$i] < $thr->[3]) { $int ++; }
    }
    if ($num > $thrNum) {
      $total = $total - $thrNum;
    }
    else {
      $total = $thrNum;
    }
    my $perc;
    $perc = ($num - $thrNum) / $total if $total > 0;
    $res1->[$i] = [$id, $name, $thr->[0], $mean,
      $mean - $thr->[0], $perc, $min, $max, $sd, $thrNum, $high, $int, $low];
  }
  return $res1;
}

sub getDynamicRangeVInfo {
  my ($self, $id, $groups) = @_;
  my $maxPoints = 20; # -100 to 1, 0, 1, to 100
  my $gmin = 1.0;
  my $gmax = 16.0;
  my $gstep = 0.2;
  my $gnum = int(($gmax - $gmin)/$gstep);
  my $expr = $self->getExprData($id);
  my $thrdata = $self->getThrData($id);
  my $pval = $self->{'thrStatPval'}->{$id};
  my $res = [];
  my $ng = 1;
  if (defined $groups) {
    $ng = scalar(@$groups);
    for (my $i = 0; $i < scalar(@$groups); $i++) {
      my $g = $groups->[$i];
      foreach my $j (@{$g->[2]}) {
        next if (!defined $expr->[$j] || $expr->[$j] eq "");
        push @$res, [$i+1, $expr->[$j], $j];
      }
    }
  }
  else {
    for (my $i = $self->getStart(); $i < $self->getNum(); $i++) {
      next if (!defined $expr->[$i] || $expr->[$i] eq "");
      push @$res, [1, $expr->[$i], $i];
    }
  }
  my $vinfo = [];
  for (my $i = 0; $i < $ng; $i++) {
    my @keys = map { $_->[2] } grep { $_->[0] == ($i+1) } @$res;
    my @expr = map { $expr->[$_] } @keys;
    my @sortedExpr = sort {$a <=> $b} @expr;
    my $thr = $thrdata->[0];
    my ($min, $max, $mean, $sd, $num, $thrNum) = &U::getRangeInfo(\@expr, $thr);
    my @sortedExprBelowThr = grep { $_ < $thr } @sortedExpr;
    my @sortedExprAboveThr = grep { $_ >= $thr } @sortedExpr;
    push @{$vinfo->[$i]}, $id, $min;
    my @values = (0.005, 0.01, 0.025, 0.975, 0.99, 0.995);
    my @indices = map { int($_*$#sortedExpr) } @values;
    my @data = map { $sortedExpr[$_] } @indices;
    push @{$vinfo->[$i]}, @data;
    push @{$vinfo->[$i]}, $max, $mean, $sd, $thr, $num - $thrNum, $pval;
    for (my $j = 0; $j < $maxPoints; $j++) {
      if (scalar(@sortedExprBelowThr) != 0) {
        my $index = int($j*1.0*$#sortedExprBelowThr/$maxPoints);
        my $e = $sortedExprBelowThr[$index];
        push @{$vinfo->[$i]}, sprintf("%.2f", $e);
      }
      else {
        push @{$vinfo->[$i]}, 0;
      }
    }
    for (my $j = 0; $j <= $maxPoints; $j++) {
      if (scalar(@sortedExprAboveThr) != 0) {
        my $index = int($j*1.0*$#sortedExprAboveThr/$maxPoints);
        my $e = $sortedExprAboveThr[$index];
        if ($j == 0) {
          $e = $thr;
        }
        push @{$vinfo->[$i]}, sprintf("%.2f", $e);
      }
      else {
        push @{$vinfo->[$i]}, 0;
      }
    }
    my $min_below = $sortedExprBelowThr[0];
    my $max_below = $sortedExprBelowThr[$#sortedExprBelowThr];
    my $last_index = 0;
    for (my $j = 0; $j < $maxPoints; $j++) {
      if (scalar(@sortedExprBelowThr) != 0) {
        my $e = $min_below + $j*($max_below - $min_below)/$maxPoints;
        my $index = &U::BinarySearch(\@sortedExprBelowThr, $e);
        push @{$vinfo->[$i]}, sprintf("%d", $index - $last_index);
        $last_index = $index;
      }
      else {
        push @{$vinfo->[$i]}, 1;
      }
    }
    my $min_above = $sortedExprAboveThr[0];
    my $max_above = $sortedExprAboveThr[$#sortedExprAboveThr];
    for (my $j = 0; $j <= $maxPoints; $j++) {
      if (scalar(@sortedExprAboveThr) != 0) {
        my $e = $min_above + $j*($max_above - $min_above)/$maxPoints;
        my $index = &U::BinarySearch(\@sortedExprAboveThr, $e);
        my $new_index = $index + scalar(@sortedExprBelowThr);
        push @{$vinfo->[$i]}, sprintf("%d", $new_index - $last_index);
        $last_index = $new_index;
      }
      else {
        push @{$vinfo->[$i]}, 1;
      }
    }
    my $last_index = 0;
    for (my $j = 0; $j <= $gnum; $j++) {
      if (scalar(@sortedExpr) != 0) {
        my $e = $gmin + $j*$gstep;
        my $index = &U::BinarySearch(\@sortedExpr, $e);
        if ($index < 0) {
          $index = 0;
        }
        push @{$vinfo->[$i]}, sprintf("%d", $index - $last_index);
        $last_index = $index;
      }
      else {
        push @{$vinfo->[$i]}, 0;
      }
    }
  }
  return $vinfo;
}

sub plotStepMiner {
  my ($self, $outfile, $id, $groups, $w, $h, $lw, $bg, $thr, $expr, $cmd) = @_;
  $w = 6.4 if (!defined $w);
  $h = 4.8 if (!defined $h);
  $lw = 0.5 if (!defined $lw);
  $bg = "#444444" if (!defined $bg);
  $expr = $self->getExprData($id) if (!defined $expr);
  $thr = $self->getThrData($id) if (!defined $thr);
  $cmd = "python" if (!defined $cmd);
  my $xn = $self->getName($id);
  my ($phr, $phw);
  my $pid = open2($phr, $phw, $cmd);
  #$phw = \*STDOUT;
  print $phw <<END
import matplotlib
matplotlib.use('agg')
import re
from pylab import *
from numpy import *
fig = figure(figsize=($w,$h))
ax = fig.add_axes([0.11, 0.11, 0.78, 0.78])
END
;
  if (!defined $groups) {
    my $data = [];
    for (my $i = $self->getStart(); $i < $self->getNum(); $i++) {
      next if (!defined $expr->[$i] || $expr->[$i] eq "");
      push @$data, $expr->[$i];
    }
    print $phw "data = [", join(",", @$data), "]\n";
    print $phw <<END
c = '$bg'
data.sort()
ax.plot(range(0,len(data)),data, color=c, ls='None', marker='+', mew=1.1, ms=4, mec=c)
END
;
  }
  else {
    my $data = [];
    my $hash = {};
    for (my $i = 0; $i < scalar(@$groups); $i++) {
      my $g = $groups->[$i];
      foreach my $j (@{$g->[2]}) {
        next if (!defined $expr->[$j] || $expr->[$j] eq "");
        push @$data, [$i+1, $expr->[$j], $g->[1], $j];
        $hash->{$j} = 1;
      }
    }
    for (my $i = $self->getStart(); $i < $self->getNum(); $i++) {
      next if (!defined $expr->[$i] || $expr->[$i] eq "");
      next if (defined $hash->{$i});
      push @$data, [0, $expr->[$i], $bg, $i];
    }
    my @sorted = sort { $a->[1] <=> $b->[1] } @$data;
    for (my $i = 0; $i < scalar(@$data); $i++) {
      $sorted[$i]->[3] = $i;
    }
    for (my $i = 0; $i <= scalar(@$groups); $i++) {
      print $phw "data = [", join(",", map { $_->[1] } grep { $_->[0] == $i } @sorted), "]\n";
      print $phw "idx = [", join(",", map { $_->[3] } grep { $_->[0] == $i } @sorted), "]\n";
      if ($i > 0) {
        print $phw "clr = '", $groups->[$i - 1]->[1], "'\n";
      }
      else {
        print $phw "clr = '$bg'\n";
      }
      print $phw "ax.plot(idx,data, color=clr, ls='None', marker='+', mew=1.1, ms=4, mec=clr)\n";
    }
    print $phw "data = [", join(",", map { $_->[1] } @sorted), "]\n";
  }
  print $phw <<END
i = len([x for x in data if x < $thr->[0]])
m1 = mean([x for x in data if x < $thr->[0]])
m2 = mean([x for x in data if x >= $thr->[0]])
ax.plot([0, i], [m1, m1], linewidth=$lw, color='b')
ax.plot([i, len(data)], [m2, m2], linewidth=$lw, color='b')
ax.plot([i, i], [m1, m2], linewidth=$lw, color='b')
ax.plot([0, len(data)], [$thr->[0], $thr->[0]], linewidth=$lw, color='r')
ax.plot([0, len(data)], [$thr->[2], $thr->[2]], linewidth=$lw, color='c')
ax.plot([0, len(data)], [$thr->[3], $thr->[3]], linewidth=$lw, color='c')
ax.axis([0, len(data), min(data)-0.5, max(data)+0.5])
ax.set_xlabel('Sorted arrays', fontsize=10)
ax.set_ylabel('$id:  $xn expression', fontsize=10)
fig.savefig('$outfile', dpi=100)
END
;
  close($phr);
  close($phw);
  sleep(1);
  chmod 0666, $outfile;
}

sub heatmap1 {
  my ($outfile, $expr, $rows, $columns, $ph, %params) = @_;
  my $nrow = scalar(@$rows);
  my $ncol = scalar(@$columns);
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
  delete $params{'width'};
  delete $params{'height'};
  delete $params{'dpi'};
  delete $params{'wd'};
  delete $params{'ht'};

  my ($phr, $phw);
  if (defined $ph) {
    ($phr, $phw) = @$ph;
  }
  else {
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  my $devstr = "";
  if (defined $outfile) {
    my @l = split("\\.", $outfile);
    if ($l[$#l] eq "png") {
      $devstr = "png(filename=\"$outfile\", width=$width, height=$height, pointsize=15)";
    }
    if ($l[$#l] eq "ps") {
      $devstr = "postscript(file=\"$outfile\", width=$wd, height=$ht, pointsize=11)";
    }
    if ($l[$#l] eq "pdf") {
      $devstr = "pdf(file=\"$outfile\", width=$wd, height=$ht, pointsize=11)";
    }
  }
  my $tmpfile = "__data.txt";
  open($ofh, ">$tmpfile") || die "Can't open $tmpfile\n";
  print $ofh join("\t", @$columns), "\n";
  my $hash = {};
  for (my $i = 0; $i < $nrow; $i++) {
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
  
  print $phw "library(gplots)\n";
  my $br = 50;
  if (defined $params{'breaks'}) { $br = $params{'breaks'}; }
  print $phw "$devstr\n";
  print $phw "cat('BEGIN', '\\n')\n";
  print $phw "data <- read.table('__data.txt', header=TRUE, sep='\\t')\n";
  print $phw "data <- as.matrix(data)\n";
  print $phw "heatmap.2(data, breaks=$br, Rowv=FALSE, Colv=FALSE, symm=FALSE, dendrogram='none', labCol=NULL, cexCol=1.0, col=bluered, trace='none'";
  foreach my $k (keys %params) {
    next if $k eq 'breaks';
    print $phw ",$k=", $params{$k};
  }
  print $phw ")\n";
  print $phw "dev.off()\n";
  close($phw);
  my $h;
  while ($h = <$phr>) {
    last if ($h =~ /^BEGIN/);
  }
  if ($h =~ /^BEGIN/) {
    return 1;
  }
  close($phr);
  return 0;
}

sub heatmap {
  my ($outfile, $expr, $rows, $columns, $ph, %params) = @_;
  my $nrow = scalar(@$rows);
  my $ncol = scalar(@$columns);
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
  delete $params{'width'};
  delete $params{'height'};
  delete $params{'dpi'};
  delete $params{'wd'};
  delete $params{'ht'};

  my ($phr, $phw);
  if (defined $ph) {
    ($phr, $phw) = @$ph;
  }
  else {
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  my $devstr = "";
  if (defined $outfile) {
    my @l = split("\\.", $outfile);
    if ($l[$#l] eq "png") {
      $devstr = "png(filename=\"$outfile\", width=$width, height=$height, pointsize=15)";
    }
    if ($l[$#l] eq "ps") {
      $devstr = "postscript(file=\"$outfile\", width=$wd, height=$ht, pointsize=11)";
    }
    if ($l[$#l] eq "pdf") {
      $devstr = "pdf(file=\"$outfile\", width=$wd, height=$ht, pointsize=11)";
    }
  }
  print $phw "library(gplots)\n";
  print $phw "data <- matrix(c(\n";
  for (my $i = 0; $i < $nrow; $i++) {
    print $phw join(",", @{$expr->[$i]});
    if ($i < ($nrow - 1)) {
        print $phw ",\n";
    }
    else {
        print $phw "\n";
    }
  }
  print $phw "), $nrow, byrow=T)\n";
  my $br = 50;
  if (defined $params{'breaks'}) { $br = $params{'breaks'}; }
  print $phw "$devstr\n";
  print $phw "cat('BEGIN', '\\n')\n";
  print $phw "rownames(data) <- c('", join("','", @$rows),"')\n";
  print $phw "colnames(data) <- c('", join("','", @$columns),"')\n";
  print $phw "heatmap.2(data, breaks=$br, Rowv=FALSE, symm=TRUE, dendrogram='none', labCol=NULL, cexCol=1.0, col=bluered, trace='none'";
  foreach my $k (keys %params) {
    next if $k eq 'breaks';
    print $phw ",$k=", $params{$k};
  }
  print $phw ")\n";
  print $phw "dev.off()\n";
  close($phw);
  my $h;
  while ($h = <$phr>) {
    last if ($h =~ /^BEGIN/);
  }
  if ($h =~ /^BEGIN/) {
    return 1;
  }
  close($phr);
  return 0;
}

sub compareIds {
  my ($self, $id1, $id2) = @_;
  my $data1 = $self->getExprData($id1);
  my $data2 = $self->getExprData($id2);
  if (!defined $data1 || !defined $data2) {
    return 0;
  }
  my $thr1 = $self->getThrData($id1);
  my $thr2 = $self->getThrData($id2);
  my $count1 = 0;
  my $count2 = 0;
  for (my $i = $self->getStart(); $i <= $self->getEnd(); $i++) {
    next if (!defined $data1->[$i] || $data1->[$i] eq "");
    next if (!defined $data2->[$i] || $data2->[$i] eq "");
    next if ($data1->[$i] < $thr1->[3]);
    next if ($data2->[$i] < $thr2->[3]);
    if ($data1->[$i] >= $data2->[$i]) {
        $count1 ++;
    }
    else {
        $count2 ++;
    }
  }
  if ($count1 >= $count2) {
    return 1;
  }
  return 0;
}

sub sortIds {
  my ($self, $list) = @_;
  my @ids = @$list;
  if (scalar(@ids) <= 0) {
    return undef;
  }
  if (scalar(@ids) == 1) {
    return [$ids[0]];
  }
  my @sortedIds = sort { $self->compareIds($b, $a) } @ids;
  return [@sortedIds];
}

sub getBestId {
  my ($self, $name) = @_;
  my $list = $self->getIDs($name);
  my $sortedIds = $self->sortIds($list);
  if (!defined $sortedIds || scalar(@$sortedIds) <= 0) {
    return undef;
  }
  return $sortedIds->[0];
}

sub printInfo {
  my ($h, $groups, $rest) = @_;
  my @ids = keys(%{$h->{"idHash"}});
  my $hd ="AffyID name thr mean mean-thr perc min max sd thrNum hi int lo";
  print join("\t", split(" ", $hd)), "\n";
  for (my $j = 0; $j < scalar(@ids) ; $j++) {
    my $id2 = $ids[$j];
    my $res = $h->getDynamicRangeInfo($id2, $groups);
    print join("\t", map { $res->[0]->[$_] } 0 .. 1), "\t";
    print join("\t", map { sprintf("%.2f", $_) } map { $res->[0]->[$_] } 2 .. 8), "\t";
    print join("\t", map { $res->[0]->[$_] } 9 .. 12), "\n";
    if (($j % 1000) == 0) {
        print STDERR "$j\n";
    }
  }
}

sub printVInfo {
  my ($h, $groups, $rest) = @_;
  my $maxPoints = 20; # -100 to 1, 0, 1, to 100
  my $gmin = 1.0;
  my $gmax = 16.0;
  my $gstep = 0.2;
  my $gnum = int(($gmax - $gmin)/$gstep);
  my @h1 = map { sprintf("Perc -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
  my @h2 = map { sprintf("Perc %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
  my @h3 = map { sprintf("Expr -%.2f", ($maxPoints - $_)*1.0/$maxPoints) } 0 .. ($maxPoints - 1);
  my @h4 = map { sprintf("Expr %.2f", $_*1.0/$maxPoints) } 0 .. $maxPoints;
  my @h5 = map { sprintf("Abs %.2f", $_*$gstep + $gmin) } 0 .. $gnum;
  print join("\t", $h->{'headers'}->[0], "min");
  print "\t", join("\t", "0.5\%", "1\%", "2.5\%", "97.5\%", "99\%", "99.5\%");
  print "\t", join("\t", "max", "mean", "sd", "thr", "thrNum", "pval");
  print "\t", join("\t", @h1), "\t", join("\t", @h2);
  print "\t", join("\t", @h3), "\t", join("\t", @h4);
  print "\t", join("\t", @h5), "\n";
  my @ids = keys(%{$h->{"idHash"}});
  for (my $j = 0; $j < scalar(@ids) ; $j++) {
    my $id2 = $ids[$j];
    my $res = $h->getDynamicRangeVInfo($id2, $groups);
    print join("\t", @{$res->[0]}), "\n";
    if (($j % 1000) == 0) {
        print STDERR "$j\n";
    }
  }
}

sub printBooleanImplication {
  my ($h, $ifile, $n, $num, $index, $st, $p, $id3,
        $thrx0, $thrx2, $thry0, $thry2, $hithr, $sdthr, $drthr) = @_;
  $hithr = 100 if (!defined $hithr);
  $sdthr = 1 if (!defined $sdthr);
  $drthr = 6 if (!defined $drthr);
  my $ihash = &U::getHash($ifile);
  my @headers = @{$ihash->{'AffyID'}};
  my $ihh = {};
  for (my $i = 0; $i < scalar(@headers); $i++) {
    $ihh->{$headers[$i]} = $i;
  }
  
  my $l = $h->getIDs($n);
  for (my $i = 0; $i < scalar(@$l); $i++) {
    print "$n\[$i] = ", $l->[$i], "\n";
  }
  my $id = $l->[0];
  if ($num < 0) {
    $id = $h->getBestId($n);
  }
  else {
    $id = $l->[$num];
  }
  print "Using $n\[$num] = ", $id, "\n";
  if (defined $id3 && $id3 ne '-') {
    $id3 = $h->getBestId($id3);
    my $bs = $h->getBooleanStats($id, $id3, undef, $thrx0, $thrx2, $thry0, $thry2);
    my $n3 = $h->getName($id3);
    print join("\t", $id3, $n3, map { sprintf("%.2f", $_) } map { ($bs->[2]->[$_], $bs->[3]->[$_]) } 0 .. 3), "\n";
    print "DynR = [id,name,thr,mean,mean-thr,perc,min,max,sd,thrNum,hi,int,lo]\n";
    my $res = $h->getDynamicRangeInfo($id3);
    print "DynR = [", join(",", $id3, $n3, (map { sprintf("%.2f", $res->[0]->[$_]) } 2 .. 8), map { $res->[0]->[$_] } 9 .. 12), "]\n";
    my $res = $h->getDynamicRangeInfo($id);
    print "DynR = [", join(",", $id, $n, (map { sprintf("%.2f", $res->[0]->[$_]) } 2 .. 8), map { $res->[0]->[$_] } 9 .. 12), "]\n";
    return;
  }
  my @ids = keys(%{$h->{"idHash"}});
  print join("\t", "ID", "Name", map { ("ST[$_]", "p[$_]") } 0 .. 3), "\n";
  my $idx = 0;
  for (my $j = 0; $j < scalar(@ids) ; $j++) {
    my $id2 = $ids[$j];
    next if (defined $hithr && $ihash->{$id2}->[$ihh->{'hi'}] < $hithr);
    next if (defined $sdthr && $ihash->{$id2}->[$ihh->{'sd'}] < $sdthr);
    my $dr = $ihash->{$id2}->[$ihh->{'max'}] - $ihash->{$id2}->[$ihh->{'min'}];
    next if (defined $drthr && $dr < $drthr);
    my $n2 = $h->getName($id2);
    my $bs = $h->getBooleanStats($id, $id2, undef, $thrx0, $thrx2, $thry0, $thry2);
    my $str = join("\t", $id2, $n2, map { sprintf("%.2f", $_) } 
        map { ($bs->[2]->[$_], $bs->[3]->[$_]) } 0 .. 3);
    if (defined $index && defined $st && defined $p) {
      if ($index <= 3 && $bs->[2]->[$index] > $st && $bs->[3]->[$index] <= $p) {
        print $str, "\n";
      }
      elsif ($index == 4 && $bs->[2]->[1] > $st && $bs->[3]->[1] <= $p
          && $bs->[2]->[2] > $st && $bs->[3]->[2] <= $p) {
        print $str, "\n";
      }
      elsif ($index == 5 && $bs->[2]->[0] > $st && $bs->[3]->[0] <= $p
          && $bs->[2]->[3] > $st && $bs->[3]->[3] <= $p) {
        print $str, "\n";
      }
    }
    else {
      print $str, "\n";
    }
    if (($idx % 1000) == 0) {
        print STDERR "$idx\n";
    }
    $idx++;
  }
}

sub printAllSurvival {
  my ($h, $ct, $maxt, @rest) = @_;
  my @ids = keys(%{$h->{"idHash"}});
  my @keys = ("pvalue", "exp(coef)", 
      "lower .95", "upper .95", "cp", "c", "groups");
  print "ArrayID\tName\tThr\t", join("\t", @keys), "\n";
  my ($phr, $phw);
  my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  for (my $j = 0; $j < scalar(@ids) ; $j++) {
    my $id = $ids[$j];
    my $n = $h->getName($id);
    my $thr = $h->getThrData($id);
    print "$id\t$n\t", sprintf("%.2f", $thr->[0]);
    my $stage1 = $h->getArraysAll($id, "thr0", "lo");
    my $stage2 = $h->getArraysAll($id, "thr0", "hi", $id, "thr2", "lo");
    my $stage3 = $h->getArraysAll($id, "thr2", "hi");
    my $groups = [ ["$n high :", "green", $stage3],
       ["$n med  :", "blue",  $stage2],
       ["$n low  :", "red",   $stage1] ];
    my $outfile = "survival.png";
    my $res = $h->plotSurvival($outfile, $groups, $ct, $maxt, [$phr, $phw]);
    foreach my $k (@keys) {
      if ($k ne "groups") {
        $res->{$k}->[0] =~ s/\t.*//g;
        if ($k ne "c" && $k ne "cp" && $k ne "pvalue" && $res->{$k}->[0] ne "") {
          print sprintf("\t%.2f", $res->{$k}->[0]);
        }
        else {
          print "\t", $res->{$k}->[0];
        }
      }
      else {
        print "\t", join("\t", @{$res->{$k}});
      }
    }
    print "\n";
  }
  close($phr);
  close($phw);
}

sub histogramArray {
  my ($self, $ofile, $arrnum, %params) = @_;
  my $break = 100;
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
  my $code = "";
  $code = $params{'code'} if (defined $params{'code'});
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
  print $phw "l <- c(";
  my $fh = $self->{'fh'};
  seek($fh, 0, 0);
  my $head = <$fh>;
  my $index = 0;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $v = $list[$arrnum];
    next if (!defined $v || $v eq "");
    if ($params{'logx'}) {
      if ($v > 0) {
        $v = log($v)/log(2);
      }
      else {
        next;
      }
    }
    if ($index <= 0) {
      print $phw $v;
      $index++;
    }
    else {
      print $phw ", ". $v;
      $index++;
    }
    if (($index % 10000) == 0) {
        print STDERR $index, "\n";
    }
  }
  print $phw ");\n";
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
$code
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
$code
h <- hist(l, breaks=$break, freq=TRUE, main="$main", xlab="value", ylab="count")
n <- length(h\$counts)
plot(h\$breaks[1:n], h\$counts, main="$main", xlab="value", ylab="count"$pa)
cat(h\$breaks, "\\n")
cat(h\$counts, "\\n")
dev.off()
END
;
  }
  close($phw) if (!defined $params{'ph'});
  my $breaks = <$phr>;
  my $counts = <$phr>;
  close($phr) if (!defined $params{'ph'});
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

sub histogramID {
  my ($self, $ofile, $id, %params) = @_;
  my $break = 100;
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
  my $code = "";
  $code = $params{'code'} if (defined $params{'code'});
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
  my $e = $self->getExprData($id);
  my $thr = $self->getThrData($id);
  $thr = $params{'thr'} if (defined $params{'thr'});
  my $start = $self->getStart();
  my $end = $self->getEnd();
  my $order = [ $start .. $end ];
  if (defined $params{'group'}) {
    my $hash = {};
    foreach my $g (@{$params{'group'}}) {
      foreach my $i (@{$g->[2]}) {
        if (!defined $hash->{$i}) {
          push @$order, $i;
        }
        $hash->{$i} = 1;
      }
    }
  }
  print $phw "l <- c(", join(",", map { $e->[$_] } @$order), ");\n";
  print $phw "thr <- c(", join(",", map { $thr->[$_] } 0 .. 3), ");\n";
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
$code
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
$code
h <- hist(l, breaks=$break, freq=TRUE, main="$main", xlab="value", ylab="count")
n <- length(h\$counts)
x <- h\$breaks[1:n]
y <- h\$counts
r1 <- x < thr[3]
r2 <- x >= thr[3] & x < thr[1]
r3 <- x >= thr[1] & x < thr[4]
r4 <- x >= thr[4]
col <- rep("red", n)
col[r1] <- "green"
col[r2] <- "cyan"
col[r3] <- "blue"
col[r4] <- "red"
plot(x, y, col=col, main="$main", xlab="value", ylab="count"$pa)
abline(v=thr[1], col="red")
abline(v=thr[3], col="cyan")
abline(v=thr[4], col="cyan")
cat(h\$breaks, "\\n")
cat(h\$counts, "\\n")
lo <- lowess(x, y)
lo <- data.frame(x=x, y=y)
lines(lo, col="blue", main="$main", xlab="value", ylab="count"$pa)
dev.off()
END
;
  }
  close($phw) if (!defined $params{'ph'});
  my $breaks = <$phr>;
  my $counts = <$phr>;
  close($phr) if (!defined $params{'ph'});
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

sub histogramArrays {
  my ($self, $ofile, $arrays, %params) = @_;
  my $break = 100;
  $break = $params{'break'} if (defined $params{'break'});
  my $width = 640;
  my $height = 480;
  my $dpi = 100;
  my $freq = "TRUE";
  my $idhash;
  $width = $params{'width'} if (defined $params{'width'});
  $height = $params{'height'} if (defined $params{'height'});
  $dpi = $params{'dpi'} if (defined $params{'dpi'});
  $freq = $params{'freq'} if (defined $params{'freq'});
  $idhash = $params{'idhash'} if (defined $params{'idhash'});
  my $wd = $width/$dpi;
  my $ht = $height/$dpi;
  $wd = $params{'wd'} if (defined $params{'wd'});
  $ht = $params{'ht'} if (defined $params{'ht'});
  my $code = "";
  $code = $params{'code'} if (defined $params{'code'});
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
  print $phw "$devstr\n";
  foreach my $arrnum (@$arrays) {
    print $phw "l <- c(";
    my $fh = $self->{'fh'};
    seek($fh, 0, 0);
    my $head = <$fh>;
    my $index = 0;
    while (<$fh>) {
      s/[\r\n]//g;
      my @list = split("\t");
      next if (defined $idhash && !defined $idhash->{$list[0]});
      my $v = $list[$arrnum];
      next if (!defined $v || $v eq "");
      if ($params{'logx'}) {
        if ($v > 0) {
          $v = log($v)/log(2);
        }
        else {
          next;
        }
      }
      if ($index <= 0) {
        print $phw $v;
        $index++;
      }
      else {
        print $phw ", ". $v;
        $index++;
      }
      if (($index % 10000) == 0) {
        print STDERR $index, "\n";
      }
    }
    print $phw ");\n";
    my $main = "$ofile";
    if (defined $params{'title'}) {
      $main = $params{'title'};
    }
    my $pa = "";
    if (defined $params{'pa'}) {
      $pa = $params{'pa'};
    }
    if ($arrnum eq $arrays->[0]) {
      print $phw <<END
$code
h <- hist(l, breaks=$break, freq=$freq, main="$main", xlab="value", ylab="count", plot=F)
n <- length(h\$counts)
#lo <- lowess(h\$breaks[2:n], h\$counts[2:n])
lo <- data.frame(x=h\$breaks[2:n], y=h\$counts[2:n])
plot(lo, type="l", col=$arrnum, main="$main", xlab="value", ylab="count"$pa)
m <- sum(lo\$x * t(lo\$y))/sum(lo\$y)
abline(v=m, col=$arrnum)
END
;
    }
    else {
      print $phw <<END
$code
h <- hist(l, breaks=$break, freq=$freq, main="$main", xlab="value", ylab="count", plot=F)
n <- length(h\$counts)
#lo <- lowess(h\$breaks[2:n], h\$counts[2:n])
lo <- data.frame(x=h\$breaks[2:n], y=h\$counts[2:n])
lines(lo, col=$arrnum, main="$main", xlab="value", ylab="count"$pa)
m <- sum(lo\$x * t(lo\$y))/sum(lo\$y)
abline(v=m, col=$arrnum)
END
;
    }
  }
  if (0) {
  my @hdrs = map { $self->{'headers'}->[$_] } @$arrays;
  @hdrs = @{$params{'labels'}} if (defined $params{'labels'});
  print $phw "nm <- c(\"", join("\",\"", @hdrs), "\")\n"; 
  print $phw "clr <- c(", join(",", @$arrays), ")\n"; 
  print $phw <<END
legend("topright", nm, col=clr, lwd=3, inset=0.02)
END
;
  }
  close($phw) if (!defined $params{'ph'});
  close($phr) if (!defined $params{'ph'});
}

sub getGeneSignature {
  my ($h, $ids) = @_;
  my $arrExpr = [];
  foreach my $id (@$ids) {
    my $expr = $h->getExprData($id);
    push @$arrExpr, $expr;
  }
  my $res = ["ID", "GeneSignature"];
  for (my $i = $h->getStart(); $i <= $h->getEnd(); $i++) {
    my @data = map { $_->[$i] } @$arrExpr;
    my @data = grep { defined $_ && $_ ne "" } @data;
    $res->[$i] = &U::mean(\@data);
  }
  return $res;
}

sub printBooleanStats {
  my ($h, $ifile, $lfile, $gfile, $st, $p, $drthr, $pGroups) = @_;
  $st = 3 if (!defined $st);
  $p = 0.1 if (!defined $p);
  my $ghash = &U::getHash($gfile);
  my @ids = keys(%{$ghash});
  my $ihash = &U::getHash($ifile);
  my @headers = @{$ihash->{'AffyID'}};
  my $ihh = {};
  for (my $i = 0; $i < scalar(@headers); $i++) {
    $ihh->{$headers[$i]} = $i;
  }
  my @ids = grep { defined $ihash->{$_} } @ids;
  print STDERR "Num IDs = ", scalar(@ids), "\n";
  my $lhash = &U::getHash($lfile);
  my @reqids = grep { defined $lhash->{$_} } @ids;
  print STDERR "Num outer IDs = ", scalar(@reqids), "\n";
  print join("\t", "ID", "Name", "nhi", "sd", "dr", map { ("$_:index", "$_:st", "$_:p") } @ids), "\n";
  my $idx = 0;
  foreach my $id (@reqids) {
    $idx++;
    print STDERR "$id $idx\n";
    my $nhi = $ihash->{$id}->[$ihh->{'hi'}];
    my $sd = $ihash->{$id}->[$ihh->{'sd'}];
    my $dr = $ihash->{$id}->[$ihh->{'max'}] - $ihash->{$id}->[$ihh->{'min'}];
    next if (defined $drthr && $dr < $drthr);
    print join("\t", $id, $h->getName($id), $nhi, $sd, $dr);
    for (my $j = 0; $j < scalar(@ids) ; $j++) {
      my $id2 = $ids[$j];
      my $bs = $h->getBooleanStats($id, $id2, $pGroups);
      my $bsindex = &getBSindex($bs, $st, $p);
      print "\t", join("\t", @$bsindex);
    }
    print "\n";
  }
}

sub getBSindex {
  my ($bs, $st, $p) = @_;
  my $count = 0;
  my $bindex = 0;
  my $bst = 0;
  my $bp = 0;
  for (my $index = 0; $index < 4; $index++) {
    if ($bs->[2]->[$index] > $st && $bs->[3]->[$index] <= $p) {
      $count++;
      ($bindex, $bst, $bp) = ($index, $bs->[2]->[$index], $bs->[3]->[$index]);
    }
  }
  if ($count == 1) {
    return [$bindex, $bst, $bp];
  }
  elsif ($count == 2 && $bs->[2]->[1] > $st && $bs->[3]->[1] <= $p
      && $bs->[2]->[2] > $st && $bs->[3]->[2] <= $p) {
    return [4, &U::min([$bs->[2]->[1], $bs->[2]->[2]]),
        &U::max([$bs->[3]->[1], $bs->[3]->[2]])];
  }
  elsif ($count == 2 && $bs->[2]->[0] > $st && $bs->[3]->[0] <= $p
      && $bs->[2]->[3] > $st && $bs->[3]->[3] <= $p) {
    return [5, &U::min([$bs->[2]->[0], $bs->[2]->[3]]),
        &U::max([$bs->[3]->[0], $bs->[3]->[3]])];
  }
  return [-1, 0, 0];
}

sub getBvCode {
  my ($val, $thr) = @_;
  if (!defined $val || $val eq "") {
    return " ";
  }
  elsif ($val < $thr->[2]) {
    return "0";
  }
  elsif ($val >= $thr->[3]) {
    return "2";
  }
  else {
    return "1";
  }
}

sub printBvFile {
  my ($self, $rest) = @_;
  my $fh = $self->{'fh'};
  seek($fh, 0, 0);
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @header = split("\t", $head);
  my $start = $self->getStart();
  my $end = $self->getEnd();
  print join("\t", "ArrayID", "Name", "BitVector"), "\n";
  my $res = [];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $id = $list[0];
    my $name = $list[1];
    my $thr = $self->getThrData($id);
    my $bv = join("", map { &getBvCode($list[$_], $thr) } $start .. $end);
    print join("\t", $id, $name, $bv), "\n";
  }
}

sub cluster {
  my ($outfile, $expr, $rows, $columns, $ph, %params) = @_;
  my $nrow = scalar(@$rows);
  my $ncol = scalar(@$columns);
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
  delete $params{'width'};
  delete $params{'height'};
  delete $params{'dpi'};
  delete $params{'wd'};
  delete $params{'ht'};

  my ($phr, $phw);
  if (defined $ph) {
    ($phr, $phw) = @$ph;
  }
  else {
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  my $devstr = "";
  if (defined $outfile) {
    my @l = split("\\.", $outfile);
    if ($l[$#l] eq "png") {
      $devstr = "png(filename=\"$outfile\", width=$width, height=$height, pointsize=15)";
    }
    if ($l[$#l] eq "ps") {
      $devstr = "postscript(file=\"$outfile\", width=$wd, height=$ht, pointsize=11)";
    }
    if ($l[$#l] eq "pdf") {
      $devstr = "pdf(file=\"$outfile\", width=$wd, height=$ht, pointsize=11)";
    }
  }
  my $tmpfile = "__data.txt";
  open($ofh, ">$tmpfile") || die "Can't open $tmpfile\n";
  print $ofh join("\t", @$columns), "\n";
  my $hash = {};
  for (my $i = 0; $i < $nrow; $i++) {
    my $id = $rows->[$i];
    if (!defined $id || $id eq "") { $id = "---"; }
    if (defined $hash->{$id}) {
      $hash->{$id} ++;
      $id = $id.".".$hash->{$id};
    }
    else {
      $hash->{$id} = 0;
    }
    print $ofh join("\t", $id, @{$expr->[$i]}), "\n";
  }
  close($ofh);
  
  if (1) {
    print $phw "library(gplots)\n";
    my $br = 50;
    if (defined $params{'breaks'}) { $br = $params{'breaks'}; }
    print $phw "$devstr\n";
    print $phw "cat('BEGIN', '\\n')\n";
    print $phw "data <- read.table('__data.txt', header=TRUE, sep='\\t')\n";
    print $phw "data <- as.matrix(data)\n";
    if (defined $params{'precode'}) { print $phw $params{'precode'}; }
    print $phw "heatmap.2(data, breaks=$br, col=bluered, trace='none'";
    foreach my $k (keys %params) {
      next if $k eq 'breaks';
      next if $k eq 'code';
      next if $k eq 'out';
      next if $k eq 'precode';
      print $phw ",$k=", $params{$k};
    }
    print $phw ")\n";
    if (defined $params{'code'}) { print $phw $params{'code'}; }
    print $phw "dev.off()\n";
    print $phw "cat('END', '\\n')\n";
  }
  close($phw);
  my $h;
  while ($h = <$phr>) {
    last if ($h =~ /^END/);
  }
  if ($h =~ /^END/) {
    return 1;
  }
  close($phr);
  return 0;
}

sub getFloatStr {
  my $v = shift;
  return sprintf("%.3g", $v);
}

sub getBHdiff {
  my ($h, $g1, $g2, %param) = @_;

  my $alpha = 0.2;
  my @ids = keys(%{$h->{"idHash"}});

  $alpha = $param{'alpha'} if (defined $param{'alpha'});
  @ids = @{$param{'ids'}} if (defined $param{'ids'});

  my $res = [];
  foreach my $id (@ids) {
    my $expr = $h->getExprData($id);
    my $name = $h->getName($id);
    my @l = split(" /// ", $name);
    my $ea = [map {$expr->[$_]} @$g1];
    my $eb = [map {$expr->[$_]} @$g2];
    my $ma = &U::mean($ea);
    my $mb = &U::mean($eb);
    my $p = &U::ttest($ea, $eb);
    push @$res, [$p, $ma, $mb, $expr->[0], $l[0], $expr->[1]];
  }
  my @sres = sort { $a->[0] <=> $b->[0] } @$res;
  my $num = scalar(@sres);
  my @rank = 1 .. $num;
  my @BH = map { ($alpha * $_) / $num } @rank;
  my @FDR = grep { ($sres[$_]->[0] < $BH[$_]) } 0 .. $#sres;
  my @up = grep { $sres[$_]->[1] > $sres[$_]->[2] } @FDR;
  my @down = grep { $sres[$_]->[1] <= $sres[$_]->[2] } @FDR;
  my $tnum = scalar(@up);
  if ($tnum < scalar(@down)) {
    $tnum = scalar(@down);
  }
  my $sig = 1;
  if ($tnum < 1) {
    $tnum = 10;
    $sig = 0;
    @up = grep { $sres[$_]->[1] > $sres[$_]->[2] } 0 .. $#sres;
    @down = grep { $sres[$_]->[1] <= $sres[$_]->[2] } 0 .. $#sres;
  }

  my $results = [ $sig,
        scalar(@up), [ map { [$sres[$up[$_]]->[3], $sres[$up[$_]]->[4],
        &getFloatStr($sres[$up[$_]]->[0]), &getFloatStr($BH[$up[$_]])] } 0 .. $#up ],
        scalar(@down), [ map { [$sres[$down[$_]]->[3], $sres[$down[$_]]->[4],
        &getFloatStr($sres[$down[$_]]->[0]), &getFloatStr($BH[$down[$_]])] } 0 .. $#down ]
    ];
  return $results;
}

sub printBHdiff1 {
  my $res = shift;
  my $tnum = $res->[1];
  if ($tnum < $res->[3]) {
    $tnum = $res->[3];
  }
  my $sig = 1;
  if ($res->[0] eq 0) {
    print "No significance\n";
    $sig = 0;
    $tnum = 20;
  }
  if ($sig == 1) {
    if ($res->[1] > 0) {
      foreach my $i (0 .. ($res->[1] - 1)) {
        print join("\t", @{$res->[2]->[$i]}), "\tUp\n";
      }
    }
    if ($res->[3] > 0) {
      foreach my $i (0 .. ($res->[3] - 1)) {
        print join("\t", @{$res->[4]->[$i]}), "\tDown\n";
      }
    }
  }
  else {
    foreach my $i (0 .. ($tnum - 1)) {
      print join("\t", @{$res->[2]->[$i]}), "\tUp\n";
    }
    foreach my $i (0 .. ($tnum - 1)) {
      print join("\t", @{$res->[4]->[$i]}), "\tDown\n";
    }
  }
}

sub printBHdiff {
  my $res = shift;
  my $tnum = $res->[1];
  if ($tnum < $res->[3]) {
    $tnum = $res->[3];
  }
  my $sig = 1;
  if ($res->[0] eq 0) {
    print "No significance\n";
    $tnum = 20;
  }
  foreach my $i (0 .. ($tnum - 1)) {
    if ($res->[1] > 0) {
      print join("\t", @{$res->[2]->[$i]});
    }
    else {
      print join("\t", map { "" } 1 .. 4);
    }
    print "\t";
    if ($res->[3] > 0) {
      print join("\t", @{$res->[4]->[$i]});
    }
    else {
      print join("\t", map { "" } 1 .. 4);
    }
    print "\n";
  }
}

sub samrAnalysis {
  my ($h, $g1, $g2, $ph, %params) = @_;
  my $delta = 0.1;
  $delta = $params{'delta'} if (defined $params{'delta'});
  my $fdr = 100;
  $fdr = $params{'fdr'} if (defined $params{'fdr'});
  if (defined $params{'fdr'}) { $delta = 0; }

  my $n1 = scalar(@$g1);
  my $n2 = scalar(@$g2);
  use File::Temp;
  my $ft = File::Temp->new(
      UNLINK   => 0,
      TEMPLATE => 'tmp-XXXX',
      SUFFIX => '-expr.txt'
      );
  my $outfile = $ft->filename;
  open(my $fh, ">$outfile") || die "Can't open $outfile\n";
  my @hdr = map { $h->{"headers"}->[$_] } (@$g1, @$g2);
  print $fh "ID\tName\t", join("\t", @hdr), "\n";
  my @ids = keys(%{$h->{"idHash"}});
  my @columnNum = map { 0 } 3 .. ($n1 + $n2 + 2);
  my $totalNum = 0;
  foreach my $id (@ids) {
    my $e = $h->getExprData($id);
    my @gr1 = map { $e->[$_] } @$g1;
    my @gr2 = map { $e->[$_] } @$g2;
    if (defined $params{"removeSparse"}) {
        my $ne1 = scalar(grep { $_ eq "" } @gr1);
        my $ne2 = scalar(grep { $_ eq "" } @gr2);
        my $nm = &U::max([$ne1/$n1, $ne2/$n2]);
        if ($nm > $params{"removeSparse"}) {
          next;
        }
    }
    print $fh join("\t", $e->[0], $e->[1], @gr1), "\t", join("\t", @gr2), "\n";
    $totalNum++;
    if (defined $params{"removeColumn"}) {
      foreach (grep { $gr1[$_] eq "" } 0 .. $#gr1) { $columnNum[$_] ++; }
      foreach (grep { $gr2[$_] eq "" } 0 .. $#gr2) { $columnNum[$n1 + $_] ++; }
    }
  }
  close($fh);
  if (defined $params{"removeColumn"}) {
    foreach my $i (0 .. ($n1 + $n2 - 1)) {
      if (($columnNum[$i]/$totalNum) > $params{"removeColumn"}) {
        print STDERR join("\t", $i, $h->{'hhash'}->{$hdr[$i]}, 
            $hdr[$i], ($columnNum[$i]/$totalNum)), "\n";
      }
    }
  }
  my ($phr, $phw);
  if (defined $ph) {
    ($phr, $phw) = @$ph;
  }
  else {
    my $pid = open2($phr, $phw, "R --slave 2>/dev/null");
  }
  if (!defined $phw) {
    unlink($outfile);
    return;
  }
print $phw <<END
library(samr)
file <- "$outfile"
END
;
  print $phw "ord <- c(", join(",", 3 .. ($n1 + $n2 + 2)), ")\n";
  print $phw "cls <- c(", join(",", (map { 1 } 1 .. $n1), (map {2} 1 .. $n2)), ")\n";
print $phw <<END
header.first <- scan(file, what="", comment.char='#', nlines=1, sep="\\t")
header.second <- scan(file, what="", comment.char='#', skip=1, nlines=1, sep="\\t")
table.data <- read.table(file, comment.char='#', quote="\\"", skip=1, sep="\\t")
xdim <- dim(table.data)
gid <- as.character(table.data[,1])
names <- as.character(table.data[,2])
gdata <- as.matrix(table.data[,ord])
set.seed(100)

data=list(x=gdata,y=cls, geneid=gid,genenames=names, logged2=TRUE)

samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)

png(filename="$outfile\.results.png", width=640, height=480, pointsize=15)
opar <- par(no.readonly = TRUE)
# create an array of plots
par(mfrow=c(1,2))

delta=$delta
samr.plot(samr.obj,delta)

par(opar)
dev.off()

delta.table <- samr.compute.delta.table(samr.obj)
dt <- delta.table[,c(1,5)]
dt

cat("Significant :\\n")
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, data, delta.table, all.genes=F)
cat("Significant Genes :\\n")
up = NULL
if (!is.null(siggenes.table\$genes.up)) {
  up <- siggenes.table\$genes.up[,8] < $fdr
}
down = NULL
if (!is.null(siggenes.table\$genes.lo)) {
  down <- siggenes.table\$genes.lo[,8] < $fdr
}
cat("--BEGIN\\n")
write.table(t(c("Up Regulated Genes", sum(up))), quote=F, file="",
    sep="\\t", na = " ", row.names=F, col.names=F, append=F)
if (!is.null(up) & sum(up) != 0) {
write.table(siggenes.table\$genes.up[up,], quote=F, file="",
    sep="\\t", na = " ", row.names=F, col.names=T, append=T)
}
write.table(t(c("Down Regulated Genes", sum(down))), quote=F, file="",
    sep="\\t", na = " ", row.names=F, col.names=F, append=T)
if (!is.null(down) & sum(down) != 0) {
write.table(siggenes.table\$genes.lo[down,], quote=F, file="",
    sep="\\t", na = " ", row.names=F, col.names=T, append=T)
}
cat("--END\\n")
END
;
  my $start = 0;
  while (my $in = <$phr>) {
    last if ($in =~ /^--END/);
    if ($start == 1) {
      print $in;
    }
    $start = 1 if ($in =~ /^--BEGIN/);
  }

  if (!defined $params{"keep_tmps"}) {
    unlink($outfile);
    unlink("$outfile\.results.png");
  }
}

sub printRankFile {
  my ($h, $g1, $g2, $ofile, $rest) = @_;
  open(my $fh, ">$ofile") || die "Can't open $ofile\n";
  my @ids = keys(%{$h->{"idHash"}});
  foreach my $id (@ids) {
    my $e = $h->getExprData($id);
    my $m1 = &U::mean([ map { $e->[$_] } @$g1]);
    my $m2 = &U::mean([ map { $e->[$_] } @$g2]);
    print $fh join("\t", $e->[0], $m1 - $m2), "\n";
  }
  close($fh);
}

sub selectLinear {
   my ($h, $id1, $id2, $d, $p1, $p2, $rest) = @_;
   $d = 0.5 if (!defined $d);
   my $e1 = $h->getExprData($id1);
   my $e2 = $h->getExprData($id2);
   my $start = $h->getStart();
   my $end = $h->getEnd();
   my $ex1 = [map {$e1->[$_]} $start .. $end];
   my $ex2 = [map {$e2->[$_]} $start .. $end];
   my $min1 = &U::min($ex1);
   my $max1 = &U::max($ex1);
   my $min2 = &U::min($ex2);
   my $max2 = &U::max($ex2);
   $p1 = [$min1, $min2] if (!defined $p1);
   $p2 = [$max1, $max2] if (!defined $p2);
   if (($p2->[0] - $p1->[0]) == 0) {
     print STDERR "select Linear div by zero\n";
     exit(1);
   }
   $slope = ($p2->[1] - $p1->[1])/($p2->[0] - $p1->[0]);
   my $res = [];
   foreach my $i ($start .. $end) {
     next if ($e1->[$i] > $p2->[0]);
     next if ($e2->[$i] > $p2->[1]);
     next if ($e1->[$i] < $p1->[0]);
     next if ($e2->[$i] < $p1->[1]);
     my $y = $p1->[1] + $slope * ($e1->[$i] - $p1->[0]);
     if (abs($y - $e2->[$i]) <= $d) {
        push @$res, $i;
     }
   }
   return $res;
}

sub plotSingle {
  my ($outfile, $expr, $groups, %params) = @_;
  my $w = 9;
  my $h = 3;
  my $bg = "#444444";
  my $cmd = "python";
  my $xn = "Gene";
  my $axes = "[0.1, 0.15, 0.85, 0.78]";
  $w = $params{'w'} if (defined $params{'w'});
  $h = $params{'h'} if (defined $params{'h'});
  $bg = $params{'bg'} if (defined $params{'bg'});
  $cmd = $params{'cmd'} if (defined $params{'cmd'});
  $xn = $params{'xn'} if (defined $params{'xn'});
  $axes = $params{'axes'} if (defined $params{'axes'});
  my ($phr, $phw);
  if (defined $params{'ph'}) {
    ($phr, $phw) = @{$params{'ph'}};
  }
  else {
    my $pid = open2($phr, $phw, $cmd);
  }
  print $phw <<END
import matplotlib
matplotlib.use('agg')
import re
from pylab import *
from numpy import *
fig = figure(figsize=($w,$h))
ax = fig.add_axes($axes)
END
;
  if (!defined $groups) {
    my $data = [];
    for (my $i = 2; $i < scalar(@$expr); $i++) {
      next if (!defined $expr->[$i] || $expr->[$i] eq "");
      push @$data, $expr->[$i];
    }
    print $phw "data = [", join(",", @$expr), "]\n";
    print $phw <<END
c = '$bg'
ax.plot(range(1,len(data)+1),data, color=c, ls='None', marker='o', mew=1.1, ms=4, mec=c)
END
;
  }
  else {
    my $data = [];
    my $clr = [];
    for (my $i = 0; $i < scalar(@$groups); $i++) {
      my $g = $groups->[$i];
      my $d = [];
      my $n = scalar(@$data) + 1;
      foreach my $j (@{$g->[2]}) {
        next if (!defined $expr->[$j] || $expr->[$j] eq "");
        push @$data, $expr->[$j];
        push @$d, $expr->[$j];
        push @$clr, "'$g->[1]'";
      }
      print $phw "d = [", join(",", @$d), "]\n";
      print $phw "c = '$g->[1]'\n";
      print $phw <<END
ax.plot(range($n,len(d)+$n),d, c=c, ls='None', marker='o', mew=1.1, ms=4, mec=c)
END
;
    }
    print $phw "data = [", join(",", @$data), "]\n";
  }
  print $phw <<END
ax.axis([0, len(data)+1, min(data)-0.5, max(data)+0.5])
ax.set_xlabel('Array Index', fontsize=10)
ax.set_ylabel('$xn expression', fontsize=10)
fig.savefig('$outfile', dpi=100)
END
;
  if (!defined $params{'ph'}) {
    close($phr);
    close($phw);
  }
  sleep(1);
  chmod 0666, $outfile;
}

sub plotErrorbar {
  my ($outfile, $expr, $lo, $hi, $groups, %params) = @_;
  my $w = 9;
  my $h = 3;
  my $bg = "#444444";
  my $cmd = "python";
  my $xn = "Gene";
  my $axes = "[0.1, 0.15, 0.85, 0.78]";
  $w = $params{'w'} if (defined $params{'w'});
  $h = $params{'h'} if (defined $params{'h'});
  $bg = $params{'bg'} if (defined $params{'bg'});
  $cmd = $params{'cmd'} if (defined $params{'cmd'});
  $xn = $params{'xn'} if (defined $params{'xn'});
  $axes = $params{'axes'} if (defined $params{'axes'});
  my ($phr, $phw);
  if (defined $params{'ph'}) {
    ($phr, $phw) = @{$params{'ph'}};
  }
  else {
    my $pid = open2($phr, $phw, $cmd);
  }
  print $phw <<END
import matplotlib
matplotlib.use('agg')
import re
from pylab import *
from numpy import *
fig = figure(figsize=($w,$h))
ax = fig.add_axes($axes)
END
;
  if (!defined $groups) {
    my $d= [];
    my $d_lo = [];
    my $d_hi = [];
    my $n = 1;
    for (my $i = 2; $i < scalar(@$expr); $i++) {
      next if (!defined $expr->[$i] || $expr->[$i] eq "");
      push @$d, $expr->[$i];
      push @$d_lo, $lo->[$i];
      push @$d_hi, $hi->[$i];
    }
    print $phw "d = asarray([", join(",", @$d), "])\n";
    print $phw "d_lo = asarray([", join(",", @$d_lo), "])\n";
    print $phw "d_hi = asarray([", join(",", @$d_hi), "])\n";
    print $phw <<END
c = '$bg'
ax.errorbar(range($n,len(d)+$n),d, yerr=[d-d_lo, d_hi - d], c=c, ls='None', marker='o', mew=1.1, ms=4, mec=c)
END
;
    print $phw "data = d\n";
    print $phw "data_lo = d_lo\n";
    print $phw "data_hi = d_hi\n";
  }
  else {
    my $data = [];
    my $data_lo = [];
    my $data_hi = [];
    my $clr = [];
    for (my $i = 0; $i < scalar(@$groups); $i++) {
      my $g = $groups->[$i];
      my $d = [];
      my $d_lo = [];
      my $d_hi = [];
      my $n = scalar(@$data) + 1;
      foreach my $j (@{$g->[2]}) {
        next if (!defined $expr->[$j] || $expr->[$j] eq "");
        push @$data, $expr->[$j];
        push @$data_lo, $lo->[$j];
        push @$data_hi, $hi->[$j];
        push @$d, $expr->[$j];
        push @$d_lo, $lo->[$j];
        push @$d_hi, $hi->[$j];
        push @$clr, "'$g->[1]'";
      }
      print $phw "d = asarray([", join(",", @$d), "])\n";
      print $phw "d_lo = asarray([", join(",", @$d_lo), "])\n";
      print $phw "d_hi = asarray([", join(",", @$d_hi), "])\n";
      print $phw "c = '$g->[1]'\n";
      print $phw <<END
ax.errorbar(range($n,len(d)+$n),d, yerr=[d-d_lo, d_hi-d], c=c, ls='None', marker='o', mew=1.1, ms=4, mec=c)
END
;
    }
    print $phw "data = [", join(",", @$data), "]\n";
    print $phw "data_lo = [", join(",", @$data_lo), "]\n";
    print $phw "data_hi = [", join(",", @$data_hi), "]\n";
  }
  print $phw <<END
ax.axis([0, len(data)+1, min(data_lo)-0.5, max(data_hi)+0.5])
ax.set_xlabel('Array Index', fontsize=10)
ax.set_ylabel('$xn expression', fontsize=10)
fig.savefig('$outfile', dpi=100)
END
;
  if (!defined $params{'ph'}) {
    close($phr);
    close($phw);
  }
  sleep(1);
  chmod 0666, $outfile;
}

sub mergeHegemon {
  my ($self, $h, $if, $cf, $ef, $idxf, @rest) = @_;
  my $fh = $self->{'fh'};
  seek($fh, 0, 0);
  open(my $oexprfh, ">$ef") || die "Can't open $ef\n";
  open(my $oidxfh, ">$idxf") || die "Can't open $idxf\n";
  print $oidxfh "ProbeID\tPtr\tName\tDescription\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @header = split("\t", $head);
  my $order = [];
  my @hdr = @{$h->{'headers'}};
  foreach my $i (2 .. $#hdr) {
    if (!defined $self->{'hhash'}->{$hdr[$i]}) {
        push @$order, $i;
    }
  }
  my @headers = (@header, map {$hdr[$_]} @$order);
  $head = join("\t", @headers)."\n";
  print $oexprfh $head;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $ptr = tell($oexprfh);
    my @info = @{$self->{"idHash"}->{$list[0]}};
    my $e = $h->getExprData($list[0]);
    print $oidxfh join("\t", $list[0], $ptr, @info[1,2]), "\n";
    print $oexprfh join("\t", @list, map {$e->[$_]} @$order), "\n";
  }
  close($oexprfh);
  close($oidxfh);
  my @chdrs;
  foreach my $k (keys %{$self->{'survivalHdrs'}}) {
    $chdrs[$self->{'survivalHdrs'}->{$k}] = $k;
  }
  my @chdrs_2;
  my $order = [];
  foreach my $k (keys %{$h->{'survivalHdrs'}}) {
    if (!defined $self->{'survivalHdrs'}->{$k}) {
      push @chdrs_2, $k;
      push @$order, $h->{'survivalHdrs'}->{$k};
    }
  }
  open(my $ocfh, ">$cf") || die "Can't open $cf\n";
  open(my $oifh, ">$if") || die "Can't open $if\n";
  print $oifh "ArrayID\tArrayHeader\tClinicalhHeader\n";
  print $ocfh join("\t", @chdrs, @chdrs_2), "\n";
  foreach my $i (2 .. $#headers) {
    my $arr = $headers[$i];
    my @list = map { "" } (@chdrs, @chdrs_2);
    $list[0] = $arr;
    if (defined $h->{'survivalHash'}->{$arr}) {
      foreach my $k (@chdrs) {
        next if ($k eq "ArrayID");
        if (defined $h->{'survivalHdrs'}->{$k}) {
          my $j1 = $self->{'survivalHdrs'}->{$k};
          my $j2 = $h->{'survivalHdrs'}->{$k};
          $list[$j1] = $h->{'survivalHash'}->{$arr}->[$j2 - 1];
        }
      }
      foreach my $k (0 .. $#chdrs_2) {
        my $j2 = $h->{'survivalHdrs'}->{$chdrs_2[$k]};
        $list[$k + scalar(@chdrs)] = $h->{'survivalHash'}->{$arr}->[$j2 - 1];
      }
    }
    if (defined $self->{'survivalHash'}->{$arr}) {
      foreach my $k (@chdrs) {
        next if ($k eq "ArrayID");
        my $j1 = $self->{'survivalHdrs'}->{$k};
        $list[$j1] = $self->{'survivalHash'}->{$arr}->[$j1 - 1];
      }
    }
    print $ocfh join("\t", @list), "\n";
    my $j = $self->{'survivalHdrs'}->{"c Title"};
    print $oifh join("\t", $arr, $arr, $list[$j]), "\n";
  }
  close($oifh);
  close($ocfh);
}

sub subset {
  my ($self, $pre, $order, @rest) = @_;
  my $if = "$pre\-ih.txt";
  my $cf = "$pre\-survival.txt";
  my $ef = "$pre\-expr.txt";
  my $idxf = "$pre\-idx.txt";
  my $fh = $self->{'fh'};
  seek($fh, 0, 0);
  open(my $oexprfh, ">$ef") || die "Can't open $ef\n";
  open(my $oidxfh, ">$idxf") || die "Can't open $idxf\n";
  print $oidxfh "ProbeID\tPtr\tName\tDescription\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @headers = split("\t", $head);
  $head = join("\t", map { $headers[$_] } (0, 1, @$order))."\n";
  print $oexprfh $head;
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    my $ptr = tell($oexprfh);
    my @info = @{$self->{"idHash"}->{$list[0]}};
    print $oidxfh join("\t", $list[0], $ptr, @info[1,2]), "\n";
    print $oexprfh join("\t", map {$list[$_]} (0, 1, @$order)), "\n";
  }
  close($oexprfh);
  close($oidxfh);
  my @chdrs;
  foreach my $k (keys %{$self->{'survivalHdrs'}}) {
    $chdrs[$self->{'survivalHdrs'}->{$k}] = $k;
  }
  open(my $ocfh, ">$cf") || die "Can't open $cf\n";
  open(my $oifh, ">$if") || die "Can't open $if\n";
  print $oifh "ArrayID\tArrayHeader\tClinicalhHeader\n";
  print $ocfh join("\t", @chdrs), "\n";
  foreach my $i (@$order) {
    my $arr = $headers[$i];
    my @list = map { "" } (@chdrs);
    $list[0] = $arr;
    if (defined $self->{'survivalHash'}->{$arr}) {
      foreach my $k (@chdrs) {
        next if ($k eq "ArrayID");
        my $j1 = $self->{'survivalHdrs'}->{$k};
        $list[$j1] = $self->{'survivalHash'}->{$arr}->[$j1 - 1];
      }
    }
    print $ocfh join("\t", @list), "\n";
    my $j = $self->{'survivalHdrs'}->{"c Title"};
    print $oifh join("\t", $arr, $arr, $list[$j]), "\n";
  }
  close($oifh);
  close($ocfh);
}

sub getObjSurvival {
  my $self = shift;
  my $start = $self->getStart();
  my $end = $self->getEnd();
  my $obj = {};
  $obj->{'hash'} = {};
  $obj->{'headers'} = [ map { $self->{'headers'}->[$_] } $start .. $end ];
  $obj->{'cheaders'} = [];
  $obj->{'narray'} = [];
  $obj->{'carray'} = [];

  foreach my $k (keys %{$self->{'survivalHdrs'}}) {
    my $i = $self->{'survivalHdrs'}->{$k};
    next if ($i < 3);
    my ($c, $v) = split(" ", $k, 2);
    $obj->{'cheaders'}->[$i - 3] = $v;
    if ($c eq "c") {
      push @{$obj->{'carray'}}, $i - 3;
    }
    if ($c eq "n") {
      push @{$obj->{'narray'}}, $i - 3;
    }
  }

  $obj->{'narray'} = [ sort { $a <=> $b } @{$obj->{'narray'}} ];
  $obj->{'carray'} = [ sort { $a <=> $b } @{$obj->{'carray'}} ];

  foreach my $arr (@{$obj->{'headers'}}) {
    my $time = "";
    my $status = "";
    my @l = map { "" } @{$obj->{'cheaders'}};
    if (defined $self->{'survivalHash'}->{$arr}) {
       $time = $self->{'survivalHash'}->{$arr}->[0];
       $status = $self->{'survivalHash'}->{$arr}->[1];
       @l = map { $self->{'survivalHash'}->{$arr}->[$_ + 2] } 0 .. $#{$obj->{'cheaders'}};
    }
    $obj->{'hash'}->{$arr} = [$time, $status, [@l]];
  }

  return $obj;
}

sub getCoxArrays {
  my ($self, %params) = @_;
  my $time = $self->getTimeData();
  my $status = $self->getStatusData();
  return &getCoxArraysR($time, $status, %params);
}

sub getCoxArraysR {
  my ($time, $status, %params) = @_;
  my $res = [];
  for(my $i = 2; $i < scalar(@$time); $i++) {
    next if (!defined $time->[$i] || $time->[$i] eq "");
    next if (!defined $status->[$i] || $status->[$i] eq "");
    my $empty = 0;
    foreach my $k (keys(%params)) {
      if (!defined $params{$k}->[$i] || $params{$k}->[$i] eq "") {
        $empty = 1;
        last;
      }
    }
    next if ($empty == 1);
    push @$res, $i;
  }
  return $res;
}

sub setupArrayListData {
  my ($h, $ofile, $groups, %params) = @_;
  open(my $ofh, ">$ofile") || die "Can't write $ofile\n";
  foreach my $i (0 .. $#{$groups}) {
    foreach my $j (@{$groups->[$i]->[2]}) {
      print $ofh join("\t", $h->{'headers'}->[$j], $i, $groups->[$i]->[0]), "\n";
    }
  }
  close($ofh);
}

sub getJavaDiff {
  my ($h, $ef, $lf, $dir) = @_;
  my $cmd = "java -cp $dir tools.Hegemon diff $ef $lf";
  #print $cmd, "\n";
  my $pid = open2($phr, $phw, $cmd);
  my $start = 0;
  my $res = [];
  while (my $in = <$phr>) {
    $in =~ s/[\r\n]//g;
    my @list = split("\t", $in);
    push @$res, [@list];
  }
  return $res;
}

sub getBoxStats {
  my ($h, $id, $pG, $ph) = @_;
  my $ghash = $h->getGroupsData($id, $pG);
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
    my $xa = $ghash->[$i-1];
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
    my $xa = $ghash->[$i-1];
    if (scalar(@$xa) <= 0) {
      $res .= "ex <- c()\n";
    }
    else {
      $res .= "ex <- x[[$i]]\n";
    }
    for($j= ($i + 1); $j <= $num; $j++){
    my $xb = $ghash->[$j-1];
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
  my ($h, $file, $x_id, $x_name, $pG, %params) = @_;
  my $ghash = $h->getGroupsData($x_id, $pG);
  my $bplots = [];
  foreach my $g (@$ghash) {
    my $res = &U::box_plot_values($g);
    push @$bplots, $res;
  }
  my $binfo = $h->getBoxStats($x_id, $pG);
  my $xldata = [ 0 .. $#{$pG} ];
  if (defined $params{"sort"}) {
    $xldata = [ sort { 
      scalar(@{$ghash->[$a]}) <=> scalar(@{$ghash->[$b]}) 
      } 0 .. $#{$pG} ];
  }
  open(my $ofh, ">$file") || die "Can't open $file\n";
  print $ofh "\\begin{tikzpicture}\n";
  my $index = 1;
  foreach my $i (@$xldata) {
    my $clr = &U::getPScolor("clr$i", $pG->[$i]->[1]);
    print $ofh "$clr\n";
    if (scalar(@{$ghash->[$i]}) > 0) {
      my $c = &U::n2a($i);
      print $ofh "\\pgfplotstableread{%\n";
      foreach my $v (@{$ghash->[$i]}) {
        print $ofh "$index $v\n";
      }
      print $ofh "}\\g$c"."data%\n";
    }
    $index++;
  }

  my $yl = &U::myescape("$x_id: $x_name");
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
    if (defined $params{"points"} && scalar(@{$ghash->[$i]}) > 0) {
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

1;
