# Author: Debashis Sahoo <dsahoo@ucsd.edu>
package Network;

use POSIX;
use Bit::Vector;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my %params = @_;
  my $self  = {};
  $self->{'file'}  = undef;
  $self->{'cachesize'}  = 1000;
  $self->{'blocksize'}  = 5000;
  $self->{'mode'}  = "r";

  $self->{'list'}  = {};
  $self->{'fh'}  = undef;
  $self->{'cache'}  = {};
  $self->{'num'}  = 0;
  $self->{'numBits'}  = 0;
  $self->{'numBytes'}  = 0;
  $self->{'startPtr'} = 0;

  foreach (keys(%params)) {
    $self->{$_} = $params{$_};
  }

  bless ($self, $class);
  return $self;
}

sub init {
  my $self = shift;
  my $file = $self->{'file'};
  if ($self->{'mode'} eq "r") {
    my $fh;
    open($fh, "<$file") || die "Can't open $file\n";
    binmode($fh);
    $self->{'fh'} = $fh;
  }
  if ($self->{'mode'} eq "rw") {
    my $fh;
    open($fh, "+<$file") || die "Can't open $file\n";
    binmode($fh);
    $self->{'fh'} = $fh;
  }
  if ($self->{'mode'} eq "w") {
    my $fh;
    open($fh, ">$file") || die "Can't open $file\n";
    binmode($fh);
    $self->{'fh'} = $fh;
  }
}

sub readHeader {
  my $self = shift;
  my $fh = $self->{'fh'};
  read($fh, $buffer, 1); 
  my $magic = unpack("C", $buffer);
  $self->{'magic'} = $magic; 
  read($fh, $buffer, 1); 
  my $major = unpack("C", $buffer);
  $self->{'major'} = $major; 
  read($fh, $buffer, 1); 
  my $minor = unpack("C", $buffer);
  $self->{'minor'} = $minor; 
}

sub readList {
  my ($self, $num) = @_;
  my $fh = $self->{'fh'};
  my $buffer;
  seek($fh, 3 + $num * 4, 0);
  read($fh, $buffer, 4);
  my $ptr = unpack("I", $buffer);
  seek($fh, $ptr, 0);
  read($fh, $buffer, 4);
  my $len = unpack("I", $buffer);
  read($fh, $buffer, $len);
  my $name = $buffer;
  read($fh, $buffer, 4);
  my $size = unpack("I", $buffer);
  my $res = [];
  for (my $i =0; $i < $size; $i ++) {
    read($fh, $buffer, 4);
    my $len = unpack("I", $buffer);
    read($fh, $buffer, $len);
    my $n = $buffer;
    push @$res, $n;
  }
  $self->{'list'}->{$name} = $res;
  return $res;
}

sub getLow {
  my $self = shift;
  return $self->readList(0);
}

sub getHigh {
  my $self = shift;
  return $self->readList(1);
}

sub getBalanced {
  my $self = shift;
  return $self->readList(2);
}

sub getNum {
  my $self = shift;
  return $self->{'num'};
}

sub getNumBits {
  my $self = shift;
  return $self->{'numBits'};
}

sub getNumBytes {
  my $self = shift;
  return $self->{'numBytes'};
}

sub getStartPtr {
  my $self = shift;
  return $self->{'startPtr'};
}

sub getMatrixEnd {
  my $self = shift;
  return $self->{'startPtr'} + $self->{'num'} * $self->{'numBytes'};
}

sub setMatrixEnd {
  my $self = shift;
  my $ptr = $self->getMatrixEnd();
  my $fh = $self->{'fh'};
  seek($fh, $ptr, 0);
}

sub readNetwork {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $buffer;

  seek($fh, 3 + 3 * 4, 0);
  read($fh, $buffer, 4);
  my $ptr = unpack("I", $buffer);
  read($fh, $buffer, 4);
  my $num = unpack("I", $buffer);
  read($fh, $buffer, 4);
  my $numBits = unpack("I", $buffer);
  my $numBytes = floor($num * $numBits/8) + 1;
  seek($fh, $ptr, 0);
  $self->{'num'} = $num;
  $self->{'numBits'} = $numBits;
  $self->{'numBytes'} = $numBytes;
  $self->{'startPtr'} = $ptr;
  #print "$num, $numBits, $numBytes, $ptr\n";
}

sub readNetworkFile {
  my $self = shift;
  $self->init();
  $self->readHeader();
  $self->readList(0);
  $self->readList(1);
  $self->readList(2);
  $self->readNetwork();
}

sub hashCode {
  my ($self, $a, $b) = @_;
  my $hash = $a * $self->{'numBytes'} + floor($b * $self->{'numBits'}/8);
  return floor($hash/$self->{'blocksize'});
}

sub loadBlock {
  my ($self, $hash)= @_;
  if (scalar(keys(%{$self->{'cache'}})) >= $self->{'cachesize'}) {
    $self->flushCache();
  }
  if (!defined $self->{'cache'}->{$hash}) {
    print STDERR "L $hash\n";
    my $pos = $self->{'startPtr'} + $hash * $self->{'blocksize'};
    my $buffer;
    my $fh = $self->{'fh'};
    seek($fh, $pos, 0);
    read($fh, $buffer, $self->{'blocksize'});
    $self->{'cache'}->{$hash} = $buffer;
  }
}

sub flushCache {
  my $self = shift;
  if ($self->{'mode'} eq "w" || $self->{'mode'} eq "rw") {
    foreach (sort keys(%{$self->{'cache'}})) {
      my $fh = $self->{'fh'};
      my $buffer = $self->{'cache'}->{$_};
      my $pos = $self->{'startPtr'} + $_ * $self->{'blocksize'};
      if ($pos != tell($fh)) {
        seek($fh, $pos, 0);
      }
      print $fh $buffer;
    }
    print STDERR "Cache clear\n";
  }
  $self->{'cache'} = {};
}

sub close {
  my $self = shift;
  $self->flushCache();
  my $fh = $self->{'fh'};
  close($fh);
}

sub writeHeader {
  my $self = shift;
  my $fh = $self->{'fh'};
  $self->{'magic'} = 0x55;
  $self->{'major'} = 1;
  $self->{'minor'} = 0;
  print $fh pack("C", $self->{'magic'});
  print $fh pack("C", $self->{'major'});
  print $fh pack("C", $self->{'minor'});
  for (my $i = 0; $i < 10; $i++) {
    $self->writeInt(0);
  }
}

sub writeHeader_1_1 {
  my $self = shift;
  my $fh = $self->{'fh'};
  $self->{'magic'} = 0x55;
  $self->{'major'} = 1;
  $self->{'minor'} = 1;
  print $fh pack("C", $self->{'magic'});
  print $fh pack("C", $self->{'major'});
  print $fh pack("C", $self->{'minor'});
  for (my $i = 0; $i < 10; $i++) {
    $self->writeInt(0);
  }
}

sub writeInt {
  my ($self, $val) = @_;
  my $fh = $self->{'fh'};
  my $str = pack("L", $val);
  print $fh  $str;
}

sub readInt {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $buffer;
  read($fh, $buffer, 4);
  my $val = unpack("I", $buffer);
  return $val;
}

sub writeLengthPrefixString {
  my ($self, $str) = @_;
  my $fh = $self->{'fh'};
  $self->writeInt(length($str));
  print $fh  $str;
}

sub writeList {
  my ($self, $num, $tag, $list) = @_;
  my $fh = $self->{'fh'};
  my $ptr = tell($fh);
  seek($fh, 3 + $num * 4, 0);
  $self->writeInt($ptr);
  seek($fh, $ptr, 0);
  $self->writeLengthPrefixString($tag);
  my $size = scalar(@{$list});
  $self->writeInt($size);
  for (my $i =0; $i < $size; $i ++) {
    $self->writeLengthPrefixString($list->[$i]);
  }
  $self->{'list'}->{$tag} = $list;
}

sub startMatrix {
  my ($self, $num, $numBits) = @_;
  my $numBytes = floor($num * $numBits/8) + 1;
  my $fh = $self->{'fh'};
  $self->{'num'} = $num;
  $self->{'numBits'} = $numBits;
  $self->{'numBytes'} = $numBytes;
  $self->{'startPtr'} = tell($fh);
  seek($fh, 3 + 3 * 4, 0);
  $self->writeInt($self->{'startPtr'});
  $self->writeInt($self->{'num'});
  $self->writeInt($self->{'numBits'});
  seek($fh, $self->{'startPtr'}, 0);
}

sub writeCode {
  my ($self, $a, $b, $code) = @_;
  my $hash = $self->hashCode($a, $b);
  my $buffer;
  $self->loadBlock($hash);
  $buffer = $self->{'cache'}->{$hash};
  my $numBytes = $self->{'numBytes'};
  my $numBits = $self->{'numBits'};
  my $num = $self->{'num'};
  my $blocksize = $self->{'blocksize'};
  my $cachesize = $self->{'cachesize'};
  my $byte_offset = ($a * $numBytes + floor($b * $numBits/8)) % $blocksize;
  my $bit_offset = ($b * $numBits) % 8;
  my $mask = (1 << $numBits) - 1;
  if (($bit_offset + $numBits) > 8) {
    if (($byte_offset+1) < $blocksize) {
      my $val0 = unpack "c", substr($buffer, $byte_offset, 1);
      $val0 = $val0 & 0xff;
      my $val1 = unpack "c", substr($buffer, $byte_offset+1, 1);
      $val1 = $val1 & 0xff;
      my $val = $val0 | ($val1 << 8);
      $val = (($val & ~($mask << $bit_offset)) | (($code & $mask) << $bit_offset));
      substr($buffer, $byte_offset, 1) = pack "c", ($val & 0xff);
      substr($buffer, $byte_offset+1, 1) = pack "c", ($val >> 8);
    }
    else {
      my $val0 = unpack "c", substr($buffer, $byte_offset, 1);
      $val0 = $val0 & 0xff;
      $self->loadBlock($hash+1);
      $buffer = $self->{'cache'}->{$hash+1};
      my $val1 = unpack "c", substr($buffer, 0, 1);
      $val1 = $val1 & 0xff;
      my $val = $val0 | ($val1 << 8);
      $val = (($val & ~($mask << $bit_offset)) | (($code & $mask) << $bit_offset));
      substr($buffer, 0, 1) = pack "c", ($val >> 8);
      $self->{'cache'}->{$hash+1} = $buffer;
      $self->loadBlock($hash);
      $buffer = $self->{'cache'}->{$hash};
      substr($buffer, $byte_offset, 1) = pack "c", ($val & 0xff);
    }
  }
  else {
    my $val = unpack "c", substr($buffer, $byte_offset, 1);
    $val = $val & 0xff;
    $val = (($val & ~($mask << $bit_offset)) | (($code & $mask) << $bit_offset));
    substr($buffer, $byte_offset, 1) = pack "c", ($val & 0xff);
    #my $pos = $startPtr + $hash * $blocksize;
    #print "$val $code $bit_offset $byte_offset $pos\n";
  }
  $self->{'cache'}->{$hash} = $buffer;
}

sub readCode {
  my ($self, $a, $b) = @_;
  my $hash = $self->hashCode($a, $b);
  my $buffer;
  $self->loadBlock($hash);
  $buffer = $self->{'cache'}->{$hash};
  my $numBytes = $self->{'numBytes'};
  my $numBits = $self->{'numBits'};
  my $num = $self->{'num'};
  my $blocksize = $self->{'blocksize'};
  my $cachesize = $self->{'cachesize'};
  my $byte_offset = ($a * $numBytes + floor($b * $numBits/8)) % $blocksize;
  my $bit_offset = ($b * $numBits) % 8;
  my $mask = (1 << $numBits) - 1;
  my $code = 0;
  if (($bit_offset + $numBits) > 8) {
    if (($byte_offset+1) < $blocksize) {
      my $val0 = unpack "c", substr($buffer, $byte_offset, 1);
      $val0 = $val0 & 0xff;
      my $val1 = unpack "c", substr($buffer, $byte_offset+1, 1);
      $val1 = $val1 & 0xff;
      my $val = $val0 | ($val1 << 8);
      $code = ($val >> $bit_offset) & $mask;
    }
    else {
      my $val0 = unpack "c", substr($buffer, $byte_offset, 1);
      $val0 = $val0 & 0xff;
      $self->loadBlock($hash+1);
      $buffer = $self->{'cache'}->{$hash+1};
      my $val1 = unpack "c", substr($buffer, 0, 1);
      $val1 = $val1 & 0xff;
      my $val = $val0 | ($val1 << 8);
      $code = ($val >> $bit_offset) & $mask;
    }
  }
  else {
    my $val = unpack "c", substr($buffer, $byte_offset, 1);
    $val = $val & 0xff;
    $code = ($val >> $bit_offset) & $mask;
    #my $pos = $startPtr + $hash * $blocksize;
    #print "$val $code $bit_offset $byte_offset $pos\n";
  }
  return $code;
}

sub print_1_0 {
  my $self = shift;
  my $numBytes = $self->{'numBytes'};
  my $numBits = $self->{'numBits'};
  my $num = $self->{'num'};
  my $blocksize = $self->{'blocksize'};
  my $cachesize = $self->{'cachesize'};
  print "Matrix [ n=$num, bit=$numBits, row=$numBytes bytes]\n";
  foreach my $name ("low", "high", "balanced") {
    my $size = scalar(@{$self->{'list'}->{$name}});
    print "$name ($size):";
    for (my $i =0; $i < $size; $i ++) {
      if (($i % 5) == 0) {
        print "\n";
      }
      my $id = $self->{'list'}->{$name}->[$i];
      print $id, ", ";
    }
    print "\n";
  }
  print "Matrix [ n=$num, bit=$numBits, row=$numBytes bytes]\n";
  my $balanced = $self->{'list'}->{"balanced"};
  for (my $i = 0; $i < $self->{'num'}; $i++) {
    print STDERR $i, "\n";
    for (my $j = 0; $j < $self->{'num'}; $j++) {
      my $code = $self->readCode($i, $j);
      if ($code > 0) {
        print "$code\t$i\t$j\t".$balanced->[$i]."\t".$balanced->[$j]."\n";
      }
    }
  }
}

# 0 - No relation
# 1 - $i low -> $j high
# 2 - $i low -> $j low
# 3 - $i high -> $j high
# 4 - $i high -> $j low
# 5 - Equivalent
# 6 - Opposite
sub readCode_1_1 {
  my ($self, $i, $j) = @_;
  my $b0 = $self->readCode(2 * $i,   2 * $j);
  my $b1 = $self->readCode(2 * $i,   2 * $j+1);
  my $b2 = $self->readCode(2 * $i+1, 2 * $j);
  my $b3 = $self->readCode(2 * $i+1, 2 * $j+1);
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

sub hasStats {
  my $self = shift;
  my $buffer;
  my $fh = $self->{'fh'};
  seek($fh, 3 + 6 * 4, 0);
  read($fh, $buffer, 4);
  my $ptr = unpack("I", $buffer);
  return $ptr;
}

sub printStats {
  my $self = shift;
  my $balanced = $self->{'list'}->{"balanced"};
  print "Stats : ", $self->getMatrixEnd(), "\n";
  $self->setMatrixEnd();
  for (my $i = 0; $i < $self->{'num'}/2; $i++) {
    print $balanced->[$i];
    for (my $j = 0; $j < 7; $j++) {
      print "\t".$self->readInt();
    }
    for (my $j = 0; $j < 3; $j++) {
      $self->readInt();
    }
    print "\n";
  }
}

sub print_1_1 {
  my $self = shift;
  my $numBytes = $self->{'numBytes'};
  my $numBits = $self->{'numBits'};
  my $num = $self->{'num'};
  my $blocksize = $self->{'blocksize'};
  my $cachesize = $self->{'cachesize'};
  print "Matrix [ n=$num, bit=$numBits, row=$numBytes bytes]\n";
  foreach my $name ("low", "high", "balanced") {
    my $size = scalar(@{$self->{'list'}->{$name}});
    print "$name ($size):";
    for (my $i =0; $i < $size; $i ++) {
      if (($i % 5) == 0) {
        print "\n";
      }
      my $id = $self->{'list'}->{$name}->[$i];
      print $id, ", ";
    }
    print "\n";
  }
  print "Matrix [ n=$num, bit=$numBits, row=$numBytes bytes]\n";
  my $balanced = $self->{'list'}->{"balanced"};
  for (my $i = 0; $i < $self->{'num'}/2; $i++) {
    print STDERR $i, "\n";
    for (my $j = 0; $j < $self->{'num'}/2; $j++) {
      my $code = $self->readCode_1_1($i, $j);
      if ($code > 0) {
        print "$code\t$i\t$j\t".$balanced->[$i]."\t".$balanced->[$j]."\n";
      }
    }
  }
  if ($self->hasStats() != 0) {
    $self->printStats();
  }
}

sub print {
  my $self = shift;
  my $version = $self->{'major'}.".".$self->{'minor'};
  print "Version : $version\n";
  if ($version eq "1.0") {
    return $self->print_1_0();
  }
  if ($version eq "1.1") {
    return $self->print_1_1();
  }
}

sub readNCode_1_0 {
  my ($self, $network, $i, $j, $vec) = @_;
  my $offset = $self->{'numBytes'}*8 - ($j+1)*$self->{'numBits'};
  $vec->Interval_Copy($network->[$i],0,$offset,$self->{'numBits'});
}

sub writeNCode_1_0 {
  my ($self, $network, $i, $j, $code) = @_;
  my $offset = $self->{'numBytes'}*8 - ($j+1)*$self->{'numBits'};
  $network->[$i]->Interval_Copy($code, $offset, 0, $self->{'numBits'});
}

sub readNCode_1_1 {
  my ($self, $network, $i, $j, $vec) = @_;
  my $len = $self->{'numBits'}*2;
  my $offset = $self->{'num'} - ($j+1)*$len;
  #print "len=$len offset=$offset\n";
  $vec->Interval_Copy($network->[2*$i],0,$offset,$len);
  $vec->Interval_Copy($network->[2*$i+1],$len,$offset,$len);
}

sub writeNCode_1_1 {
  my ($self, $network, $i, $j, $code) = @_;
  my $len = $self->{'numBits'}*2;
  my $offset = $self->{'num'} - ($j+1)*$len;
  $network->[2*$i]->Interval_Copy($code, $offset, 0, $len);
  $network->[2*$i+1]->Interval_Copy($code, $offset, $len, $len);
}

sub readNCode {
  my ($self, $network, $i, $j, $vec) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  #print "version=$version\n";
  if ($version eq "1.0") {
    return $self->readNCode_1_0($network, $i, $j, $vec);
  }
  if ($version eq "1.1") {
    return $self->readNCode_1_1($network, $i, $j, $vec);
  }
}

sub writeNCode {
  my ($self, $network, $i, $j, $code) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    return $self->writeNCode_1_0($network, $i, $j, $code);
  }
  if ($version eq "1.1") {
    return $self->writeNCode_1_1($network, $i, $j, $code);
  }
}

sub contraPositive_1_0 {
  my $code = shift;
  my $c = $code->to_Bin();
  if ($c eq "010") { $code->from_Bin("110"); }
  if ($c eq "110") { $code->from_Bin("010"); }
  return $code;
}

sub contraPositive_1_1 {
  my $code = shift;
  my $len = shift;
  $code->Reverse($code);
  $code->Transpose($len,$len,$code,$len,$len);
  return $code;
}

sub fillLowerTriangle_1_0 {
  my $self = shift;
  my $network = $self->readAll();
  my $code = Bit::Vector->new($self->{'numBits'});
  for (my $i = 0; $i < $self->{'num'}; $i++) {
    print STDERR $i, "\n";
    for (my $j = $i; $j < $self->{'num'}; $j++) {
      $self->readNCode($network, $i, $j, $code);
      if (!$code->is_empty()) {
        $self->writeNCode($network, $j, $i, &contraPositive_1_0($code));
      }
    }
  }
  $self->readNetwork();
  $self->startMatrix($self->{'num'}, $self->{'numBits'});
  $self->writeNetwork($network);
}

sub fillLowerTriangle_1_1 {
  my $self = shift;
  my $network = $self->readAll();
  my $len = $self->{'numBits'}*2;
  my $code = Bit::Vector->new($self->{'numBits'}*2*2);
  for (my $i = 0; $i < $self->{'num'}/2; $i++) {
    print STDERR $i, "\n";
    for (my $j = $i; $j < $self->{'num'}/2; $j++) {
      $self->readNCode($network, $i, $j, $code);
      if (!$code->is_empty()) {
        $self->writeNCode($network, $j, $i, &contraPositive_1_1($code, $len));
      }
    }
  }
  $self->readNetwork();
  $self->startMatrix($self->{'num'}, $self->{'numBits'});
  $self->writeNetwork($network);
}

sub fillLowerTriangle {
  my $self = shift;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    return $self->fillLowerTriangle_1_0();
  }
  if ($version eq "1.1") {
    return $self->fillLowerTriangle_1_1();
  }
}

sub getBitVector {
  my $str = shift;
  my $bitstr = unpack "b*", $str;
  $vec = Bit::Vector->new_Bin(length($bitstr),$bitstr);
  return $vec;
}

sub readMatrix {
  my $self = shift;
  $self->readNetwork();
  my $fh = $self->{'fh'};
  my $numBytes = $self->{'numBytes'};
  my $num = $self->{'num'};
  my $network = [];
  my $buffer;
  for (my $i = 0; $i < $num; $i++) {
    print STDERR $i, "\n";
    read($fh, $buffer, $numBytes);
    if (length($buffer) < $numBytes) {
      my $str = "\0" x ($numBytes - length($buffer));
      $buffer = $buffer.$str;
    }
    $network->[$i] = $buffer;
  }
  return $network;
}

sub convertBitVectorMatrix {
  my ($self, $network) = @_;
  my $num = $self->{'num'};
  for (my $i = 0; $i < $num; $i++) {
    print STDERR $i, "\n";
    $network->[$i] = &getBitVector($network->[$i]);
  }
  return $network;
}

sub readAll_1_0 {
  my $self = shift;
  my $network = $self->readMatrix();
  $self->convertBitVectorMatrix($network);
  return $network;
}

sub readMatrix_1_1 {
  my $self = shift;
  $self->readNetwork();
  my $fh = $self->{'fh'};
  my $numBytes = $self->{'numBytes'};
  my $num = $self->{'num'};
  my $network = [];
  my $buffer;
  for (my $i = 0; $i < $num; $i++) {
    print STDERR $i, "\n";
    read($fh, $buffer, $numBytes);
    $buffer = unpack "b*", $buffer;
    if (length($buffer) < $num) {
      my $str = "0" x ($num - length($buffer));
      $buffer = $buffer.$str;
    }
    else {
      $buffer = substr($buffer, 0, $num);
    }
    my $vec = Bit::Vector->new_Bin(length($buffer),$buffer);
    if ($vec->Size() != $num) {
      print STDERR "Error in reading ", $vec->Size(), " != $num\n";
      exit;
    }
    $network->[$i] = $vec;
  }
  return $network;
}

sub readAll_1_1 {
  my $self = shift;
  my $network = $self->readMatrix_1_1();
  return $network;
}

sub readAll {
  my $self = shift;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    return $self->readAll_1_0();
  }
  if ($version eq "1.1") {
    return $self->readAll_1_1();
  }
}

sub writeNetwork_1_0 {
  my ($self, $network) = @_;
  my $fh = $self->{'fh'};
  my $num = $self->{'num'};
  if ($num != scalar(@{$network})) {
    print STDERR "Number of element in the network and the balanced list are different";
    print STDERR "$num != ".scalar(@{$network})."\n";
    exit(1);
  }
  my $i = 0;
  foreach my $row (@{$network}) {
    print STDERR "$i\n";
    my $str = $row->to_Bin();
    my $len = length($str);
    if ($len > ($self->{'numBytes'} * 8)) {
      $str = substr($str, 0, ($self->{'numBytes'} * 8));
    }
    else {
      my $remaining = ($self->{'numBytes'} * 8) - $len;
      $str = $str . ("0" x $remaining);
    }
    $len = length($str);
    my $buffer = pack("b$len", $str);
    print $fh $buffer;
    $i++;
  }
}

sub writeNetwork_1_1 {
  my ($self, $network) = @_;
  my $fh = $self->{'fh'};
  my $num = $self->{'num'};
  if ($num != scalar(@{$network})) {
    print STDERR "Number of element in the network and the balanced list are different";
    print STDERR "$num != ".scalar(@{$network})."\n";
    exit(1);
  }
  my $i = 0;
  foreach my $row (@{$network}) {
    print STDERR "$i\n";
    my $str = $row->to_Bin();
    my $len = length($str);
    if ($len > ($self->{'numBytes'} * 8)) {
      $str = substr($str, 0, ($self->{'numBytes'} * 8));
    }
    else {
      my $remaining = ($self->{'numBytes'} * 8) - $len;
      $str = $str . ("0" x $remaining);
    }
    $len = length($str);
    my $buffer = pack("b$len", $str);
    print $fh $buffer;
    $i++;
  }
}

sub writeNetwork {
  my ($self, $network) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    return $self->writeNetwork_1_0($network);
  }
  if ($version eq "1.1") {
    return $self->writeNetwork_1_1($network);
  }
}

sub write_1_0 {
  my ($self, $low, $high, $balanced, $network, $numBits) = @_;
  $self->init();
  $self->writeHeader();
  $self->writeList(0, "low", $low);
  $self->writeList(1, "high", $high);
  $self->writeList(2, "balanced", $balanced);
  my $num = scalar(@{$balanced});
  $self->startMatrix($num, $numBits);
  $self->writeNetwork($network);
}

sub write_1_1 {
  my ($self, $low, $high, $balanced, $network, $numBits) = @_;
  $self->init();
  $self->writeHeader_1_1();
  $self->writeList(0, "low", $low);
  $self->writeList(1, "high", $high);
  $self->writeList(2, "balanced", $balanced);
  my $num = scalar(@{$balanced});
  $self->startMatrix($num*2, $numBits);
  $self->writeNetwork($network);
}

sub write {
  my ($self, $low, $high, $balanced, $network, $numBits) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    return $self->write_1_0($low, $high, $balanced, $network, $numBits);
  }
  if ($version eq "1.1") {
    return $self->write_1_1($low, $high, $balanced, $network, $numBits);
  }
  return $self->write_1_0($low, $high, $balanced, $network, $numBits);
}

sub readIndex {
  my $indexFile = shift;
  my $idhash = {} ;
  my $namehash = {} ;
  my $revidhash = {} ;
  my $revnamehash = {} ;
  my $revptrhash = {} ;
  open(FL, "<$indexFile") || die "Can't open $indexFile\n";
  <FL>; # Header
  my $index = 0;
  while (<FL>) {
    my ($id, $name, $ptr, $desc) = split("\t", $_);
    my $idu = uc($id);
    my $nameu = uc($name);
    if (!defined $idhash->{$idu}) {
      $idhash->{$idu} = [];
    }
    push @{$idhash->{$idu}}, $index;
    if (!defined $namehash->{$nameu}) {
      $namehash->{$nameu} = [];
    }
    push @{$namehash->{$nameu}}, $index;
    $revptrhash->{$index} = $ptr;
    $revidhash->{$index} = $id;
    $revnamehash->{$index} = $name;
    $index++;
  }
  close(FL);
  $indexFile =~ s/.idx$/.thr/g;
  my $thrHash = &readThrFile($indexFile);
  $indexFile =~ s/.thr$/.pcl/g;
  return [$idhash, $namehash, $revidhash, $revnamehash, $revptrhash, $thrHash,
         $indexFile];
}

sub getLink {
  my ($gene, $org) = @_;
  if ($org eq "Mm") {
    return "http://smd.stanford.edu/cgi-bin/source/sourceResult?criteria=$gene&choice=Gene&option=Name&organism=Mm";
  }
  elsif ($org eq "Dm") {
    return "http://flybase.org/cgi-bin/uniq.html?context=$gene&species=Dmel&db=fbgn&caller=quicksearch";
  }
  elsif ($org eq "Sgd") {
    return "http://www.genedb.org/genedb/Dispatcher?formType=navBar&organism=cerevisiae&name=$gene&desc=yes&submit=Search";
  }
  elsif ($org eq "Ath") {
    return "http://www.arabidopsis.org/servlets/TairObject?name=$gene&type=gene";
  }
  elsif ($org eq "Affy") {
    return "https://www.affymetrix.com/LinkServlet?probeset=$gene";
  }
  else {
    return "http://smd.stanford.edu/cgi-bin/source/sourceResult?criteria=$gene&choice=Gene&option=Name&organism=Hs";
  }
}

sub getSymbol {
  my ($index, $hashTable, $balanced) = @_;
  my $low = 0;
  if ($index =~ /b/) {
    $low = 1;
  }
  $index =~ s/b//g;
  my $probeid = $balanced->[$index];
  my $n = $hashTable->[0]->{uc($probeid)}->[0];
  my $name = $hashTable->[3]->{$n};
  if ($low == 0) {
    return "$name\-$probeid";
  }
  else {
    return "$name\_b-$probeid";
  }
}

sub readThrFile {
  my $thrFile = shift;
  my $hash = {};
  open(FL, "< $thrFile") || die "Can't open $thrFile\n";
  while (<FL>) {
    s/[\r\n]//g;
    my ($index, $thr1, $stat, $thr0, $thr2, @rest) = split("\t");
    $hash->{$index} = [$thr1, $stat, $thr0, $thr2];
  }
  return $hash;
}

sub printListHtml {
  my ($list, $n, $title, $org, $hashTable, $listid, $fh) = @_;
  if (scalar(@{$list}) <= 0) {
    return;
  }
  if (!defined $fh) {
    $fh = \*STDOUT;
  }
  my $symhash = $hashTable->[3];
  my $filePtrHash = $hashTable->[4];
  my $thrHash = $hashTable->[5];
  my $idhash = $hashTable->[2];
  my $file = $hashTable->[6];
  my $phstr = $hashTable->[7];
  my $hash = {};
  for (my $i = 0; $i < scalar(@{$list}); $i++) {
    my $probeid = uc($listid->[$list->[$i]]);
    my $nid = $hashTable->[0]->{$probeid}->[0];
    my $name = $hashTable->[3]->{$nid};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    push @{$hash->{$name}}, $nid;
  }
  print $fh "<li> <font size=+1> $title </font> (",scalar(@{$list}),")<ul>\n";
  print $fh "<li> <table border=0><tr>\n";
  my $num = 5;
  my $index = 0;
  foreach my $gene (sort keys(%{$hash})) {
    foreach my $g (@{$hash->{$gene}}) {
      next if ($gene eq "---");
      next if ($gene =~ /^\s*$/);
      $gene =~ s/\/.*//g;
      if  ( ($index % $num) == 0 ) {
        print $fh "</tr><tr>\n";
      }
      my $name = $gene;
      if (length($name) > 15) {
        $name = substr($name, 0, 6)."...".substr($name,length($name)-6, 6);
      }
      my $link = &getLink($name, $org);
      print $fh "<td> <a target=\"_blank\" href=\"$link\"> $name </a></td>\n";
      my $m = $g;
      my $id = $idhash->{$m};
      my $x = $filePtrHash->{$n};
      my $y = $filePtrHash->{$m};
      my ($thrx1, $stat, $thrx0, $thrx2) = @{$thrHash->{$n}};
      my ($thry1, $stat, $thry0, $thry2) = @{$thrHash->{$m}};
      my $idname = $id;
      my $idlink = &getLink($idname, "Affy");
      if (length($idname) > 15) {
        $idname = substr($idname, 0, 6)."...".substr($idname,length($idname)-6 , 6);   
      }
      print $fh "<td><a target=\"_blank\" href=\"http://hegemon.ucsd.edu/microarray/test/plot.php?id=2&file=$file&x=$x&y=$y&thrx0=$thrx0&thrx1=$thrx1&thrx2=$thrx2&thry0=$thry0&thry1=$thry1&thry2=$thry2&$phstr\"> <img src=\"http://hegemon.ucsd.edu/microarray/test/sp.gif\" border=0> </img> </a>";
      print $fh "<a target=\"_blank\" href=\"$idlink\"> $idname </a> </td>\n";
      $index++;
    }
  }
  print $fh "</tr></table></li>\n";
  print $fh "</ul></li>\n";
}

sub printHTML_1_0 {
  my ($self, $indexFile, $org, $probeid) = @_;
  $probeid = uc($probeid);
  my $file = $self->{'file'};
  my $hashTable = &readIndex($indexFile);
  if (!defined $hashTable->[0]->{$probeid}) {
    print "Cannot fine probe $probeid\n";
    return;
  }
  my $low = $self->getLow();
  my $high = $self->getHigh();
  my $balanced = $self->getBalanced();

  my $lowhash = {};
  my $highhash = {};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$low}); $i++) {
    $lowhash->{$low->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$high}); $i++) {
    $highhash->{$high->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }

print <<END;
<html>
<head>
<SCRIPT SRC="http://genedesk.ucsd.edu/home/public/mktree.js"
LANGUAGE="JavaScript"></SCRIPT>
<LINK REL="stylesheet"
HREF="http://genedesk.ucsd.edu/home/public/mktree.css">
</head>
<body>
<center>
<h1> Boolean interactions of genes </h1>
<br>
<A href="#" onClick="expandTree('tree1'); return
false;">ExpandAll</A>&nbsp;&nbsp;&nbsp;
<A href="#" onClick="collapseTree('tree1'); return
false;">CollapseAll</A>&nbsp;&nbsp;&nbsp;
<A href="#" onClick="expandTreeDepth('tree1',2); return false;">ExpandDepth 2
</A>

</center>
<ul class="mktree" id="tree1">
END

  my $index = $hashTable->[0]->{$probeid}->[0];
  my $name = $hashTable->[3]->{$index};
  my $id = $hashTable->[2]->{$index};
  my $namelink = &getLink($name, $org);
  my $idlink = &getLink($id, "Affy");
  if (defined $lowhash->{$id}) {
    print "<li> <font size=+2> <a href=\"$namelink\"> $name </a> </font> - <a href=\"$idlink\"> $id </a> (Always Low)\n";
    print "</body> </html>\n";
    return;
  }
  if (defined $highhash->{$id}) {
    print "<li> <font size=+2> <a href=\"$namelink\"> $name </a> </font> - <a href=\"$idlink\"> $id </a> (Always High)\n";
    print "</body> </html>\n";
    return;
  }
  if (!defined $balancedhash->{$id}) {
    print "<li> <font size=+2> <a href=\"$namelink\"> $name </a> </font> - <a href=\"$idlink\"> $id </a> (Bad dynamic range)\n";
    print "</body> </html>\n";
    return;
  }
  my $ptr = $self->getStartPtr();
  my $num = $self->getNum();
  my $numBits = $self->getNumBits();
  my $numBytes = $self->getNumBytes();
  my $relations = {};
  my $numrel = 0;
  my $i = $balanced->{$id};
  #print "$num $numBits $numBytes $ptr ".$balanced->{$id}."\n";
  for (my $j = 0; $j < $num; $j++) {
    my $code = $self->readCode($i, $j);
    if ($code > 0) {
      if (!defined $relations->{$code}) {
        $relations->{$code} = [];
      }
      push @{$relations->{$code}}, $j;
      #print "$code ".$balanced->{$id}." $i\n";
      $numrel++;
    }
  }
  print "<li> <font size=+2> <a target=\"_blank\" href=\"$namelink\"> $name </a> </font> - <a target=\"_blank\" href=\"$idlink\"> $id </a> (", $numrel ,")<ul>\n";
  &printListHtml($relations->{5}, $index, "Equivalent", $org, $hashTable, $balanced);
  &printListHtml($relations->{6}, $index, "Opposite", $org, $hashTable, $balanced);
  &printListHtml($relations->{2}, $index, "($name low -> B low)", $org, $hashTable, $balanced);
  &printListHtml($relations->{1}, $index, "($name low -> B high)", $org, $hashTable, $balanced);
  &printListHtml($relations->{4}, $index, "($name high -> B low)", $org, $hashTable, $balanced);
  &printListHtml($relations->{3}, $index, "($name high -> B high)", $org, $hashTable, $balanced);
  print "</ul> </body> </html>\n";
}

sub printHTML_1_1 {
  my ($self, $indexFile, $org, $probeid, $phstr) = @_;
  $probeid = uc($probeid);
  my $file = $self->{'file'};
  my $hashTable = &readIndex($indexFile);
  push @{$hashTable}, $phstr;
  if (!defined $hashTable->[0]->{$probeid}) {
    print "Cannot fine probe $probeid\n";
    return;
  }
  my $low = $self->getLow();
  my $high = $self->getHigh();
  my $balanced = $self->getBalanced();

  my $lowhash = {};
  my $highhash = {};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$low}); $i++) {
    $lowhash->{$low->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$high}); $i++) {
    $highhash->{$high->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }

print <<END;
<html>
<head>
<SCRIPT SRC="http://genedesk.ucsd.edu/home/public/mktree.js"
LANGUAGE="JavaScript"></SCRIPT>
<LINK REL="stylesheet"
HREF="http://genedesk.ucsd.edu/home/public/mktree.css">
</head>
<body>
<center>
<h1> Boolean interactions of genes </h1>
<br>
<A href="#" onClick="expandTree('tree1'); return
false;">ExpandAll</A>&nbsp;&nbsp;&nbsp;
<A href="#" onClick="collapseTree('tree1'); return
false;">CollapseAll</A>&nbsp;&nbsp;&nbsp;
<A href="#" onClick="expandTreeDepth('tree1',2); return false;">ExpandDepth 2
</A>

</center>
<ul class="mktree" id="tree1">
END

  my $index = $hashTable->[0]->{$probeid}->[0];
  my $name = $hashTable->[3]->{$index};
  my $id = $hashTable->[2]->{$index};
  my $cname = &getSplitName($name, 20);
  my $genename = $name;
  $genename =~ s/\/.*//g;
  my $namelink = &getLink($genename, $org);
  my $idlink = &getLink($id, "Affy");
  if (defined $lowhash->{$id}) {
    print "<li> <font size=+2> <a href=\"$namelink\"> $genename </a> </font> - <a href=\"$idlink\"> $id </a> (Always Low)\n";
    print "</body> </html>\n";
    return;
  }
  if (defined $highhash->{$id}) {
    print "<li> <font size=+2> <a href=\"$namelink\"> $genename </a> </font> - <a href=\"$idlink\"> $id </a> (Always High)\n";
    print "</body> </html>\n";
    return;
  }
  if (!defined $balancedhash->{$id}) {
    print "<li> <font size=+2> <a href=\"$namelink\"> $genename </a> </font> - <a href=\"$idlink\"> $id </a> (Bad dynamic range)\n";
    print "</body> </html>\n";
    return;
  }
  my $ptr = $self->getStartPtr();
  my $num = $self->getNum();
  my $numBits = $self->getNumBits();
  my $numBytes = $self->getNumBytes();
  my $relations = {};
  my $numrel = 0;
  my $i = $balancedhash->{$id};
  #print "$num $numBits $numBytes $ptr ".$balancedhash->{$id}."\n";
  for (my $j = 0; $j < $num/2; $j++) {
    my $code = $self->readCode_1_1($i, $j);
    if ($code > 0) {
      if (!defined $relations->{$code}) {
        $relations->{$code} = [];
      }
      push @{$relations->{$code}}, $j;
      #print "$code ".$balancedhash->{$id}." $i\n";
      $numrel++;
    }
  }
  print "<li> <font size=+2> <a target=\"_blank\" href=\"$namelink\"> $genename </a> </font> - <a target=\"_blank\" href=\"$idlink\"> $id </a> (", $numrel ,")<ul>\n";
  &printListHtml($relations->{5}, $index, "Equivalent", $org, $hashTable, $balanced);
  &printListHtml($relations->{6}, $index, "Opposite", $org, $hashTable, $balanced);
  &printListHtml($relations->{2}, $index, "($genename low -> B low)", $org, $hashTable, $balanced);
  &printListHtml($relations->{1}, $index, "($genename low -> B high)", $org, $hashTable, $balanced);
  &printListHtml($relations->{4}, $index, "($genename high -> B low)", $org, $hashTable, $balanced);
  &printListHtml($relations->{3}, $index, "($genename high -> B high)", $org, $hashTable, $balanced);
  print "</ul> </body> </html>\n";
}

sub printHTML {
  my ($self, $indexFile, $org, $probeid, $ph) = @_;
  $probeid = uc($probeid);
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    return $self->printHTML_1_0($indexFile, $org, $probeid);
  }
  if ($version eq "1.1") {
    return $self->printHTML_1_1($indexFile, $org, $probeid, $ph);
  }
}

sub printIndex {
  my ($self, $indexFile) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  my $hashTable = &readIndex($indexFile);
  my $hash = {};
  my $balanced = $self->{'list'}->{"balanced"};
  my $low = $self->{'list'}->{"low"};
  my $high = $self->{'list'}->{"high"};
  for (my $i = 0; $i < scalar(@{$balanced}); $i++) {
    my $probeid = uc($balanced->[$i]);
    my $n = $hashTable->[0]->{$probeid}->[0];
    my $name = $hashTable->[3]->{$n};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    my $chr = substr($name, 0, 1);
    if (!defined $hash->{$chr}) {
        $hash->{$chr} = 0;
    }
    $hash->{$chr} ++;
  }
  my $lownum = scalar(@{$low});
  my $highnum = scalar(@{$high});
  my $balancednum = scalar(@{$balanced});
  print "low\t$lownum\n";
  print "high\t$highnum\n";
  print "balanced\t$balancednum\n";
  foreach (sort keys(%{$hash})) {
    print "$_\t".$hash->{$_}."\n";
  }
}

sub printGeneListAndLinks {
  my ($self, $indexFile, $tag, $org) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  my $hashTable = &readIndex($indexFile);
  my $list = $self->{'list'}->{$tag};
  my $filePtrHash = $hashTable->[4];
  my $thrHash = $hashTable->[5];
  my $idhash = $hashTable->[2];
  my $file = $hashTable->[6];
  my $hash = {};
  for (my $i = 0; $i < scalar(@{$list}); $i++) {
    my $probeid = uc($list->[$i]);
    my $n = $hashTable->[0]->{$probeid}->[0];
    my $name = $hashTable->[3]->{$n};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    push @{$hash->{$name}}, $n;
  }
  foreach my $name (sort keys(%{$hash})) {
    foreach my $n (@{$hash->{$name}}) {
      my $id = $idhash->{$n};
      my $ilink = &getLink($id, "Affy");
      my $link = &getLink($name, $org);
      my $x = $filePtrHash->{$n};
      my ($thrx1, $stat, $thrx0, $thrx2) = @{$thrHash->{$n}};
      my $plink = "http://hegemon.ucsd.edu/microarray/test/plot.php?id=8&file=$file&x=$x&thrx0=$thrx0&thrx1=$thrx1&thrx2=$thrx2";
      print "$name\t$id\t$link\t$ilink\t$plink\n";
    }
  }
}

sub printGeneAListAndLinks {
  my ($self, $indexFile, $tag, $org) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    print "Only Relation format 1.1 is supported\n";
    exit;
  }
  my $hashTable = &readIndex($indexFile);
  my $balanced = $self->{'list'}->{"balanced"};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{uc($balanced->[$i])} = $i;
  }
  $self->printGeneTagListAndLinks($hashTable, $balancedhash, $tag, $org);
}

sub printGeneTagListAndLinks {
  my ($self, $hashTable, $listHash, $tag, $org) = @_;
  my $filePtrHash = $hashTable->[4];
  my $thrHash = $hashTable->[5];
  my $idhash = $hashTable->[2];
  my $file = $hashTable->[6];
  my $hash = {};
  $self->setMatrixEnd();
  my $list = $self->{'list'}->{"balanced"};
  for (my $i = 0; $i < scalar(@{$list}); $i++) {
    my @rels;
    $self->readInt();
    for (my $j = 1; $j < 7; $j++) {
      push @rels, $self->readInt();
    }
    for (my $j = 0; $j < 3; $j++) {
      $self->readInt();
    }
    my $probeid = uc($list->[$i]);
    next if (!defined $listHash->{$probeid});
    my $n = $hashTable->[0]->{$probeid}->[0];
    my $name = $hashTable->[3]->{$n};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    my $chr = substr($name, 0, 1);
    if ($chr eq $tag) {
      push @{$hash->{$name}}, [$n, @rels];
    }
    elsif ($tag eq "All") {
      push @{$hash->{$name}}, [$n, @rels];
    }
  }
  foreach my $name (sort keys(%{$hash})) {
    foreach my $rel (@{$hash->{$name}}) {
      my $n = $rel->[0];
      my $id = $idhash->{$n};
      my $ilink = &getLink($id, "Affy");
      my $link = &getLink($name, $org);
      print "$name\t$id\t";
      for (my $j = 1; $j < 7; $j++) {
        print $rel->[$j]."\t";
      }
      print "$link\t$ilink\n";
    }
  }
}

sub printSubnet {
  my ($self, $indexFile, $list, $org) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    print "Only Relation format 1.1 is supported\n";
    exit;
  }
  my $hashTable = &readIndex($indexFile);
  my $genes = {};
  my $fh;
  open($fh, "<$list") || die "Can't open $list\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split(/[^\w-]+/);
    foreach my $w (@list) {
      $w = uc($w);
      next if ($w =~ /^-*$/);
      next if ($w eq "UNKNOWN");
      if (defined $hashTable->[0]->{$w}) {
        foreach my $n (@{$hashTable->[0]->{$w}}) {
          $genes->{$n} = 1;
        }
      }
      if (defined $hashTable->[1]->{$w}) {
        foreach my $n (@{$hashTable->[1]->{$w}}) {
          $genes->{$n} = 1;
        }
      }
    }
  }
  close($fh);
  my $low = $self->getLow();
  my $high = $self->getHigh();
  my $balanced = $self->getBalanced();

  my $lowhash = {};
  my $highhash = {};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$low}); $i++) {
    $lowhash->{$low->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$high}); $i++) {
    $highhash->{$high->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }
  my $lowgenes = {};
  my $highgenes = {};
  my $badgenes = {};
  my $listHash = {};
  my $lHash = {};
  my ($numbal, $numlow, $numhigh, $numbad) = (0, 0, 0, 0);
  foreach my $k (keys(%{$genes})) {
    my $probeid = $hashTable->[2]->{$k};
    my $name = $hashTable->[3]->{$k};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    if (defined $lowhash->{$probeid}) {
      push @{$lowgenes->{$name}}, $k;
      $numlow++;
    }
    elsif (defined $highhash->{$probeid}) {
      push @{$highgenes->{$name}}, $k;
      $numhigh++
    }
    elsif (defined $balancedhash->{$probeid}) {
      $lHash->{uc($probeid)} = 1;
      push @{$listHash->{$name}}, $k;
      $numbal++;
    }
    else {
      push @{$badgenes->{$name}}, $k;
      $numbad++;
    }
  }
  print "balanced\t$numbal\n";
  $self->printGeneTagListAndLinks($hashTable, $lHash, "All", $org);
  print "high\t$numhigh\n";
  foreach my $name (sort keys(%{$highgenes})) {
    foreach my $index (@{$highgenes->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $ilink = &getLink($probeid, "Affy");
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $link = &getLink($genename, $org);
      my $cname = &getSplitName($name, 20);
      print "$name\t$probeid\t$link\t$ilink\n";
    }
  }
  print "low\t$numlow\n";
  foreach my $name (sort keys(%{$lowgenes})) {
    foreach my $index (@{$lowgenes->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $ilink = &getLink($probeid, "Affy");
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $link = &getLink($genename, $org);
      my $cname = &getSplitName($name, 20);
      print "$name\t$probeid\t$link\t$ilink\n";
    }
  }
  print "bad\t$numbad\n";
  foreach my $name (sort keys(%{$badgenes})) {
    foreach my $index (@{$badgenes->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $ilink = &getLink($probeid, "Affy");
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $link = &getLink($genename, $org);
      my $cname = &getSplitName($name, 20);
      print "$name\t$probeid\t$link\t$ilink\n";
    }
  }
}

sub getSplitName {
  my ($name, $len) = @_;
  my $cname = $name;
  if (length($name) > $len) {
    $cname = "";
    for ($i = 0; $i < length($name); $i+=$len) {
      $cname = $cname . substr($name, $i, $len) . "<br/>";
    }
  }
  return $cname;
}

sub printSubgraph {
  my ($self, $indexFile, $list, $org, $outdir) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    print "Only Relation format 1.1 is supported\n";
    exit;
  }
  print "Creating directory ...\n";
  if (! -e $outdir) {
    mkdir($outdir) || die "Can't create $outdir\n";
  }
  print "Reading Symbols ...\n";
  my $hashTable = &readIndex($indexFile);
  print "Done. Time=".times()."\n";
  print "Reading list ...\n";
  my $genes = {};
  my $fh;
  open($fh, "<$list") || die "Can't open $list\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split(/[^\w-]+/);
    foreach my $w (@list) {
      $w = uc($w);
      next if ($w =~ /^-*$/);
      next if ($w eq "UNKNOWN");
      if (defined $hashTable->[0]->{$w}) {
        foreach my $n (@{$hashTable->[0]->{$w}}) {
          $genes->{$n} = 1;
        }
      }
      if (defined $hashTable->[1]->{$w}) {
        foreach my $n (@{$hashTable->[1]->{$w}}) {
          $genes->{$n} = 1;
        }
      }
    }
  }
  close($fh);
  my $low = $self->getLow();
  my $high = $self->getHigh();
  my $balanced = $self->getBalanced();

  my $lowhash = {};
  my $highhash = {};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$low}); $i++) {
    $lowhash->{$low->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$high}); $i++) {
    $highhash->{$high->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }
  my $lowgenes = {};
  my $highgenes = {};
  my $badgenes = {};
  my $listHash = {};
  my ($numbal, $numlow, $numhigh, $numbad) = (0, 0, 0, 0);
  foreach my $k (keys(%{$genes})) {
    my $probeid = $hashTable->[2]->{$k};
    my $name = $hashTable->[3]->{$k};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    if (defined $lowhash->{$probeid}) {
      push @{$lowgenes->{$name}}, $k;
      $numlow++;
    }
    elsif (defined $highhash->{$probeid}) {
      push @{$highgenes->{$name}}, $k;
      $numhigh++
    }
    elsif (defined $balancedhash->{$probeid}) {
      push @{$listHash->{$name}}, $k;
      $numbal++;
    }
    else {
      push @{$badgenes->{$name}}, $k;
      $numbad++;
    }
  }
  print "Done. Time=".times()."\n";
  print "Creating HTML files ... \n";
  my $indexfh;
  my $genesfh;
  my $gmlfh;
  open($indexfh, ">$outdir/index.html") || die "Can't open $outdir/index.html\n";
  open($genesfh, ">$outdir/genes.txt") || die "Can't open $outdir/genes.txt\n";
  open($gmlfh, ">$outdir/graph.gml") || die "Can't open $outdir/graph.gml\n";
  print $indexfh "<html> <body> <center>\n";
  print $indexfh "<h1>  Boolean interactions of genes </h1>\n";
  print $indexfh "<a href=\"genes.txt\"> List of genes </a> <br/>\n";
  print $indexfh "<a href=\"graph.gml\"> Graph in GML format</a> <br/>\n";

  print $indexfh "<table border=0>";
  print $indexfh "<tr> <td> Affy ID </td>\n";
  print $indexfh "<td> Total </td> <td> lohi </td>\n";
  print $indexfh "<td> lolo </td> <td> hihi </td> <td> hilo </td>\n";
  print $indexfh "<td> equ </td> <td> opo </td> <td> link </td> <td>Symbol</td></tr>\n";
  print $genesfh "balanced\t$numbal\n";
  print $genesfh "Affy ID\tTotal\tlohi\tlolo\thihi\thilo\tequ\topo\tSymbol\n";
  foreach my $name (sort keys(%{$listHash})) {
    foreach my $index (@{$listHash->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $relations = {};
      my $numrel = 0;
      my $i = $balancedhash->{$probeid};
  #print "$num $numBits $numBytes $ptr ".$balancedhash->{$id}."\n";
      foreach my $name1 (sort keys(%{$listHash})) {
        foreach my $index1 (@{$listHash->{$name1}}) {
          my $probeid1 = $hashTable->[2]->{$index1};
          my $j = $balancedhash->{$probeid1};
          my $code = $self->readCode_1_1($i, $j);
          if ($code > 0) {
            if (!defined $relations->{$code}) {
              $relations->{$code} = [];
            }
            push @{$relations->{$code}}, $j;
  #print "$code ".$balancedhash->{$id}." $i\n";
            $numrel++;
          }
        }
      }
      print "$probeid -> $numrel\n";
      next if ($numrel <= 0);
      my $cname = &getSplitName($name, 20);
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $namelink = &getLink($genename, $org);
      my $idlink = &getLink($probeid, "Affy");
      my $fileprobeid = $probeid;
      $fileprobeid =~ s/\//_/g;
      $fileprobeid =~ s/\:/_/g;
      print $genesfh "$probeid\t";
      print $genesfh "$numrel\t";
      print $genesfh scalar(@{$relations->{1}})+0,"\t";
      print $genesfh scalar(@{$relations->{2}})+0,"\t";
      print $genesfh scalar(@{$relations->{3}})+0,"\t";
      print $genesfh scalar(@{$relations->{4}})+0,"\t";
      print $genesfh scalar(@{$relations->{5}})+0,"\t";
      print $genesfh scalar(@{$relations->{6}})+0,"\t";
      print $genesfh "$name\n";
      print $indexfh "<tr>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$idlink\"> $probeid </a> </td>";
      print $indexfh "<td> $numrel </td>";
      print $indexfh "<td>", scalar(@{$relations->{1}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{2}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{3}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{4}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{5}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{6}})+0,"</td>";
      print $indexfh "<td> <a href=\"$fileprobeid\.html\"> link </a> </td>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$namelink\"> $cname </a> </td>";
      print $indexfh "</tr>";
      my $idfh;
      open($idfh, ">$outdir/$fileprobeid\.html") || die "Can't open $outdir/$fileprobeid\.html\n";
print $idfh <<END;
<html>
<head>
<SCRIPT SRC="http://genedesk.ucsd.edu/home/public/mktree.js"
LANGUAGE="JavaScript"></SCRIPT>
<LINK REL="stylesheet"
HREF="http://genedesk.ucsd.edu/home/public/mktree.css">
</head>
<body>
<center>
<h1> Boolean interactions of genes </h1>
<br>
<A href="#" onClick="expandTree('tree1'); return
false;">ExpandAll</A>&nbsp;&nbsp;&nbsp;
<A href="#" onClick="collapseTree('tree1'); return
false;">CollapseAll</A>&nbsp;&nbsp;&nbsp;
<A href="#" onClick="expandTreeDepth('tree1',2); return false;">ExpandDepth 2
</A>

</center>
<ul class="mktree" id="tree1">
END
      print $idfh "<li> <font size=+2> <a target=\"_blank\" href=\"$namelink\"> $genename </a> </font> - <a target=\"_blank\" href=\"$idlink\"> $probeid </a> (", $numrel ,")<ul>\n";
      &printListHtml($relations->{5}, $index, "Equivalent", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{6}, $index, "Opposite", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{2}, $index, "($genename low -> B low)", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{1}, $index, "($genename low -> B high)", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{4}, $index, "($genename high -> B low)", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{3}, $index, "($genename high -> B high)", $org, $hashTable, $balanced, $idfh);
      print $idfh "</ul> </body> </html>\n";
    }
  }
  print $indexfh "</table><br/>";

  print $indexfh "<b> High Genes </b>($numhigh)<br/>\n";
  print $indexfh "<table border=0>\n";
  print $genesfh "high\t$numhigh\n";
  foreach my $name (sort keys(%{$highgenes})) {
    foreach my $index (@{$highgenes->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $ilink = &getLink($probeid, "Affy");
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $link = &getLink($genename, $org);
      my $cname = &getSplitName($name, 20);
      print $genesfh "$probeid\t$name\n";
      print $indexfh "<tr> <td> <a target=\"_blank\" href=\"$ilink\"> $probeid </a> </td>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$link\"> $cname </a> </td> </tr>\n";
    }
  }
  print $indexfh "</table><br/>\n";
  print $indexfh "<b> Low Genes </b>($numlow)<br/>\n";
  print $indexfh "<table border=0>\n";
  print $genesfh "low\t$numlow\n";
  foreach my $name (sort keys(%{$lowgenes})) {
    foreach my $index (@{$lowgenes->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $ilink = &getLink($probeid, "Affy");
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $link = &getLink($genename, $org);
      my $cname = &getSplitName($name, 20);
      print $genesfh "$probeid\t$name\n";
      print $indexfh "<tr> <td> <a target=\"_blank\" href=\"$ilink\"> $probeid </a> </td>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$link\"> $cname </a> </td> </tr>\n";
    }
  }
  print $indexfh "</table><br/>\n";
  print $indexfh "<b> Bad Dynamic Range Genes </b>($numbad)\n";
  print $indexfh "<table border=0>\n";
  print $genesfh "bad\t$numbad\n";
  foreach my $name (sort keys(%{$badgenes})) {
    foreach my $index (@{$badgenes->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $ilink = &getLink($probeid, "Affy");
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $link = &getLink($genename, $org);
      my $cname = &getSplitName($name, 20);
      print $genesfh "$probeid\t$name\n";
      print $indexfh "<tr> <td> <a target=\"_blank\" href=\"$ilink\"> $probeid </a> </td>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$link\"> $cname </a> </td> </tr>\n";
    }
  }
  print $indexfh "</table><br/>\n";
  print $indexfh "</center></body></html>\n";
  close($indexfh);
  close($genesfh);
  close($gmlfh);
  print "Done. Time=".times()."\n";
}

sub groupGenes {
  my ($self, $indexFile, $list, $org) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    print "Only Relation format 1.1 is supported\n";
    exit;
  }
  print "Reading Symbols ...\n";
  my $hashTable = &readIndex($indexFile);
  print "Done. Time=".times()."\n";
  print "Reading list ...\n";
  my $genes = {};
  my $fh;
  open($fh, "<$list") || die "Can't open $list\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split(/[^\w-]+/);
    foreach my $w (@list) {
      $w = uc($w);
      next if ($w =~ /^-*$/);
      next if ($w eq "UNKNOWN");
      if (defined $hashTable->[0]->{$w}) {
        foreach my $n (@{$hashTable->[0]->{$w}}) {
          $genes->{$n} = 1;
        }
      }
      if (defined $hashTable->[1]->{$w}) {
        foreach my $n (@{$hashTable->[1]->{$w}}) {
          $genes->{$n} = 1;
        }
      }
    }
  }
  close($fh);
  my $low = $self->getLow();
  my $high = $self->getHigh();
  my $balanced = $self->getBalanced();

  my $lowhash = {};
  my $highhash = {};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$low}); $i++) {
    $lowhash->{$low->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$high}); $i++) {
    $highhash->{$high->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }
  my $lowgenes = {};
  my $highgenes = {};
  my $badgenes = {};
  my $listHash = {};
  my ($numbal, $numlow, $numhigh, $numbad) = (0, 0, 0, 0);
  foreach my $k (keys(%{$genes})) {
    my $probeid = $hashTable->[2]->{$k};
    my $name = $hashTable->[3]->{$k};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    if (defined $lowhash->{$probeid}) {
      push @{$lowgenes->{$name}}, $k;
      $numlow++;
    }
    elsif (defined $highhash->{$probeid}) {
      push @{$highgenes->{$name}}, $k;
      $numhigh++
    }
    elsif (defined $balancedhash->{$probeid}) {
      push @{$listHash->{$name}}, $k;
      $numbal++;
    }
    else {
      push @{$badgenes->{$name}}, $k;
      $numbad++;
    }
  }
  print "Done. Time=".times()."\n";
  my $graphs = {};
  foreach my $name (keys(%{$listHash})) {
    foreach my $index (@{$listHash->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $relations = {};
      my $numrel = 0;
      my $i = $balancedhash->{$probeid};
#print "$num $numBits $numBytes $ptr ".$balancedhash->{$id}."\n";
      foreach my $name1 (sort keys(%{$listHash})) {
        foreach my $index1 (@{$listHash->{$name1}}) {
          my $probeid1 = $hashTable->[2]->{$index1};
          my $j = $balancedhash->{$probeid1};
          my $code = $self->readCode_1_1($i, $j);
          if ($code > 0) {
            if (!defined $relations->{$code}) {
              $relations->{$code} = [];
            }
            push @{$relations->{$code}}, $j;
#print "$code ".$balancedhash->{$id}." $i\n";
            $numrel++;
          }
        }
      }
      $graphs->{$i} = $relations;
    } # end $index
  } # end $name
  #my $groups = &divideGroups($graphs, $hashTable, $balanced, $graphs);
  #$groups = &divideGroups($graphs, $hashTable, $balanced, $groups->{0});
  my $graph = &getHiLo_LoLoGraph($graphs, $balanced);
  #&transitiveReduce($graph->[0], $graph->[1]);
  &printSGraph("test.sadj", $graph->[0], $graph->[1], $hashTable, $balanced, 1);
  #&printGMLGraph("test.gml", $graph->[0], $graph->[1], $hashTable, $balanced);
  #&printGraphStats($graph->[0], $graph->[1], $hashTable, $balanced, $graphs);
}

#
# $type = 0 -> raw numbers
# $type = 1 -> symbols
#
sub printSGraph {
  my ($file, $nodes, $edges, $hashTable, $balanced, $type) = @_;
  my $fh;
  open($fh, ">$file") || die "Can't open $file\n";
  foreach my $i (keys(%{$nodes})) {
    my $len = scalar(keys(%{$edges->{$i}}));
    if ($type == 1) {
      my $sym1 = &getSymbol($i, $hashTable, $balanced);
      print $fh $sym1,"\t", $len;
      foreach my $j (keys(%{$edges->{$i}})) {
        my $sym2 = &getSymbol($j, $hashTable, $balanced);
        print $fh "\t", $sym2;
        if ($sym1 eq "CD6_b-1566448_at" && $sym2 eq "LY9_b-210370_s_at") {
          print "$i $j\n";
        }
      }
    }
    else {
      print $fh "$i ", $len;
      foreach my $j (keys(%{$edges->{$n}})) {
        print $fh " ", $j;
      }
    }
    print $fh "\n";
  }
  close($fh);
}

sub printGraphStats {
  my ($nodes, $edges, $hashTable, $balanced, $graphs) = @_;
  print "Clusters :\n";
  foreach my $i (keys(%{$nodes})) {
    my $n = &getEquivalentNodes($graphs, $i);
    if (scalar(@$n) > 1) {
      my $name = &getSymbol($i, $hashTable, $balanced);
      print "$name : ";
      foreach my $j (@$n) {
        my $name = &getSymbol($j, $hashTable, $balanced);
        print "$name, ";
      }
      print "\n";
    }
  }

  my $links = {};
  foreach my $n (keys(%{$edges})) {
    next if ($n =~ /b/);
    foreach my $j (keys(%{$edges->{$n}})) {
      next if !($j =~ /b/);
      if (!defined $links->{$n}) {
        $links->{$n} = {};
      }
      $links->{$n}->{$j} = $edges->{$n}->{$j};
    }
  }
  foreach my $i (keys(%{$links})) {
    my $name = &getSymbol($i, $hashTable, $balanced);
    print "$name : \n";
    print "\tlinks:\t";
    foreach my $j (keys(%{$links->{$i}})) {
      my $name = &getSymbol($j, $hashTable, $balanced);
      print "$name, ";
    }
    print "\n\tlo->lo:\t";
    $i =~ s/b//g;
    foreach my $j (@{$graphs->{$i}->{2}}) {
      my $name = &getSymbol($j, $hashTable, $balanced);
      print "$name, ";
    }
    print "\n\tequv:\t";
    $i =~ s/b//g;
    foreach my $j (@{$graphs->{$i}->{5}}) {
      my $name = &getSymbol($j, $hashTable, $balanced);
      print "$name, ";
    }
    print "\n\tSucc-Lo:\t";
    &printSuccessorLo($nodes, $edges, $i."b", $hashTable, $balanced);
    print "\n-------------------------\n";
  }
}

sub printSuccessorLo {
  my ($nodes, $edges, $start, $hashTable, $balanced) = @_;
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
    foreach my $j (keys(%{$edges->{$i}})) {
        if ($visited->{$j} != 1) {
          $visited->{$j} = 1;
          if ($j =~ /b/) {
            $parent->{$j} = $i;
            push @queue, $j;
            #print $j, " ";
            my $name = &getSymbol($j, $hashTable, $balanced);
            print "$name, ";
          }
        }
    }
    #print "\n";
  }
}

sub transitiveReduce {
  my ($nodes, $edges) = @_;
  &transitiveReduceHiLo($nodes, $edges);
  &transitiveReduceLoLo($nodes, $edges);
  &transitiveReduceHiHi($nodes, $edges);
}

sub transitiveReduceHiLo {
  my ($nodes, $edges) = @_;
  my @links;
  foreach my $n (keys(%{$edges})) {
    next if ($n =~ /b/);
    foreach my $j (keys(%{$edges->{$n}})) {
      next if !($j =~ /b/);
      push @links, [$n, $j,  $edges->{$n}->{$j}];
    }
  }
  my $size = scalar(@links);
  my $index = 0;
  my $lastPerc = 0;
  my $num = 0;
  print "Size(HiLo) before: ", $size, "\n";
  foreach (@links) {
    my ($n, $j, $w) = @{$_};
    my $perc = int($index * 100/$size);
    if ( ($perc -  $lastPerc) >= 1 ) {
        $lastPerc = $perc;
        my $time = times();
        print STDERR "Percent finished : $perc Time = $time\n";
    } 
    my $path = &findPath($nodes, $edges, $n, $j, 2, 2000);
    if (scalar(@{$path}) > 2) {
      delete $edges->{$n}->{$j};
      $num ++;
    }
    $index ++;
  }
  print "Size after: ", $size -$num, "\n";
}

sub transitiveReduceLoLo {
  my ($nodes, $edges) = @_;
  my @links;
  foreach my $n (keys(%{$edges})) {
    next if !($n =~ /b/);
    foreach my $j (keys(%{$edges->{$n}})) {
      next if !($j =~ /b/);
      push @links, [$n, $j,  $edges->{$n}->{$j}];
    }
  }
  my $size = scalar(@links);
  my $index = 0;
  my $lastPerc = 0;
  my $num = 0;
  print "Size(LoLo) before: ", $size, "\n";
  foreach (@links) {
    my ($n, $j, $w) = @{$_};
    my $perc = int($index * 100/$size);
    if ( ($perc -  $lastPerc) >= 1 ) {
        $lastPerc = $perc;
        my $time = times();
        print STDERR "Percent finished : $perc Time = $time\n";
    } 
    my $path = &findPath($nodes, $edges, $n, $j, 2, 2000);
    if (scalar(@{$path}) > 2) {
      delete $edges->{$n}->{$j};
      $num ++;
    }
    $index ++;
  }
  print "Size after: ", $size -$num, "\n";
}

sub transitiveReduceHiHi {
  my ($nodes, $edges) = @_;
  my @links;
  foreach my $n (keys(%{$edges})) {
    next if ($n =~ /b/);
    foreach my $j (keys(%{$edges->{$n}})) {
      next if ($j =~ /b/);
      push @links, [$n, $j,  $edges->{$n}->{$j}];
    }
  }
  my $size = scalar(@links);
  my $index = 0;
  my $lastPerc = 0;
  my $num = 0;
  print "Size(HiHi) before: ", $size, "\n";
  foreach (@links) {
    my ($n, $j, $w) = @{$_};
    my $perc = int($index * 100/$size);
    if ( ($perc -  $lastPerc) >= 1 ) {
        $lastPerc = $perc;
        my $time = times();
        print STDERR "Percent finished : $perc Time = $time\n";
    } 
    my $path = &findPath($nodes, $edges, $n, $j, 2, 2000);
    if (scalar(@{$path}) > 2) {
      delete $edges->{$n}->{$j};
      $num ++;
    }
    $index ++;
  }
  print "Size after: ", $size -$num, "\n";
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

sub printGMLGraph {
  my ($file, $nodes, $edges, $hashTable, $balanced, $props) = @_;

  my $fh;
  open($fh, ">$file") || die "Can't open $file\n";
  print $fh "graph [\n";
  print $fh "directed 1\n";
  print $fh "id 1\n";
  foreach my $n (keys(%{$nodes})) {
    my $l = &getSymbol($n, $hashTable, $balanced);
    print $fh "node [\n";
    print $fh "\tid ".$nodes->{$n}."\n\tlabel \"$l\"\n";
    my $nc = "#3cbbc6";
    if (defined $props) {
      my $name = &getSymbol($n, $hashTable, $balanced);
      if (defined $props->{$name}) {
        $nc = "#bb0000";
      }
    }
    print $fh <<END;
\tgraphics [
\t\tw       30.0
\t\th       30.0
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
      print $fh "edge [\n";
      print $fh "\tsource ".$nodes->{$n}."\n\ttarget ".$nodes->{$j}."\n";
      print $fh <<END;
\tgraphics [
\t\twidth   2
\t\ttype    "line"
\t\tfill    "#3b3b3b"
\t\tarrow   "last"
\t]
END
      print $fh "]\n";
    }
  }
  print $fh "]\n";
  close($fh);
}

sub getHiLoGraph {
  my $graphs = shift;
  my $nodes = {};
  my $edges = {};
  my $index = 0;
  foreach my $i (keys(%{$graphs})) {
    my $n = &getSuperNode($graphs, $i);
    if (scalar(@{$graphs->{$n}->{4}}) > 0) {
      if (scalar(@{$graphs->{$n}->{3}}) == 0) {
        if (!defined $nodes->{$n}) {
            $nodes->{$n} = $index;
            $index++;
        }
        foreach $j (@{$graphs->{$n}->{4}}) {
          my $m = &getSuperNode($graphs, $j);
          if (!defined $nodes->{$m}) {
            $nodes->{$m} = $index;
            $index++;
          }
          if (!defined $edges->{$n}) {
            $edges->{$n} = {};
          } 
          $edges->{$n}->{$j} = 1;
          if (!defined $edges->{$m}) {
            $edges->{$m} = {};
          } 
          $edges->{$m}->{$n} = 1;
        }
      }
    }
  }
  return [$nodes, $edges];
}

sub getEquivalentNodes {
  my ($graphs, $i) = @_;
  my @equiv = @{$graphs->{$i}->{5}};
  my $hash = {};
  $hash->{$i} = 1;
  foreach (@equiv) {
    $hash->{$_} = 1;
    foreach (@{$graphs->{$_}->{5}}) {
      $hash->{$_} = 1;
    }
  }
  my @list = sort { 
        scalar(@{$graphs->{$b}->{5}}) <=> scalar(@{$graphs->{$a}->{5}}) ||
        $a <=> $b } keys(%{$hash});
  return [@list];
}

sub getSuperNode {
  my ($graphs, $i) = @_;
  return $i;
  my @equiv = @{$graphs->{$i}->{5}};
  my $hash = {};
  $hash->{$i} = 1;
  foreach (@equiv) {
    $hash->{$_} = 1;
    foreach (@{$graphs->{$_}->{5}}) {
      $hash->{$_} = 1;
    }
  }
  my @list = sort { 
        scalar(@{$graphs->{$b}->{5}}) <=> scalar(@{$graphs->{$a}->{5}}) ||
        $a <=> $b } keys(%{$hash});
  return $list[0];
}

sub getHiLo_LoLoGraph {
  my $graphs = shift;
  my $balanced = shift;
  my $nodes = {};
  my $edges = {};
  my $index = 0;
  foreach my $i (keys(%{$graphs})) {
    my $n = &getSuperNode($graphs, $i);
    if (!defined $nodes->{$n}) {
      $nodes->{$n} = $index;
      $index++;
    }
    if (!defined $nodes->{$n."b"}) {
      $nodes->{$n."b"} = $index;
      $index++;
    }
    foreach $j (@{$graphs->{$n}->{4}}) {
      my $m = &getSuperNode($graphs, $j);
      if (!defined $nodes->{$m}) {
        $nodes->{$m} = $index;
        $index++;
      }
      if (!defined $nodes->{$m."b"}) {
        $nodes->{$m."b"} = $index;
        $index++;
      }
      if (!defined $edges->{$n}) {
        $edges->{$n} = {};
      } 
      $edges->{$n}->{$m."b"} = 1;
      if (!defined $edges->{$m}) {
        $edges->{$m} = {};
      } 
      $edges->{$m}->{$n."b"} = 1;
    }
    foreach $j (@{$graphs->{$n}->{2}}) {
      my $m = &getSuperNode($graphs, $j);
      #if ($balanced->[$n] eq "1566448_at" && $balanced->[$m] eq "210370_s_at") {
      #  print "CD6 1566448_at entered 2 : $i ",  $balanced->[$i],"\n";
      #  print "LY9 210370_s_at entered 2 : $j ",  $balanced->[$j],"\n";
      #}
      if (!defined $nodes->{$m}) {
        $nodes->{$m} = $index;
        $index++;
      }
      if (!defined $nodes->{$m."b"}) {
        $nodes->{$m."b"} = $index;
        $index++;
      }
      if (!defined $edges->{$n."b"}) {
        $edges->{$n."b"} = {};
      } 
      $edges->{$n."b"}->{$m."b"} = 1;
      if (!defined $edges->{$m}) {
        $edges->{$m} = {};
      } 
      $edges->{$m}->{$n} = 1;
      #if ($n == 3342 && $m == 12083) {
      #  print "entered 2 : $n $i ",  $balanced->[$i],"\n";
      #  print "entered 2 : $m $j ",  $balanced->[$j],"\n";
      #}
    }
    foreach $j (@{$graphs->{$n}->{3}}) {
      my $m = &getSuperNode($graphs, $j);
      #if ($balanced->[$n] eq "1566448_at" && $balanced->[$m] eq "210370_s_at") {
      #  print "CD6 1566448_at entered 3 : $i ",  $balanced->[$i],"\n";
      #  print "LY9 210370_s_at entered 3 : $j ",  $balanced->[$j],"\n";
      #}
      if (!defined $nodes->{$m}) {
        $nodes->{$m} = $index;
        $index++;
      }
      if (!defined $nodes->{$m."b"}) {
        $nodes->{$m."b"} = $index;
        $index++;
      }
      if (!defined $edges->{$n}) {
        $edges->{$n} = {};
      } 
      $edges->{$n}->{$m} = 1;
      if (!defined $edges->{$m."b"}) {
        $edges->{$m."b"} = {};
      } 
      $edges->{$m."b"}->{$n."b"} = 1;
      #if ($m == 3342 && $n == 12083) {
      #  print "entered 3 : $n $i ",  $balanced->[$i],"\n";
      #  print "entered 3 : $m $j ",  $balanced->[$j],"\n";
      #}
    }
  } #end $i
  return [$nodes, $edges];
}

sub divideGroups {
  my ($graphs, $hashTable, $balanced, $list) = @_;
  my $groups = {};
  $groups->{0} = {};
  $groups->{1} = {};
  foreach my $i (keys(%{$list})) {
    #my $probeid = $balanced->[$i];
    #my $index = $hashTable->[0]->{uc($probeid)}->[0];
    #my $name = $hashTable->[3]->{$index};
    if (scalar(@{$graphs->{$i}->{4}}) > 0) {
      #print "$name($probeid) :\n";
      my $c = &findGroup($groups, $graphs, $i);
      if ($c == 0) {
        #print "Found in 0 \n";
        foreach my $j (@{$graphs->{$i}->{4}}) {
          #my $probeid = $balanced->[$j];
          #my $index = $hashTable->[0]->{uc($probeid)}->[0];
          #my $name = $hashTable->[3]->{$index};
          #print "$name,";
          if (!defined $groups->{0}->{$j} && defined $list->{$j}) {
            $groups->{1}->{$j} = 1;
          }
        }
        #print "\n";
      }
      elsif ($c == 1) {
        #print "Found in 1 \n";
        foreach my $j (@{$graphs->{$i}->{4}}) {
          #my $probeid = $balanced->[$j];
          #my $index = $hashTable->[0]->{uc($probeid)}->[0];
          #my $name = $hashTable->[3]->{$index};
          #print "$name,";
          if (!defined $groups->{1}->{$j} && defined $list->{$j}) {
            $groups->{0}->{$j} = 1;
          }
        }
        #print "\n";
      }
      else {
        #print "Not found\n";
        $groups->{0}->{$i} = 1;
        foreach my $j (@{$graphs->{$i}->{4}}) {
          #my $probeid = $balanced->[$j];
          #my $index = $hashTable->[0]->{uc($probeid)}->[0];
          #my $name = $hashTable->[3]->{$index};
          #print "$name,";
          if (!defined $groups->{0}->{$j} && defined $list->{$j}) {
            $groups->{1}->{$j} = 1;
          }
        }
        #print "\n";
      }
    }
  }
  foreach my $i (keys(%{$groups})) {
    print "$i : ";
    foreach my $j (keys(%{$groups->{$i}})) {
      my $probeid = $balanced->[$j];
      my $index = $hashTable->[0]->{uc($probeid)}->[0];
      my $name = $hashTable->[3]->{$index};
      print "$name($probeid),";
      #print "$probeid,";
    }
    print "\n";
  }
  return $groups;
}

sub findGroup {
  my ($groups, $graphs, $i) = @_;
  if (defined $groups->{0}->{$i}) {
    return 0;
  }
  if (defined $groups->{1}->{$i}) {
    return 1;
  }
  foreach my $j (@{$graphs->{$i}->{4}}) {
    if (defined $groups->{0}->{$j}) {
      return 1;
    }
    if (defined $groups->{1}->{$j}) {
      return 0;
    }
  } 
  return -1;
}

sub getCodes {
  my ($n, $id) = @_;
  my $balanced = $n->getBalanced();
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }
  my $num = $n->getNum();
  my $i = $balancedhash->{$id};
  my $relations = {};
  if (!defined $i) {
    return $relations;
  }
  for (my $j = 0; $j < $num/2; $j++) {
    my $code = $n->readCode_1_1($i, $j);
    if ($code > 0) {
      if (!defined $relations->{$code}) {
        $relations->{$code} = [];
      }
      push @{$relations->{$code}}, $balanced->[$j];
    }
  }
  return $relations;
}

sub getRelations {
  my($self, $indexFile, $probeid) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version ne "1.1") {
    print STDERR "$version: Only version 1.1 supported\n";
    exit(1);
  }
   
  $probeid = uc($probeid);
  my $file = $self->{'file'};
  my $hashTable = &readIndex($indexFile);
  if (!defined $hashTable->[0]->{$probeid}) {
    print STDERR "Cannot fine probe $probeid\n";
    exit(1);
  }
  my $low = $self->getLow();
  my $high = $self->getHigh();
  my $balanced = $self->getBalanced();

  my $lowhash = {};
  my $highhash = {};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$low}); $i++) {
    $lowhash->{$low->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$high}); $i++) {
    $highhash->{$high->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }
  my $index = $hashTable->[0]->{$probeid}->[0];
  my $name = $hashTable->[3]->{$index};
  my $id = $hashTable->[2]->{$index};
  my $cname = &getSplitName($name, 20);
  my $genename = $name;
  $genename =~ s/\/.*//g;
  if (defined $lowhash->{$id}) {
    return [0];
  }
  if (defined $highhash->{$id}) {
    return [1];
  }
  if (!defined $balancedhash->{$id}) {
    return [3];
  }
  my $ptr = $self->getStartPtr();
  my $num = $self->getNum();
  my $numBits = $self->getNumBits();
  my $numBytes = $self->getNumBytes();
  my $relations = {};
  my $numrel = 0;
  my $i = $balancedhash->{$id};
  #print "$num $numBits $numBytes $ptr ".$balancedhash->{$id}."\n";
  for (my $j = 0; $j < $num/2; $j++) {
    my $code = $self->readCode_1_1($i, $j);
    if ($code > 0) {
      if (!defined $relations->{$code}) {
        $relations->{$code} = [];
      }
      push @{$relations->{$code}}, $balanced->[$j];
      $numrel++;
    }
  }
  return [2, $relations];
}

sub printText {
  my($self, $indexFile, $probeid, $type) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version ne "1.1") {
    print "Only version 1.1 supported\n";
    exit(1);
  }
   
  $probeid = uc($probeid);
  my $file = $self->{'file'};
  my $hashTable = &readIndex($indexFile);
  if (!defined $hashTable->[0]->{$probeid}) {
    print "Cannot fine probe $probeid\n";
    return;
  }
  my $low = $self->getLow();
  my $high = $self->getHigh();
  my $balanced = $self->getBalanced();

  my $lowhash = {};
  my $highhash = {};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$low}); $i++) {
    $lowhash->{$low->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$high}); $i++) {
    $highhash->{$high->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }
  my $index = $hashTable->[0]->{$probeid}->[0];
  my $name = $hashTable->[3]->{$index};
  my $id = $hashTable->[2]->{$index};
  my $cname = &getSplitName($name, 20);
  my $genename = $name;
  $genename =~ s/\/.*//g;
  if (defined $lowhash->{$id}) {
    print "$id\t$genename\t(Always Low)\n";
    return;
  }
  if (defined $highhash->{$id}) {
    print "$id\t$genename\t(Always High)\n";
    return;
  }
  if (!defined $balancedhash->{$id}) {
    print "$id\t$genename\t(Bad dynamic range)\n";
    return;
  }
  my $ptr = $self->getStartPtr();
  my $num = $self->getNum();
  my $numBits = $self->getNumBits();
  my $numBytes = $self->getNumBytes();
  my $relations = {};
  my $numrel = 0;
  my $i = $balancedhash->{$id};
  #print "$num $numBits $numBytes $ptr ".$balancedhash->{$id}."\n";
  for (my $j = 0; $j < $num/2; $j++) {
    my $code = $self->readCode_1_1($i, $j);
    if ($code > 0) {
      if (!defined $relations->{$code}) {
        $relations->{$code} = [];
      }
      push @{$relations->{$code}}, $j;
      #print "$code ".$balancedhash->{$id}." $i\n";
      $numrel++;
    }
  }
  my @status = ("No relation", "($genename low -> B high)",
  "($genename low -> B low)", "($genename high -> B high)",
  "($genename high -> B low)", "Equivalent", "Opposite");
  my $list = $relations->{$type};  
  print "$id\t$genename\t$numrel\t", scalar(@{$list}), "\t", $status[$type], "\n";
  for (my $i = 0; $i < scalar(@{$list}); $i++) {
    my $probeid = $balanced->[$list->[$i]];
    my $nid = $hashTable->[0]->{uc($probeid)}->[0];
    my $name = $hashTable->[3]->{$nid};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    print "$probeid\t$name\n";
  }
}

sub getIDHash {
    my $file = shift;
    return undef if (!defined $file);
    my $fh;
    my $res = {};
    open($fh, "<$file") || die "Can't open $file\n";
    while (<$fh>) {
        s/[\r\n]//g;
        my @list = split("\t");
        my $name = $list[0];
        $name =~ s/\s//g;
        if (!($name =~ /^\s*$/)) {
          #print $name, "\n";
          $res->{$name} = 1;
        }
    }
    close($fh);
    return $res;
}

sub unigeneToAffy {
    my ($annotation, $hash) = @_;
    return undef if (!defined $annotation);
    my $fh;
    my $res = [];
    open($fh, "<$annotation") || die "Can't open $annotation\n";
    while (<$fh>) {
        s/[\r\n]//g;
        my @list = split("\",");
        my $id = $list[0];
        $id =~ s/[\s\"]//g;
        my $name1 = $list[9];
        $name1 =~ s/[\s\"]//g;
        my $name2 = $list[10];
        $name2 =~ s/[\s\"]//g;
        if (defined $hash->{$name1} || defined $hash->{$name2}) {
            #print "$name1\_$name2\n";
            push @{$res}, $id;
        }
    }
    close($fh);
    return $res;
}

sub getNameHashFromIdx {
  my $indexFile = shift;
  return undef if (!defined $indexFile);
  my $namehash = {} ;
  open(FL, "<$indexFile") || die "Can't open $indexFile\n";
  <FL>; # Header
  my $index = 0;
  while (<FL>) {
    my ($id, $name, $ptr, $desc) = split("\t", $_);
    $namehash->{$id} = $name;
  }
  close(FL);
  return $namehash;
}

sub getNameHashFromIdxHegemon {
  my $indexFile = shift;
  return undef if (!defined $indexFile);
  my $namehash = {} ;
  open(FL, "<$indexFile") || die "Can't open $indexFile\n";
  <FL>; # Header
  my $index = 0;
  while (<FL>) {
    my ($id, $ptr, $name, $desc) = split("\t", $_);
    $namehash->{$id} = $name;
  }
  close(FL);
  return $namehash;
}

sub getNameHashFromAnnotation {
    my $file = shift;
    return undef if (!defined $file);
    my $fh;
    my $res = {};
    open($fh, "<$file") || die "Can't open $file\n";
    while (<$fh>) {
        s/[\r\n]//g;
        my @list = split(/,"/);
        $list[0] =~ s/[\s\"]//g;
        $name = $list[14].": ".$list[13];
        $name =~ s/":/:/g;
        $name =~ s/"$//g;
        $name =~ s/:.*$//g;
        $res->{$list[0]} = $name;
    }
    close($fh);
    return $res;
}

# Mouse430_2.na28.ortholog.csv, HG-U133_Plus_2
sub getOrthologHash {
    my ($file, $org) = @_;
    return undef if (!defined $file);
    my $loop = 0;
    if ($file =~ /$org/) {
        $loop = 1;
    }
    my $fh;
    my $res = {};
    open($fh, "<$file") || die "Can't open $file\n";
    while (<$fh>) {
        s/[\r\n]//g;
        my @list = split(/,"/);
        $list[0] =~ s/[\s\"]//g;
        $list[2] =~ s/[\s\"]//g;
        $list[3] =~ s/[\s\"]//g;
        if ($loop == 1) {
            if (!defined $res->{$list[0]}) {
                $res->{$list[0]} = [];
                push @{$res->{$list[0]}}, lc($list[0]);
            }
        }
        elsif ($list[3] eq $org) {
            if (!defined $res->{$list[0]}) {
                $res->{$list[0]} = [];
            }
            push @{$res->{$list[0]}}, lc($list[2]);
        }
    }
    close($fh);
    return $res;
}

sub getPhenotype {
  my ($pcl, $phfile, $id) = @_;
  my $fh;
  open($fh, "<$pcl") || die "Can't open $pcl\n";
  my $head = <$fh>;
  chop $head;
  close($fh);
  my @header = split("\t", $head);
  if ($id eq "All") {
    my $res = [map { 1 } @header];
    $res->[0] = $id;
    return $res;
  }
  open($fh, "<$phfile") || die "Can't open $phfile\n";
  my $res = [map { 0 } @header];
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    if ($list[0] eq $id) {
      $res = [@list];
      last;
    }
  }
  return $res;
}

sub getDynamicRangeThreshold {
  my ($file, $file2, $type, $value) = @_;
  my $fh;
  my $fh2;
  open($fh, "<$file") || die "Can't open $file\n";
  open($fh2, "<$file2") || die "Can't open $file2\n";
  my $head = <$fh>;
  $head =~ s/[\r\n]//g;
  my @header = split("\t", $head);
  my $hash = {};
  for (my $i = 1; $i < scalar(@header); $i++) {
    $hash->{$header[$i]} = $i;
  }
  my $head2 = <$fh2>;
  $head2 =~ s/[\r\n]//g;
  my @header2 = split("\t", $head2);
  my $hash2 = {};
  for (my $i = 1; $i < scalar(@header2); $i++) {
    $hash2->{$header2[$i]} = $i;
  }
  my $index;
  if ($type eq "Perc") {
    $index = $hash->{sprintf("$type %.2f", $value)};
    if (!defined $index) {
      die "Can't access $type $value in the dynamic range\n";
    }
  }
  my $res = {};
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split("\t");
    $_ = <$fh2>;
    s/[\r\n]//g;
    my @list2 = split("\t");
    my $min = $list[1];
    my $max = $list[$hash->{"Perc 1.00"}];
    my $thr = $list[$hash->{"Perc 0.00"}];
    my ($v, $p, $a, $n);
    $n = $list2[$hash2->{"Expr 0.00"}];
    if ($type eq "Perc") {
      $v = $list[$index];
      $p = sprintf("%.2f", $value);
      if ($v > $thr) {
        $a = sprintf("%.2f", ($v - $thr) / ($max - $thr));
      }
      else {
        $a = sprintf("%.2f", ($v - $thr) / ($thr - $min));
      }
    }
    else {
      $a = sprintf("%.2f", $value);
      $p = $list[$hash->{"Expr $a"}];
      if ($value > 0) {
        $v = $thr + $value * ($max - $thr);
      }
      else {
        $v = $thr + $value * ($thr - $min);
      }
    }
    $res->{$list[0]} = [$v, $p, $a, $n];
  }
  close($fh);
  return $res;
}

1;
__END__

=head1 NAME

Network - A read/write interface 

Example:

  use Network;

  my $n = Network->new(file => "network.rl", mode => "r");

  $n->readNetworkFile(); # read headers

  my $low = $n->getLow(); # get the list of (Affy IDs) low genes

  my $high = $n->getHigh(); # get the list of (Affy IDs) high genes

  # get the list of genes for pairwise data
  my $balanced = $n->getBalanced();

  my $numBits = $n->getNumBits(); # get the number of bits used

  # get the network as an array of bit vectors
  # each row corresponds to a gene (same index) in the balanced list
  # An example bit vector :
  #    ith row - > 000 001 010 100 000 110
  # The relations are : (i, 1, code=4), (i, 2, code=2),
  #                     (i, 3, code=1), (i, 5, code=3)
  #  Note that the lsb and msb is reversed in the bit encoding.
  my $network = $n->readAll();

  # open a new Network for writing
  my $w = Network->new(file => "test.rl", mode => "w");

  # write using the specified lists
  $w->write($low, $high, $balanced, $network, $numBits);

  # close writing 
  $w->close();

=cut
