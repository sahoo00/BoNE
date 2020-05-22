package HR;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my %params = @_;
  my $self  = {};

  foreach (keys(%params)) {
    $self->{$_} = $params{$_};
  }

  bless ($self, $class);
  return $self;
}

sub run {
  my ($self, $str) = @_;
  print $str, "\n";
}

sub get {
  my ($self, $str) = @_;
  return "";
}

sub stop {
  my ($self, $str) = @_;
}

1;
