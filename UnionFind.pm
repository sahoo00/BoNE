# Union-find data structure
package UnionFind;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $self  = {};
  $self->{'parent'}  = {};
  $self->{'rank'}  = {};
  bless ($self, $class);
  return $self;
}

sub makeset {
  my ($self, $a) = @_;
  if (!defined $self->{'parent'}->{$a}) {
    $self->{'parent'}->{$a} = $a;
    $self->{'rank'}->{$a} = 0;
  }
}

sub union {
  my ($self, $a, $b) = @_;
  $self->linkset($self->findset($a), $self->findset($b));
}

sub linkset {
  my ($self, $a, $b) = @_;
  return if (!defined $a);
  return if (!defined $b);
  if ($self->{'rank'}->{$a} > $self->{'rank'}->{$b}) {
    $self->{'parent'}->{$b} = $a;
  }
  else {
    $self->{'parent'}->{$a} = $b;
    if ($self->{'rank'}->{$a} == $self->{'rank'}->{$b}) {
      $self->{'rank'}->{$b}++;
    }
  }
}

sub findset {
  my ($self, $a) = @_;
  if (!defined $self->{'parent'}->{$a}) {
    return undef;
  }
  if ($a ne $self->{'parent'}->{$a}) {
    $self->{'parent'}->{$a} = $self->findset($self->{'parent'}->{$a});
  }
  return $self->{'parent'}->{$a};
}

sub delete {
  my ($self, $a) = @_;
  if (!defined $self->findset($a)) {
    return undef;
  }
  my $clusters = {};
  foreach my $k (keys %{$self->{'parent'}}) {
    next if ($k eq $a);
    push @{$clusters->{$self->findset($k)}}, $k;
  }
  my $set = UnionFind->new();
  foreach my $k (keys %{$clusters}) {
    my $key = $k;
    if ($k eq $a) {
      $key = $clusters->{$k}->[0];
    }
    $set->makeset($key);
    foreach my $i (@{$clusters->{$k}}) {
      $set->makeset($i);
      $set->union($key, $i);
    }
  }
  return $set;
}

1;

__END__

=head1 NAME

UnionFind - Provides Union-Find data structure

Example:

  use UnionFind;

  my $n = UnionFind->new();
  $n->makeset($a);
  $n->union($a, $b);
  $n->findset($a);

=cut
