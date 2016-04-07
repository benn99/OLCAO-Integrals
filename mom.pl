#!use/bin/perl

##################################################################
#
# Author: Ben Walker
#
# This script generates the formatted 20x20 matrix elements for
# the momentum matrix integrals.  The overlap formalism is readily 
# extendable to orbitals higher than f-type (g-type orbitals, etc.).  
# It makes use of the "triad" formalism established by Browne & 
# Poshusta in 1962 the reference for wihch is the following:
#
# ****************************************************************
# Browne JC, Poshusta RD., "Quantum-Mechanical Integrals 
# over Gaussian Atomic Orbitals.", The Journal of Chemical Physics, 
# Vol. 36, Iss. 7, P. 19337, (1962).
# ****************************************************************
#
# The output of this program is the integral solutions in
# the form of Fortran 90 code, which is to be used in OLCAO.
##################################################################

##################################################################
# Use necessary modules
##################################################################
use strict;
use warnings;
use Env;
use POSIX;

##################################################################
# Begin main program
##################################################################

my @eTriad;
&print_Mom;
##################################################################
# End main program
##################################################################

##################################################################
# Define subroutines
##################################################################

# This subroutine initializes the triads that are to represent
# the types of orbitals.
sub triads {
  # Initialize the triads
  @{$eTriad[1]}  = qw(0 0 0); # s
  @{$eTriad[2]}  = qw(1 0 0); # x
  @{$eTriad[3]}  = qw(0 1 0); # y
  @{$eTriad[4]}  = qw(0 0 1); # z
  @{$eTriad[5]}  = qw(2 0 0); # xx
  @{$eTriad[6]}  = qw(0 2 0); # yy
  @{$eTriad[7]}  = qw(0 0 2); # zz
  @{$eTriad[8]}  = qw(1 1 0); # xy
  @{$eTriad[9]}  = qw(1 0 1); # xz
  @{$eTriad[10]} = qw(0 1 1); # yz
  @{$eTriad[11]} = qw(1 1 1); # xyz
  @{$eTriad[12]} = qw(2 1 0); # xxy
  @{$eTriad[13]} = qw(2 0 1); # xxz
  @{$eTriad[14]} = qw(0 2 1); # yyz
  @{$eTriad[15]} = qw(1 2 0); # yyx
  @{$eTriad[16]} = qw(1 0 2); # zzx
  @{$eTriad[17]} = qw(0 1 2); # zzy
  @{$eTriad[18]} = qw(3 0 0); # xxx
  @{$eTriad[19]} = qw(0 3 0); # yyy
  @{$eTriad[20]} = qw(0 0 3); # zzz
  #@{$eTriad[21]} = qw(0 1 3);
  #@{$eTriad[22]} = qw(0 3 1);
  #@{$eTriad[23]} = qw(1 0 3);
  #@{$eTriad[24]} = qw(1 3 0);
  #@{$eTriad[25]} = qw(3 0 1);
  #@{$eTriad[26]} = qw(3 1 0);
  #@{$eTriad[27]} = qw(4 0 0);
  #@{$eTriad[28]} = qw(0 4 0);
  #@{$eTriad[29]} = qw(0 0 4);
  #@{$eTriad[30]} = qw(1 1 2);
  #@{$eTriad[31]} = qw(1 2 1);
  #@{$eTriad[32]} = qw(2 1 1);
}

# This subroutine calculates binomial coefficients
sub choose {
  # Define passed parameters
  my $n = $_[0];
  my $k = $_[1];    # These are for "n choose k"

  # Define local variables
  my $l;
  my $returnValue;
  $returnValue = 1;
  foreach $l (1..$k) {
    $returnValue *= ($n - $l + 1)/$l;
  }
  return $returnValue;
}

# This subroutine computes the factorial (x!)
# of a positive integer passed to it.
sub factorial {
  # Define passed parameters
  my $inputNum = $_[0];

  # Define local variables
  my $i;
  my $returnValue;

  $returnValue = 1;

  if ($inputNum == 0 or $inputNum == 1) {
    return $returnValue;  
  }
  else {
    for ($i=1;$i<=$inputNum;$i++) {
      $returnValue *= $i;
    }
  }
  return $returnValue;
}

# This subroutine calculates the double factorial.
sub dFact {
  # Define passed parameters
  my $inputNumber = $_[0];
  # Define local variables
  my $i;
  my $product;
  my $returnValue;
  $returnValue = 1;
  if ($inputNumber == 0 or $inputNumber == 1)
    {return $returnValue;}
  else {
    for ($i = 1;((2*$i) - 1) <= ($inputNumber);$i++)
      {$returnValue *= (2*$i) - 1;}
  }
  return $returnValue;
}

# This subroutine computes a single x-, y-, or 
# z-component of the 2-center overlap integral
# according to the Cook formalism:
#
# ****************************************************************
# Cook, David B., "Handbook of Computational Quantum Chemistry."
# Dover Publications (August 22, 2005).
# ISBN-10: 0486443078
# ISBN-13: 978-0486443072
# ****************************************************************
#
# The output of this program is the integral solutions in
# the form of Fortran 90 code, which is to be used in OLCAO.
sub twoCenterSum {
  # Define passed parameters.
  my $l1 = $_[0];
  my $l2 = $_[1];
  my $dim = $_[2];
  $dim += 1;
  # Define local variables.
  my $j;
  my $string;
  my $l;
  my $s;
  $l = $l1 + $l2;
  $string = "";
  $string .= "("; # Outer paren
  foreach $j (0..floor($l/2)) {
    foreach $s (0..(2*$j)) {
      my $power1 = $l1 - $s;
      my $power2 = $l2 - 2*$j + $s;
      if (&choose($l1,$s)*&choose($l2,(2*$j - $s)) == 0) {
        $string .= "0";
      }
      else {
        $string .= &choose($l1,$s)*&choose($l2,(2*$j - $s))*&dFact(2*$j - 1)/(2**$j);
        $string .= "*";
        $string .= "(PA($dim)**$power1)*(PB($dim)**$power2)*(Y**(-$j))";
      }
      if ($j != floor($l/2) or $s != (2*$j)) {
        $string .= " + ";
      }
      else {
        $string .= "";
      }
    }
  }
  $string .= ")"; # Outer paren
  return $string;
}

# This subroutine uses the 2-center OL functionality
# to construct the X momentum matrix formulas.  These integrals 
# are weighted sums of 2-center OL integrals.  This
# uses Cook's formalism (referenced above).
sub Mom_x {
  # Define passed parameters
  my $orb1 = $_[0];
  my $orb2 = $_[1];
  # Define local variables
  my @eTriad1;
  my @eTriad2;
  my $orbital1;
  my $orbital2;
  my $l1;
  my $l2;
  my $string = "";
  my $dim;
  #my $KE_CodeFile = "KE_Intg";
  my $num;

  # Open the output file
  #open (MOX,">$KE_CodeFile") || die "Error opening file.\n";

  &triads; # Initialize triads

  foreach $orbital2 (1..$orb2) {
    foreach $orbital1 (1..$orb1) {
      @eTriad1[0..2] = @{$eTriad[$orbital1]};
      @eTriad2[0..2] = @{$eTriad[$orbital2]};
      #$string .= "Mom_x($orbital1,$orbital2) = ";
      $string .= "wo(1,$orbital1,$orbital2) = ";
      $string .= "preFactor*("; # Outer paren
      $string .= $eTriad2[0];
      $string .= "*";
      $eTriad2[0] -= 1;
      if ($eTriad2[0] >= 0) {
        $string .= "(";
        foreach $dim (0..2) {
          $l1 = $eTriad1[$dim];
          $l2 = $eTriad2[$dim];
          $string .= &twoCenterSum($l1,$l2,$dim);
          if ($dim != 2) {
            $string .= "*";
          }
        }
        $string .= ")";
        $string .= " - ";
      }
      else {
        $string .= "(0.d0) - ";
      }
      $eTriad2[0] += 2;
      $string .= "2*a2*(";
      foreach $dim (0..2) {
        $l1 = $eTriad1[$dim];
        $l2 = $eTriad2[$dim];
        $string .= &twoCenterSum($l1,$l2,$dim);
        if ($dim != 2) {
          $string .= "*";
        }
      }
      $string .= ")";      
      $string .= ")"; # Outer paren

      #**************************************************************
      # Start RegEXP portion here
      #**************************************************************
      # Eliminate terms raised to 0-power
      if ($string =~ m/\(PA\(1\)\*\*0\)/) {
        $string =~ s/\(PA\(1\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PA\(2\)\*\*0\)/) {
        $string =~ s/\(PA\(2\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PA\(3\)\*\*0\)/) {
        $string =~ s/\(PA\(3\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(1\)\*\*0\)/) {
        $string =~ s/\(PB\(1\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(2\)\*\*0\)/) {
        $string =~ s/\(PB\(2\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(3\)\*\*0\)/) {
        $string =~ s/\(PB\(3\)\*\*0\)/1/g;
      }
      if ($string =~m/\(1\/\(\(2\*Y\)\*\*0\)\)/) {
        $string =~ s/\(1\/\(\(2\*Y\)\*\*0\)\)/1/g;
      }

      # Get rid of some damn zeroes
      if ($string =~ m/ \+ 0/) {
        $string =~ s/ \+ 0//g;
      }
      if ($string =~ m/0\*\(0.d0\) - /) {
        $string =~ s/0\*\(0.d0\) - /-/g;
      }
      #if ($string =~ m/ - 0.d0/) {
      #  $string =~ s/ - 0.d0//g;
      #}
      
      # Simplify terms raised to 1-power
      if ($string =~ m/\(PB\(1\)\*\*1\)/) {
        $string =~ s/\(PB\(1\)\*\*1\)/PB\(1\)/g;
      }
      if ($string =~ m/\(PB\(2\)\*\*1\)/) {
        $string =~ s/\(PB\(2\)\*\*1\)/PB\(2\)/g;
      }
      if ($string =~ m/\(PB\(3\)\*\*1\)/) {
        $string =~ s/\(PB\(3\)\*\*1\)/PB\(3\)/g;
      }
      if ($string =~ m/\(PA\(1\)\*\*1\)/) {
        $string =~ s/\(PA\(1\)\*\*1\)/PA\(1\)/g;
      }
      if ($string =~ m/\(PA\(2\)\*\*1\)/) {
        $string =~ s/\(PA\(2\)\*\*1\)/PA\(2\)/g;
      }
      if ($string =~ m/\(PA\(3\)\*\*1\)/) {
        $string =~ s/\(PA\(3\)\*\*1\)/PA\(3\)/g;
      }
      if ($string =~ m/\(1\/\(\(2\*Y\)\*\*1\)\)/) {
        $string =~ s/\(1\/\(\(2\*Y\)\*\*1\)\)/\(1\/\(2\*Y\)\)/g;
      }

      # Get rid of redundant 1's.
      while($string =~ m/1\*1/) {
        $string =~ s/1\*1/1/;
      }
      if ($string =~ m/\(1\*\(/) {
        $string =~ s/\(1\*\(/\(\(/g;
      }
      if ($string =~ m/\)\*1\)/) {
        $string =~ s/\)\*1\)/\)\)/g;
      }
      if ($string =~ m/\*1\*/) {
        $string =~ s/\*1\*/\*/g;
      }
      if ($string =~ m/\*1 /) {
        $string =~ s/\*1 / /g;
      }
      if ($string =~ m/ 1\*/) {
        $string =~ s/ 1\*/ /g;
      }
      if ($string =~ m/1\*P/) {
        $string =~ s/1\*P/P/g;
      }
      if ($string =~ m/\*\(1\)/) {
        $string =~ s/\*\(1\)//g;
      }
      if ($string =~ m/\*\(\(1\)\*/) {
        $string =~ s/\*\(\(1\)\*/\*\(/g;
      }
      if ($string =~ m/\*\(\(1\)\)/) {
        $string =~ s/\*\(\(1\)\)//g;
      }
      if ($string =~ m/\(\(1\)\*/) {
        $string =~ s/\(\(1\)\*/\(/g;
      }

      #**************************************************************
      # End of RegEXP portion
      #**************************************************************
      
      # Now, process the output string into 80 character-long lines.
      my $portion = "";
      my $initial = 0;
      my $final = 78;
      if (length($string) > 1) {
        while ($final < length($string)) {
          $portion .= substr($string, $initial, 79);
          $portion .= "&\n";
          $portion .= "&";
          $initial += 79;
          $final += 79;
        }
        $portion .= substr($string, $initial, length($string));
      }
      $portion .= "\n\n";
      print PRI $portion;
      $string = ""; # Reset the string before starting the loops again to prevent
                    # the solutions from concatonating onto each other.
    }
  }
  return $string;
}

# This subroutine uses the 2-center OL functionality
# to construct the Y momentum matrix formulas.  These integrals 
# are weighted sums of 2-center OL integrals.  This
# uses Cook's formalism (referenced above).
sub Mom_y {
  # Define passed parameters
  my $orb1 = $_[0];
  my $orb2 = $_[1];
  # Define local variables
  my @eTriad1;
  my @eTriad2;
  my $orbital1;
  my $orbital2;
  my $l1;
  my $l2;
  my $string = "";
  my $dim;
  #my $KE_CodeFile = "KE_Intg";
  my $num;

  # Open the output file
  #open (MOX,">$KE_CodeFile") || die "Error opening file.\n";

  &triads; # Initialize triads

  foreach $orbital2 (1..$orb2) {
    foreach $orbital1 (1..$orb1) {
      @eTriad1[0..2] = @{$eTriad[$orbital1]};
      @eTriad2[0..2] = @{$eTriad[$orbital2]};
      $string .= "wo(2,$orbital1,$orbital2) = ";
      $string .= "preFactor*("; # Outer paren
      $string .= $eTriad2[1];
      $string .= "*";
      $eTriad2[1] -= 1;
      if ($eTriad2[1] >= 0) {
        $string .= "(";
        foreach $dim (0..2) {
          $l1 = $eTriad1[$dim];
          $l2 = $eTriad2[$dim];
          $string .= &twoCenterSum($l1,$l2,$dim);
          if ($dim != 2) {
            $string .= "*";
          }
        }
        $string .= ")";
        $string .= " - ";
      }
      else {
        $string .= "(0.d0) - ";
      }
      $eTriad2[1] += 2;
      $string .= "2*a2*(";
      foreach $dim (0..2) {
        $l1 = $eTriad1[$dim];
        $l2 = $eTriad2[$dim];
        $string .= &twoCenterSum($l1,$l2,$dim);
        if ($dim != 2) {
          $string .= "*";
        }
      }
      $string .= ")";      
      $string .= ")"; # Outer paren

      #**************************************************************
      # Start RegEXP portion here
      #**************************************************************
      # Eliminate terms raised to 0-power
      if ($string =~ m/\(PA\(1\)\*\*0\)/) {
        $string =~ s/\(PA\(1\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PA\(2\)\*\*0\)/) {
        $string =~ s/\(PA\(2\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PA\(3\)\*\*0\)/) {
        $string =~ s/\(PA\(3\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(1\)\*\*0\)/) {
        $string =~ s/\(PB\(1\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(2\)\*\*0\)/) {
        $string =~ s/\(PB\(2\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(3\)\*\*0\)/) {
        $string =~ s/\(PB\(3\)\*\*0\)/1/g;
      }
      if ($string =~m/\(1\/\(\(2\*Y\)\*\*0\)\)/) {
        $string =~ s/\(1\/\(\(2\*Y\)\*\*0\)\)/1/g;
      }

      # Get rid of some damn zeroes
      if ($string =~ m/ \+ 0/) {
        $string =~ s/ \+ 0//g;
      }
      if ($string =~ m/0\*\(0.d0\) - /) {
        $string =~ s/0\*\(0.d0\) - /-/g;
      }
      #if ($string =~ m/ - 0.d0/) {
      #  $string =~ s/ - 0.d0//g;
      #}
      
      # Simplify terms raised to 1-power
      if ($string =~ m/\(PB\(1\)\*\*1\)/) {
        $string =~ s/\(PB\(1\)\*\*1\)/PB\(1\)/g;
      }
      if ($string =~ m/\(PB\(2\)\*\*1\)/) {
        $string =~ s/\(PB\(2\)\*\*1\)/PB\(2\)/g;
      }
      if ($string =~ m/\(PB\(3\)\*\*1\)/) {
        $string =~ s/\(PB\(3\)\*\*1\)/PB\(3\)/g;
      }
      if ($string =~ m/\(PA\(1\)\*\*1\)/) {
        $string =~ s/\(PA\(1\)\*\*1\)/PA\(1\)/g;
      }
      if ($string =~ m/\(PA\(2\)\*\*1\)/) {
        $string =~ s/\(PA\(2\)\*\*1\)/PA\(2\)/g;
      }
      if ($string =~ m/\(PA\(3\)\*\*1\)/) {
        $string =~ s/\(PA\(3\)\*\*1\)/PA\(3\)/g;
      }
      if ($string =~ m/\(1\/\(\(2\*Y\)\*\*1\)\)/) {
        $string =~ s/\(1\/\(\(2\*Y\)\*\*1\)\)/\(1\/\(2\*Y\)\)/g;
      }

      # Get rid of redundant 1's.
      while($string =~ m/1\*1/) {
        $string =~ s/1\*1/1/;
      }
      if ($string =~ m/\(1\*\(/) {
        $string =~ s/\(1\*\(/\(\(/g;
      }
      if ($string =~ m/\)\*1\)/) {
        $string =~ s/\)\*1\)/\)\)/g;
      }
      if ($string =~ m/\*1\*/) {
        $string =~ s/\*1\*/\*/g;
      }
      if ($string =~ m/\*1 /) {
        $string =~ s/\*1 / /g;
      }
      if ($string =~ m/ 1\*/) {
        $string =~ s/ 1\*/ /g;
      }
      if ($string =~ m/1\*P/) {
        $string =~ s/1\*P/P/g;
      }
      if ($string =~ m/\*\(1\)/) {
        $string =~ s/\*\(1\)//g;
      }
      if ($string =~ m/\*\(\(1\)\*/) {
        $string =~ s/\*\(\(1\)\*/\*\(/g;
      }
      if ($string =~ m/\*\(\(1\)\)/) {
        $string =~ s/\*\(\(1\)\)//g;
      }
      if ($string =~ m/\(\(1\)\*/) {
        $string =~ s/\(\(1\)\*/\(/g;
      }

      #**************************************************************
      # End of RegEXP portion
      #**************************************************************
      
      # Now, process the output string into 80 character-long lines.
      my $portion = "";
      my $initial = 0;
      my $final = 78;
      if (length($string) > 1) {
        while ($final < length($string)) {
          $portion .= substr($string, $initial, 79);
          $portion .= "&\n";
          $portion .= "&";
          $initial += 79;
          $final += 79;
        }
        $portion .= substr($string, $initial, length($string));
      }
      $portion .= "\n\n";
      print PRI $portion;
      $string = ""; # Reset the string before starting the loops again to prevent
                    # the solutions from concatonating onto each other.
    }
  }
  return $string;
}

# This subroutine uses the 2-center OL functionality
# to construct the Z momentum matrix formulas.  These integrals 
# are weighted sums of 2-center OL integrals.  This
# uses Cook's formalism (referenced above).
sub Mom_z {
  # Define passed parameters
  my $orb1 = $_[0];
  my $orb2 = $_[1];
  # Define local variables
  my @eTriad1;
  my @eTriad2;
  my $orbital1;
  my $orbital2;
  my $l1;
  my $l2;
  my $string = "";
  my $dim;
  #my $KE_CodeFile = "KE_Intg";
  my $num;

  # Open the output file
  #open (MOX,">$KE_CodeFile") || die "Error opening file.\n";

  &triads; # Initialize triads

  foreach $orbital2 (1..$orb2) {
    foreach $orbital1 (1..$orb1) {
      @eTriad1[0..2] = @{$eTriad[$orbital1]};
      @eTriad2[0..2] = @{$eTriad[$orbital2]};
      $string .= "wo(3,$orbital1,$orbital2) = ";
      $string .= "preFactor*("; # Outer paren
      $string .= $eTriad2[2];
      $string .= "*";
      $eTriad2[2] -= 1;
      if ($eTriad2[2] >= 0) {
        $string .= "(";
        foreach $dim (0..2) {
          $l1 = $eTriad1[$dim];
          $l2 = $eTriad2[$dim];
          $string .= &twoCenterSum($l1,$l2,$dim);
          if ($dim != 2) {
            $string .= "*";
          }
        }
        $string .= ")";
        $string .= " - ";
      }
      else {
        $string .= "(0.d0) - ";
      }
      $eTriad2[2] += 2;
      $string .= "2*a2*(";
      foreach $dim (0..2) {
        $l1 = $eTriad1[$dim];
        $l2 = $eTriad2[$dim];
        $string .= &twoCenterSum($l1,$l2,$dim);
        if ($dim != 2) {
          $string .= "*";
        }
      }
      $string .= ")";      
      $string .= ")"; # Outer paren

      #**************************************************************
      # Start RegEXP portion here
      #**************************************************************
      # Eliminate terms raised to 0-power
      if ($string =~ m/\(PA\(1\)\*\*0\)/) {
        $string =~ s/\(PA\(1\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PA\(2\)\*\*0\)/) {
        $string =~ s/\(PA\(2\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PA\(3\)\*\*0\)/) {
        $string =~ s/\(PA\(3\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(1\)\*\*0\)/) {
        $string =~ s/\(PB\(1\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(2\)\*\*0\)/) {
        $string =~ s/\(PB\(2\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(3\)\*\*0\)/) {
        $string =~ s/\(PB\(3\)\*\*0\)/1/g;
      }
      if ($string =~m/\(1\/\(\(2\*Y\)\*\*0\)\)/) {
        $string =~ s/\(1\/\(\(2\*Y\)\*\*0\)\)/1/g;
      }

      # Get rid of some damn zeroes
      if ($string =~ m/ \+ 0/) {
        $string =~ s/ \+ 0//g;
      }
      if ($string =~ m/0\*\(0.d0\) - /) {
        $string =~ s/0\*\(0.d0\) - /-/g;
      }
      #if ($string =~ m/ - 0.d0/) {
      #  $string =~ s/ - 0.d0//g;
      #}
      
      # Simplify terms raised to 1-power
      if ($string =~ m/\(PB\(1\)\*\*1\)/) {
        $string =~ s/\(PB\(1\)\*\*1\)/PB\(1\)/g;
      }
      if ($string =~ m/\(PB\(2\)\*\*1\)/) {
        $string =~ s/\(PB\(2\)\*\*1\)/PB\(2\)/g;
      }
      if ($string =~ m/\(PB\(3\)\*\*1\)/) {
        $string =~ s/\(PB\(3\)\*\*1\)/PB\(3\)/g;
      }
      if ($string =~ m/\(PA\(1\)\*\*1\)/) {
        $string =~ s/\(PA\(1\)\*\*1\)/PA\(1\)/g;
      }
      if ($string =~ m/\(PA\(2\)\*\*1\)/) {
        $string =~ s/\(PA\(2\)\*\*1\)/PA\(2\)/g;
      }
      if ($string =~ m/\(PA\(3\)\*\*1\)/) {
        $string =~ s/\(PA\(3\)\*\*1\)/PA\(3\)/g;
      }
      if ($string =~ m/\(1\/\(\(2\*Y\)\*\*1\)\)/) {
        $string =~ s/\(1\/\(\(2\*Y\)\*\*1\)\)/\(1\/\(2\*Y\)\)/g;
      }

      # Get rid of redundant 1's.
      while($string =~ m/1\*1/) {
        $string =~ s/1\*1/1/;
      }
      if ($string =~ m/\(1\*\(/) {
        $string =~ s/\(1\*\(/\(\(/g;
      }
      if ($string =~ m/\)\*1\)/) {
        $string =~ s/\)\*1\)/\)\)/g;
      }
      if ($string =~ m/\*1\*/) {
        $string =~ s/\*1\*/\*/g;
      }
      if ($string =~ m/\*1 /) {
        $string =~ s/\*1 / /g;
      }
      if ($string =~ m/ 1\*/) {
        $string =~ s/ 1\*/ /g;
      }
      if ($string =~ m/1\*P/) {
        $string =~ s/1\*P/P/g;
      }
      if ($string =~ m/\*\(1\)/) {
        $string =~ s/\*\(1\)//g;
      }
      if ($string =~ m/\*\(\(1\)\*/) {
        $string =~ s/\*\(\(1\)\*/\*\(/g;
      }
      if ($string =~ m/\*\(\(1\)\)/) {
        $string =~ s/\*\(\(1\)\)//g;
      }
      if ($string =~ m/\(\(1\)\*/) {
        $string =~ s/\(\(1\)\*/\(/g;
      }

      #**************************************************************
      # End of RegEXP portion
      #**************************************************************
      
      # Now, process the output string into 80 character-long lines.
      my $portion = "";
      my $initial = 0;
      my $final = 78;
      if (length($string) > 1) {
        while ($final < length($string)) {
          $portion .= substr($string, $initial, 79);
          $portion .= "&\n";
          $portion .= "&";
          $initial += 79;
          $final += 79;
        }
        $portion .= substr($string, $initial, length($string));
      }
      $portion .= "\n\n";
      print PRI $portion;
      $string = ""; # Reset the string before starting the loops again to prevent
                    # the solutions from concatonating onto each other.
    }
  }
  return $string;
}

# This subroutine generates the matrix elements g(16,16)
sub print_g {
  # Define passed parameters
  my $orb1 = $_[0];
  my $orb2 = $_[1];
  # Define local variables.
  my $string;
  my @g;
  my $index1;
  my $index2;

  $g[1][1] =  "g(i,1,1) = +wo(i,1,1)";
  $g[1][2] =  "g(i,1,2) = +wo(i,1,2)";
  $g[1][3] =  "g(i,1,3) = +wo(i,1,3)";
  $g[1][4] =  "g(i,1,4) = +wo(i,1,4)";
  $g[1][5] =  "g(i,1,5) = +wo(i,1,8)";
  $g[1][6] =  "g(i,1,6) = +wo(i,1,9)";
  $g[1][7] =  "g(i,1,7) = +wo(i,1,10)";
  $g[1][8] =  "g(i,1,8) = +wo(i,1,5)-wo(i,1,6)";
  $g[1][9] =  "g(i,1,9) = +2*wo(i,1,7) - wo(i,1,5)-wo(i,1,6)";
  $g[1][10] = "g(i,1,10) = +wo(i,1,11)";
  $g[1][11] = "g(i,1,11) = +wo(i,1,13) - wo(i,1,14)";
  $g[1][12] = "g(i,1,12) = +wo(i,1,18) - 3*wo(i,1,15)";
  $g[1][13] = "g(i,1,13) = +3*wo(i,1,12) - wo(i,1,19)";
  $g[1][14] = "g(i,1,14) = +2*wo(i,1,20) - 3*wo(i,1,13)-3*wo(i,1,14)";
  $g[1][15] = "g(i,1,15) = +4*wo(i,1,16) - wo(i,1,18)-wo(i,1,15)";
  $g[1][16] = "g(i,1,16) = +4*wo(i,1,17) - wo(i,1,12)-wo(i,1,19)";
  $g[2][1] =  "g(i,2,1) = +wo(i,2,1)";
  $g[2][2] =  "g(i,2,2) = +wo(i,2,2)";
  $g[2][3] =  "g(i,2,3) = +wo(i,2,3)";
  $g[2][4] =  "g(i,2,4) = +wo(i,2,4)";
  $g[2][5] =  "g(i,2,5) = +wo(i,2,8)";
  $g[2][6] =  "g(i,2,6) = +wo(i,2,9)";
  $g[2][7] =  "g(i,2,7) = +wo(i,2,10)";
  $g[2][8] =  "g(i,2,8) = +wo(i,2,5) - wo(i,2,6)";
  $g[2][9] =  "g(i,2,9) = +2*wo(i,2,7) - wo(i,2,5)-wo(i,2,6)";
  $g[2][10] = "g(i,2,10) = +wo(i,2,11)";
  $g[2][11] = "g(i,2,11) = +wo(i,2,13) - wo(i,2,14)";
  $g[2][12] = "g(i,2,12) = +wo(i,2,18) - 3*wo(i,2,15)";
  $g[2][13] = "g(i,2,13) = +3*wo(i,2,12) - wo(i,2,19)";
  $g[2][14] = "g(i,2,14) = +2*wo(i,2,20) - 3*wo(i,2,13) - 3*wo(i,2,14)";
  $g[2][15] = "g(i,2,15) = +4*wo(i,2,16) - wo(i,2,18) - wo(i,2,15)";
  $g[2][16] = "g(i,2,16) = +4*wo(i,2,17) - wo(i,2,12) - wo(i,2,19)";
  $g[3][1] =  "g(i,3,1) = +wo(i,3,1)";
  $g[3][2] =  "g(i,3,2) = +wo(i,3,2)";
  $g[3][3] =  "g(i,3,3) = +wo(i,3,3)";
  $g[3][4] =  "g(i,3,4) = +wo(i,3,4)";
  $g[3][5] =  "g(i,3,5) = +wo(i,3,8)";
  $g[3][6] =  "g(i,3,6) = +wo(i,3,9)";
  $g[3][7] =  "g(i,3,7) = +wo(i,3,10)";
  $g[3][8] =  "g(i,3,8) = +wo(i,3,5) - wo(i,3,6)";
  $g[3][9] =  "g(i,3,9) = +2*wo(i,3,7) - wo(i,3,5)-wo(i,3,6)";
  $g[3][10] = "g(i,3,10) = +wo(i,3,11)";
  $g[3][11] = "g(i,3,11) = +wo(i,3,13) - wo(i,3,14)";
  $g[3][12] = "g(i,3,12) = +wo(i,3,18) - 3*wo(i,3,15)";
  $g[3][13] = "g(i,3,13) = +3*wo(i,3,12) - wo(i,3,19)";
  $g[3][14] = "g(i,3,14) = +2*wo(i,3,20) - 3*wo(i,3,13) - 3*wo(i,3,14)";
  $g[3][15] = "g(i,3,15) = +4*wo(i,3,16) - wo(i,3,18) - wo(i,3,15)";
  $g[3][16] = "g(i,3,16) = +4*wo(i,3,17) - wo(i,3,12) - wo(i,3,19)";
  $g[4][1] =  "g(i,4,1) = +wo(i,4,1)";
  $g[4][2] =  "g(i,4,2) = +wo(i,4,2)";
  $g[4][3] =  "g(i,4,3) = +wo(i,4,3)";
  $g[4][4] =  "g(i,4,4) = +wo(i,4,4)";
  $g[4][5] =  "g(i,4,5) = +wo(i,4,8)";
  $g[4][6] =  "g(i,4,6) = +wo(i,4,9)";
  $g[4][7] =  "g(i,4,7) = +wo(i,4,10)";
  $g[4][8] =  "g(i,4,8) = +wo(i,4,5) - wo(i,4,6)";
  $g[4][9] =  "g(i,4,9) = +2*wo(i,4,7) - wo(i,4,5) - wo(i,4,6)";
  $g[4][10] = "g(i,4,10) = +wo(i,4,11)";
  $g[4][11] = "g(i,4,11) = +wo(i,4,13) - wo(i,4,14)";
  $g[4][12] = "g(i,4,12) = +wo(i,4,18) - 3*wo(i,4,15)";
  $g[4][13] = "g(i,4,13) = +3*wo(i,4,12) - wo(i,4,19)";
  $g[4][14] = "g(i,4,14) = +2*wo(i,4,20) - 3*wo(i,4,13) - 3*wo(i,4,14)";
  $g[4][15] = "g(i,4,15) = +4*wo(i,4,16) - wo(i,4,18) - wo(i,4,15)";
  $g[4][16] = "g(i,4,16) = +4*wo(i,4,17) - wo(i,4,12) - wo(i,4,19)";
  $g[5][1] =  "g(i,5,1) = +wo(i,8,1)";
  $g[5][2] =  "g(i,5,2) = +wo(i,8,2)";
  $g[5][3] =  "g(i,5,3) = +wo(i,8,3)";
  $g[5][4] =  "g(i,5,4) = +wo(i,8,4)";
  $g[5][5] =  "g(i,5,5) = +wo(i,8,8)";
  $g[5][6] =  "g(i,5,6) = +wo(i,8,9)";
  $g[5][7] =  "g(i,5,7) = +wo(i,8,10)";
  $g[5][8] =  "g(i,5,8) = +wo(i,8,5) - wo(i,8,6)";
  $g[5][9] =  "g(i,5,9) = +2*wo(i,8,7) - wo(i,8,5) - wo(i,8,6)";
  $g[5][10] = "g(i,5,10) = +wo(i,8,11)";
  $g[5][11] = "g(i,5,11) = +wo(i,8,13) - wo(i,8,14)";
  $g[5][12] = "g(i,5,12) = +wo(i,8,18) - 3*wo(i,8,15)";
  $g[5][13] = "g(i,5,13) = +3*wo(i,8,12) - wo(i,8,19)";
  $g[5][14] = "g(i,5,14) = +2*wo(i,8,20) - 3*wo(i,8,13) - 3*wo(i,8,14)";
  $g[5][15] = "g(i,5,15) = +4*wo(i,8,16) - wo(i,8,18) - wo(i,8,15)";
  $g[5][16] = "g(i,5,16) = +4*wo(i,8,17) - wo(i,8,12) - wo(i,8,19)";
  $g[6][1] =  "g(i,6,1) = +wo(i,9,1)";
  $g[6][2] =  "g(i,6,2) = +wo(i,9,2)";
  $g[6][3] =  "g(i,6,3) = +wo(i,9,3)";
  $g[6][4] =  "g(i,6,4) = +wo(i,9,4)";
  $g[6][5] =  "g(i,6,5) = +wo(i,9,8)";
  $g[6][6] =  "g(i,6,6) = +wo(i,9,9)";
  $g[6][7] =  "g(i,6,7) = +wo(i,9,10)";
  $g[6][8] =  "g(i,6,8) = +wo(i,9,5) - wo(i,9,6)";
  $g[6][9] =  "g(i,6,9) = +2*wo(i,9,7) - wo(i,9,5) - wo(i,9,6)";
  $g[6][10] = "g(i,6,10) = +wo(i,9,11)";
  $g[6][11] = "g(i,6,11) = +wo(i,9,13) - wo(i,9,14)";
  $g[6][12] = "g(i,6,12) = +wo(i,9,18) - 3*wo(i,9,15)";
  $g[6][13] = "g(i,6,13) = +3*wo(i,9,12) - wo(i,9,19)";
  $g[6][14] = "g(i,6,14) = +2*wo(i,9,20) - 3*wo(i,9,13) - 3*wo(i,9,14)";
  $g[6][15] = "g(i,6,15) = +4*wo(i,9,16) - wo(i,9,18) - wo(i,9,15)";
  $g[6][16] = "g(i,6,16) = +4*wo(i,9,17) - wo(i,9,12) - wo(i,9,19)";
  $g[7][1] =  "g(i,7,1) = +wo(i,10,1)";
  $g[7][2] =  "g(i,7,2) = +wo(i,10,2)";
  $g[7][3] =  "g(i,7,3) = +wo(i,10,3)";
  $g[7][4] =  "g(i,7,4) = +wo(i,10,4)";
  $g[7][5] =  "g(i,7,5) = +wo(i,10,8)";
  $g[7][6] =  "g(i,7,6) = +wo(i,10,9)";
  $g[7][7] =  "g(i,7,7) = +wo(i,10,10)";
  $g[7][8] =  "g(i,7,8) = +wo(i,10,5) - wo(i,10,6)";
  $g[7][9] =  "g(i,7,9) = +2*wo(i,10,7) - wo(i,10,5) - wo(i,10,6)";
  $g[7][10] = "g(i,7,10) = +wo(i,10,11)";
  $g[7][11] = "g(i,7,11) = +wo(i,10,13) - wo(i,10,14)";
  $g[7][12] = "g(i,7,12) = +wo(i,10,18) - 3*wo(i,10,15)";
  $g[7][13] = "g(i,7,13) = +3*wo(i,10,12) - wo(i,10,19)";
  $g[7][14] = "g(i,7,14) = +2*wo(i,10,20) - 3*wo(i,10,13) - 3*wo(i,10,14)";
  $g[7][15] = "g(i,7,15) = +4*wo(i,10,16) - wo(i,10,18) - wo(i,10,15)";
  $g[7][16] = "g(i,7,16) = +4*wo(i,10,17) - wo(i,10,12) - wo(i,10,19)";
  $g[8][1] =  "g(i,8,1) = +wo(i,5,1)-wo(i,6,1)";
  $g[8][2] =  "g(i,8,2) = +wo(i,5,2)-wo(i,6,2)";
  $g[8][3] =  "g(i,8,3) = +wo(i,5,3)-wo(i,6,3)";
  $g[8][4] =  "g(i,8,4) = +wo(i,5,4)-wo(i,6,4)";
  $g[8][5] =  "g(i,8,5) = +wo(i,5,8)-wo(i,6,8)";
  $g[8][6] =  "g(i,8,6) = +wo(i,5,9)-wo(i,6,9)";
  $g[8][7] =  "g(i,8,7) = +wo(i,5,10) - wo(i,6,10)";
  $g[8][8] =  "g(i,8,8) = +wo(i,5,5)-wo(i,5,6) - wo(i,6,5) + wo(i,6,6)";
  $g[8][9] = "g(i,8,9) = +2*wo(i,5,7) - wo(i,5,5) - wo(i,5,6)&\n";
  $g[8][9] .= "&- 2*wo(i,6,7) + wo(i,6,5) + wo(i,6,6)\n";
  $g[8][10] = "g(i,8,10) = +wo(i,5,11) - wo(i,6,11)";
  $g[8][11] = "g(i,8,11) = +wo(i,5,13) - wo(i,5,14) - wo(i,6,13) +wo(i,6,14)";
  $g[8][12] = "g(i,8,12) = +wo(i,5,18) - 3*wo(i,5,15) - wo(i,6,18) + 3*wo(i,6,15)";
  $g[8][13] = "g(i,8,13) = +3*wo(i,5,12) - wo(i,5,19) - 3*wo(i,6,12) + wo(i,6,19)";
  $g[8][14] = "g(i,8,14) = +2*wo(i,5,20) - 3*wo(i,5,13) - 3*wo(i,5,14)&\n";
  $g[8][14] .= "&- 2*wo(i,6,20) + 3*wo(i,6,13) +3*wo(i,6,14)\n";
  $g[8][15] = "g(i,8,15) = +4*wo(i,5,16) - wo(i,5,18) - wo(i,5,15)&\n";
  $g[8][15] .= "&- 4*wo(i,6,16) + wo(i,6,18) + wo(i,6,15)\n";
  $g[8][16] = "g(i,8,16) = +4*wo(i,5,17) - wo(i,5,12) - wo(i,5,19)&\n";
  $g[8][16] .= "&- 4*wo(i,6,17) + wo(i,6,12) +wo(i,6,19)\n";
  $g[9][1] = "g(i,9,1) = +2*wo(i,7,1) - wo(i,5,1) - wo(i,6,1)";
  $g[9][2] = "g(i,9,2) = +2*wo(i,7,2) - wo(i,5,2) - wo(i,6,2)";
  $g[9][3] = "g(i,9,3) = +2*wo(i,7,3) - wo(i,5,3) - wo(i,6,3)";
  $g[9][4] = "g(i,9,4) = +2*wo(i,7,4) - wo(i,5,4) - wo(i,6,4)";
  $g[9][5] = "g(i,9,5) = +2*wo(i,7,8) - wo(i,5,8) - wo(i,6,8)";
  $g[9][6] = "g(i,9,6) = +2*wo(i,7,9) - wo(i,5,9) - wo(i,6,9)";
  $g[9][7] = "g(i,9,7) = +2*wo(i,7,10) - wo(i,5,10) - wo(i,6,10)";
  $g[9][8] = "g(i,9,8) = +2*wo(i,7,5) - 2*wo(i,7,6) - wo(i,5,5)&\n";
  $g[9][8] .= "&+ wo(i,5,6) - wo(i,6,5) + wo(i,6,6)\n";
  $g[9][9] = "g(i,9,9) = +4*wo(i,7,7) - 2*wo(i,7,5) - 2*wo(i,7,6)&\n";
  $g[9][9] .= "&- 2*wo(i,5,7) + wo(i,5,5) + wo(i,5,6) - 2*wo(i,6,7) + wo(i,6,5) + wo(i,6,6)\n";
  $g[9][10] = "g(i,9,10) = +2*wo(i,7,11) - wo(i,5,11) - wo(i,6,11)";
  $g[9][11] = "g(i,9,11) = +2*wo(i,7,13) - 2*wo(i,7,14) - wo(i,5,13)&\n";
  $g[9][11] .= "&+ wo(i,5,14) - wo(i,6,13) + wo(i,6,14)\n";
  $g[9][12] = "g(i,9,12) = +2*wo(i,7,18) - 6*wo(i,7,15) - wo(i,5,18)&\n";
  $g[9][12] .= "&+ 3*wo(i,5,15) - wo(i,6,18) + 3*wo(i,6,15)\n";
  $g[9][13] = "g(i,9,13) = +6*wo(i,7,12) - 2*wo(i,7,19) - 3*wo(i,5,12)&\n";
  $g[9][13] .= "&+ wo(i,5,19) - 3*wo(i,6,12) + wo(i,6,19)\n";
  $g[9][14] = "g(i,9,14) = +4*wo(i,7,20) - 6*wo(i,7,13) - 6*wo(i,7,14)&\n";
  $g[9][14] .= "&- 2*wo(i,5,20) + 3*wo(i,5,13) + 3*wo(i,5,14) - 2*wo(i,6,20)&\n";
  $g[9][14] .= "&+ 3*wo(i,6,13) + 3*wo(i,6,14)\n";
  $g[9][15] = "g(i,9,15) = +8*wo(i,7,16) - 2*wo(i,7,18) - 2*wo(i,7,15)&\n";
  $g[9][15] .= "&- 4*wo(i,5,16) + wo(i,5,18) + wo(i,5,15) - 4*wo(i,6,16)&\n";
  $g[9][15] .= "&+ wo(i,6,18) + wo(i,6,15)\n";
  $g[9][16] = "g(i,9,16) = +8*wo(i,7,17) - 2*wo(i,7,12) - 2*wo(i,7,19)&\n";
  $g[9][16] .= "&- 4*wo(i,5,17) + wo(i,5,12) + wo(i,5,19) - 4*wo(i,6,17)&\n";
  $g[9][16] .= "&+ wo(i,6,12) + wo(i,6,19)\n";
  $g[10][1] =  "g(i,10,1) = +wo(i,11,1)";
  $g[10][2] =  "g(i,10,2) = +wo(i,11,2)";
  $g[10][3] =  "g(i,10,3) = +wo(i,11,3)";
  $g[10][4] =  "g(i,10,4) = +wo(i,11,4)";
  $g[10][5] =  "g(i,10,5) = +wo(i,11,8)";
  $g[10][6] =  "g(i,10,6) = +wo(i,11,9)";
  $g[10][7] =  "g(i,10,7) = +wo(i,11,10)";
  $g[10][8] =  "g(i,10,8) = +wo(i,11,5) - wo(i,11,6)";
  $g[10][9] =  "g(i,10,9) = +2*wo(i,11,7) - wo(i,11,5) - wo(i,11,6)";
  $g[10][10] = "g(i,10,10) = +wo(i,11,11)";
  $g[10][11] = "g(i,10,11) = +wo(i,11,13) - wo(i,11,14)";
  $g[10][12] = "g(i,10,12) = +wo(i,11,18) - 3*wo(i,11,15)";
  $g[10][13] = "g(i,10,13) = +3*wo(i,11,12) - wo(i,11,19)";
  $g[10][14] = "g(i,10,14) = +2*wo(i,11,20) - 3*wo(i,11,13) - 3*wo(i,11,14)";
  $g[10][15] = "g(i,10,15) = +4*wo(i,11,16) - wo(i,11,18) - wo(i,11,15)";
  $g[10][16] = "g(i,10,16) = +4*wo(i,11,17) - wo(i,11,12) - wo(i,11,19)";
  $g[11][1] =  "g(i,11,1) = +wo(i,13,1) - wo(i,14,1)";
  $g[11][2] =  "g(i,11,2) = +wo(i,13,2) - wo(i,14,2)";
  $g[11][3] =  "g(i,11,3) = +wo(i,13,3) - wo(i,14,3)";
  $g[11][4] =  "g(i,11,4) = +wo(i,13,4) - wo(i,14,4)";
  $g[11][5] =  "g(i,11,5) = +wo(i,13,8) - wo(i,14,8)";
  $g[11][6] =  "g(i,11,6) = +wo(i,13,9) - wo(i,14,9)";
  $g[11][7] =  "g(i,11,7) = +wo(i,13,10) - wo(i,14,10)";
  $g[11][8] =  "g(i,11,8) = +wo(i,13,5) - wo(i,13,6)-wo(i,14,5) + wo(i,14,6)";
  $g[11][9] =  "g(i,11,9) = +2*wo(i,13,7) - wo(i,13,5) - wo(i,13,6)&\n";
  $g[11][9] .= "&- 2*wo(i,14,7) + wo(i,14,5) +wo(i,14,6)\n";
  $g[11][10] = "g(i,11,10) = +wo(i,13,11) - wo(i,14,11)";
  $g[11][11] = "g(i,11,11) = +wo(i,13,13) - wo(i,13,14) - wo(i,14,13) + wo(i,14,14)";
  $g[11][12] = "g(i,11,12) = +wo(i,13,18) - 3*wo(i,13,14) - wo(i,14,18) + 3*wo(i,14,15)";
  $g[11][13] = "g(i,11,13) = +3*wo(i,13,12) - wo(i,13,19) - 3*wo(i,14,12) + wo(i,14,19)";
  $g[11][14] = "g(i,11,14) = +2*wo(i,13,20) - 3*wo(i,13,13) - 3*wo(i,13,14)&\n";
  $g[11][14] .= "&- 2*wo(i,14,20) + 3*wo(i,14,13) + 3*wo(i,14,14)\n";
  $g[11][15] = "g(i,11,15) = +4*wo(i,13,16) - wo(i,13,18) - wo(i,13,15)&\n";
  $g[11][15] .= "&- 4*wo(i,14,16) + wo(i,14,18) + wo(i,14,15)\n";
  $g[11][16] = "g(i,11,16) = +4*wo(i,13,17) - wo(i,13,12) - wo(i,13,19)&\n";
  $g[11][16] .= "&- 4*wo(i,14,17) + wo(i,14,12) + wo(i,14,19)\n";
  $g[12][1] = "g(i,12,1) = +wo(i,18,1) - 3*wo(i,15,1)";
  $g[12][2] = "g(i,12,2) = +wo(i,18,2) - 3*wo(i,15,2)";
  $g[12][3] = "g(i,12,3) = +wo(i,18,3) - 3*wo(i,15,3)";
  $g[12][4] = "g(i,12,4) = +wo(i,18,4) - 3*wo(i,15,4)";
  $g[12][5] = "g(i,12,5) = +wo(i,18,8) - 3*wo(i,15,8)";
  $g[12][6] = "g(i,12,6) = +wo(i,18,9) - 3*wo(i,15,9)";
  $g[12][7] = "g(i,12,7) = +wo(i,18,10) - 3*wo(i,15,10)";
  $g[12][8] = "g(i,12,8) = +wo(i,18,5) - wo(i,18,6) - 3*wo(i,15,5) + 3*wo(i,15,6)";
  $g[12][9] = "g(i,12,9) = +2*wo(i,18,7) - wo(i,18,5)&\n";
  $g[12][9] .= "&- wo(i,18,6) - 6*wo(i,15,7) + 3*wo(i,15,5) + 3*wo(i,15,6)\n";
  $g[12][10] = "g(i,12,10) = +wo(i,18,11) - 3*wo(i,15,11)";
  $g[12][11] = "g(i,12,11) = +wo(i,18,13) - wo(i,18,14) - 3*wo(i,15,13) + 3*wo(i,15,14)  ";
  $g[12][12] = "g(i,12,12) = +wo(i,18,18) - 3*wo(i,18,15) - 3*wo(i,15,18) + 9*wo(i,15,15)";
  $g[12][13] = "g(i,12,13) = +3*wo(i,18,12) - wo(i,18,19) - 9*wo(i,15,12) + 3*wo(i,15,19)";
  $g[12][14] = "g(i,12,14) = +2*wo(i,18,20) - 3*wo(i,18,13) - 3*wo(i,18,14)&\n";
  $g[12][14] .= "&- 6*wo(i,15,20) + 9*wo(i,15,13) + 9*wo(i,15,14)\n";
  $g[12][15] = "g(i,12,15) = +4*wo(i,18,16) - wo(i,18,18) - wo(i,18,15)&\n";
  $g[12][15] .= "&- 12*wo(i,15,16) + 3*wo(i,15,18) + 3*wo(i,15,15)\n";
  $g[12][16] = "g(i,12,16) = +4*wo(i,18,17) - wo(i,18,12) - wo(i,18,19)&\n";
  $g[12][16] .= "&- 12*wo(i,15,17) + 3*wo(i,15,12) + 3*wo(i,15,19)\n";
  $g[13][1] = "g(i,13,1) = +3*wo(i,12,1) - wo(i,19,1)";
  $g[13][2] = "g(i,13,2) = +3*wo(i,12,2) - wo(i,19,2)";
  $g[13][3] = "g(i,13,3) = +3*wo(i,12,3) - wo(i,19,3)";
  $g[13][4] = "g(i,13,4) = +3*wo(i,12,4) - wo(i,19,4)";
  $g[13][5] = "g(i,13,5) = +3*wo(i,12,8) - wo(i,19,8)";
  $g[13][6] = "g(i,13,6) = +3*wo(i,12,9) - wo(i,19,9)";
  $g[13][7] = "g(i,13,7) = +3*wo(i,12,10) - wo(i,19,10)";
  $g[13][8] = "g(i,13,8) = +3*wo(i,12,5) - 3*wo(i,12,6) - wo(i,19,5) +wo(i,19,6)";
  $g[13][9] = "g(i,13,9) = +6*wo(i,12,7) - 3*wo(i,12,5)&\n";
  $g[13][9] .= "&- 3*wo(i,12,6) - 2*wo(i,19,7) + wo(i,19,5) +wo(i,19,6)\n";
  $g[13][10] = "g(i,13,10) = +3*wo(i,12,11) - wo(i,19,11)";
  $g[13][11] = "g(i,13,11) = +3*wo(i,12,13) - 3*wo(i,12,14) - wo(i,19,13) + wo(i,19,14)  ";
  $g[13][12] = "g(i,13,12) = +3*wo(i,12,18) - 9*wo(i,12,15) - wo(i,19,18) + 3*wo(i,19,15)";
  $g[13][13] = "g(i,13,13) = +9*wo(i,12,12) - 3*wo(i,12,19) - 3*wo(i,19,12) + wo(i,19,19)";
  $g[13][14] = "g(i,13,14) = +6*wo(i,12,20) - 9*wo(i,12,13) - 9*wo(i,12,14)&\n";
  $g[13][14] .= "&- 2*wo(i,19,20) + 3*wo(i,19,13) + 3*wo(i,19,14)\n";
  $g[13][15] = "g(i,13,15) = +12*wo(i,12,16) - 3*wo(i,12,18) - 3*wo(i,12,15)&\n";
  $g[13][15] .= "&- 4*wo(i,19,16) + wo(i,19,18) + wo(i,19,15)\n";
  $g[13][16] = "g(i,13,16) = +12*wo(i,12,17) - 3*wo(i,12,12) - 3*wo(i,12,19)&\n";
  $g[13][16] .= "&- 4*wo(i,19,17) + wo(i,19,12) + wo(i,19,19)\n";
  $g[14][1] = "g(i,14,1) = +2*wo(i,20,1) - 3*wo(i,13,1) - 3*wo(i,14,1)";
  $g[14][2] = "g(i,14,2) = +2*wo(i,20,2) - 3*wo(i,13,2) - 3*wo(i,14,2)";
  $g[14][3] = "g(i,14,3) = +2*wo(i,20,3) - 3*wo(i,13,3) - 3*wo(i,14,3)";
  $g[14][4] = "g(i,14,4) = +2*wo(i,20,4) - 3*wo(i,13,4) - 3*wo(i,14,4)";
  $g[14][5] = "g(i,14,5) = +2*wo(i,20,8) - 3*wo(i,13,8) - 3*wo(i,14,8)";
  $g[14][6] = "g(i,14,6) = +2*wo(i,20,9) - 3*wo(i,13,9) - 3*wo(i,14,9)";
  $g[14][7] = "g(i,14,7) = +2*wo(i,20,10) - 3*wo(i,13,10) - 3*wo(i,14,10)";
  $g[14][8] = "g(i,14,8) = +2*wo(i,20,5) - 2*wo(i,20,6) - 3*wo(i,13,5)&\n";
  $g[14][8] .= "&+ 3*wo(i,13,6) - 3*wo(i,14,5) +3*wo(i,14,6)\n";
  $g[14][9] = "g(i,14,9) = +4*wo(i,20,7) - 2*wo(i,20,5) - 2*wo(i,20,6)&\n";
  $g[14][9] .= "&- 6*wo(i,13,7) + 3*wo(i,13,5) + 3*wo(i,13,6)&\n";
  $g[14][9] .= "&- 6*wo(i,14,7) + 3*wo(i,14,5) + 3*wo(i,14,6)\n";
  $g[14][10] = "g(i,14,10) = +2*wo(i,20,11) - 3*wo(i,13,11) - 3*wo(i,14,11)";
  $g[14][11] = "g(i,14,11) = +2*wo(i,20,13) - 2*wo(i,20,14) - 3*wo(i,13,13)&\n";
  $g[14][11] .= "&+ 3*wo(i,13,14) - 3*wo(i,14,13) + 3*wo(i,14,14)\n";
  $g[14][12] = "g(i,14,12) = +2*wo(i,20,18) - 6*wo(i,20,15) - 3*wo(i,13,18)&\n";
  $g[14][12] .= "&+ 9*wo(i,13,15) - 3*wo(i,14,18) + 9*wo(i,14,15)\n";
  $g[14][13] = "g(i,14,13) = +6*wo(i,20,12) - 2*wo(i,20,19) - 9*wo(i,13,12)&\n";
  $g[14][13] .= "&+ 3*wo(i,13,19) - 9*wo(i,14,12) + 3*wo(i,14,19)\n";
  $g[14][14] = "g(i,14,14) = +4*wo(i,20,20) - 6*wo(i,20,13) - 6*wo(i,20,14)&\n";
  $g[14][14] .= "&- 6*wo(i,13,20) + 9*wo(i,13,13) + 9*wo(i,13,14)&\n";
  $g[14][14] .= "&- 6*wo(i,14,20) + 9*wo(i,15,13) + 9*wo(i,15,14)\n";
  $g[14][15] = "g(i,14,15) = +8*wo(i,20,16) - 2*wo(i,20,18) - 2*wo(i,20,15)&\n";
  $g[14][15] .= "&- 12*wo(i,13,16) + 3*wo(i,13,18) + 3*wo(i,13,15)&\n";
  $g[14][15] .= "&- 12*wo(i,15,16) + 3*wo(i,15,18) + 3*wo(i,15,15)\n";
  $g[14][16] = "g(i,14,16) = +8*wo(i,20,17) - 2*wo(i,20,12) - 2*wo(i,20,19)&\n";
  $g[14][16] .= "&- 12*wo(i,13,17) + 3*wo(i,13,12) + 3*wo(i,13,19)&\n";
  $g[14][16] .= "&- 12*wo(i,14,17) + 3*wo(i,14,12) + 3*wo(i,14,19)\n";
  $g[15][1] = "g(i,15,1) = +4*wo(i,16,1) - wo(i,18,1) - wo(i,15,1)";
  $g[15][2] = "g(i,15,2) = +4*wo(i,16,2) - wo(i,18,2) - wo(i,15,2)";
  $g[15][3] = "g(i,15,3) = +4*wo(i,16,3) - wo(i,18,3) - wo(i,15,3)";
  $g[15][4] = "g(i,15,4) = +4*wo(i,16,4) - wo(i,18,4) - wo(i,15,4)";
  $g[15][5] = "g(i,15,5) = +4*wo(i,16,8) - wo(i,18,8) - wo(i,15,8)";
  $g[15][6] = "g(i,15,6) = +4*wo(i,16,9) - wo(i,18,9) - wo(i,15,9)";
  $g[15][7] = "g(i,15,7) = +4*wo(i,16,10) - wo(i,18,10) - wo(i,15,10)";
  $g[15][8] = "g(i,15,8) = +4*wo(i,16,5) - 4*wo(i,16,6) - wo(i,18,5)&\n";
  $g[15][8] .= "&+ wo(i,18,6) - wo(i,15,5) + wo(i,15,6)\n";
  $g[15][9] = "g(i,15,9) = +8*wo(i,16,7) - 4*wo(i,16,5) - 4*wo(i,16,6)&\n";
  $g[15][9] .= "&- 2*wo(i,18,7) + wo(i,18,5) + wo(i,18,6)&\n";
  $g[15][9] .= "&- 2*wo(i,15,7) + wo(i,15,5)+wo(i,15,6)\n";
  $g[15][10] = "g(i,15,10) = +4*wo(i,16,11) - wo(i,18,11)-wo(i,15,11)";
  $g[15][11] = "g(i,15,11) = +4*wo(i,16,13) - 4*wo(i,16,14) - wo(i,18,13)&\n";
  $g[15][11] .= "&+ wo(i,18,14) - wo(i,15,13) + wo(i,15,14)\n";
  $g[15][12] = "g(i,15,12) = +4*wo(i,16,18) - 12*wo(i,16,15) - wo(i,18,18)&\n";
  $g[15][12] .= "&+ 3*wo(i,18,15) - wo(i,15,18) + 3*wo(i,15,15)\n";
  $g[15][13] = "g(i,15,13) = +12*wo(i,16,12) - 4*wo(i,16,19) - 3*wo(i,18,12)&\n";
  $g[15][13] .= "&+ wo(i,18,19) - 3*wo(i,15,12) + wo(i,15,19)\n";
  $g[15][14] = "g(i,15,14) = +8*wo(i,16,20) - 12*wo(i,16,13) - 12*wo(i,16,14)&\n";
  $g[15][14] .= "&- 2*wo(i,18,20) + 3*wo(i,18,13) + 3*wo(i,18,14)&\n";
  $g[15][14] .= "&- 2*wo(i,15,20) + 3*wo(i,15,13) + 3*wo(i,15,14)\n";
  $g[15][15] = "g(i,15,15) = +16*wo(i,16,16) - 4*wo(i,16,18) - 4*wo(i,16,15)&\n";
  $g[15][15] .= "&- 4*wo(i,18,16) + wo(i,18,18) + wo(i,18,15)&\n";
  $g[15][15] .= "&- 4*wo(i,15,16) + wo(i,14,18) + wo(i,14,15)";
  $g[15][16] = "g(i,15,16) = +16*wo(i,16,17) - 4*wo(i,16,12) - 4*wo(i,16,19)&\n";
  $g[15][16] .= "&- 4*wo(i,18,17) + wo(i,18,12) + wo(i,18,19)&\n";
  $g[15][16] .= "&- 4*wo(i,15,17) + wo(i,15,12) + wo(i,15,19)\n";
  $g[16][1] = "g(i,16,1) = +4*wo(i,17,1) - wo(i,12,1) - wo(i,19,1)";
  $g[16][2] = "g(i,16,2) = +4*wo(i,17,2) - wo(i,12,2) - wo(i,19,2)";
  $g[16][3] = "g(i,16,3) = +4*wo(i,17,3) - wo(i,12,3) - wo(i,19,3)";
  $g[16][4] = "g(i,16,4) = +4*wo(i,17,4) - wo(i,12,4) - wo(i,19,4)";
  $g[16][5] = "g(i,16,5) = +4*wo(i,17,8) - wo(i,12,8) - wo(i,19,8)";
  $g[16][6] = "g(i,16,6) = +4*wo(i,17,9) - wo(i,12,9) - wo(i,19,9)";
  $g[16][7] = "g(i,16,7) = +4*wo(i,17,10)  -wo(i,12,10) - wo(i,19,10)";
  $g[16][8] = "g(i,16,8) = +4*wo(i,17,5) - 4*wo(i,17,6) - wo(i,12,5)&\n";
  $g[16][8] .= "&+ wo(i,12,6) - wo(i,19,5) + wo(i,19,6)\n";
  $g[16][9] = "g(i,16,9) = +8*wo(i,17,7) - 4*wo(i,17,5) - 4*wo(i,17,6)&\n";
  $g[16][9] .= "&- 2*wo(i,12,7) + wo(i,12,5) + wo(i,12,6)&\n";
  $g[16][9] .= "&- 2*wo(i,19,7) + wo(i,19,5) + wo(i,19,6)\n";
  $g[16][10] = "g(i,16,10) = +4*wo(i,17,11) - wo(i,12,11) - wo(i,19,11)";
  $g[16][11] = "g(i,16,11) = +4*wo(i,17,13) - 4*wo(i,17,14) - wo(i,12,13)&\n";
  $g[16][11] .= "&+ wo(i,12,14) - wo(i,19,13) + wo(i,19,14)\n";
  $g[16][12] = "g(i,16,12) = +4*wo(i,17,18) - 12*wo(i,17,15) - wo(i,12,18)&\n";
  $g[16][12] .= "&+ 3*wo(i,12,15) - wo(i,19,18) + 3*wo(i,19,15)\n";
  $g[16][13] = "g(i,16,13) = +12*wo(i,17,12) - 4*wo(i,17,19) - 3*wo(i,12,12)&\n";
  $g[16][13] .= "&+ wo(i,12,19) - 3*wo(i,19,12) + wo(i,19,19)\n";
  $g[16][14] = "g(i,16,14) = +8*wo(i,17,20) - 12*wo(i,17,13) - 12*wo(i,17,14)&\n";
  $g[16][14] .= "&- 2*wo(i,12,20) + 3*wo(i,12,13) + 3*wo(i,12,14)&\n";
  $g[16][14] .= "&- 2*wo(i,19,20) + 3*wo(i,19,13) + 3*wo(i,19,14)\n";
  $g[16][15] = "g(i,16,15) = +16*wo(i,17,16) - 4*wo(i,17,18) - 4*wo(i,17,15)&\n";
  $g[16][15] .= "&- 4*wo(i,12,16) + wo(i,12,18) + wo(i,12,15)&\n";
  $g[16][15] .= "&- 4*wo(i,19,16) + wo(i,19,18) + wo(i,19,15)\n";
  $g[16][16] = "g(i,16,16) = +16*wo(i,17,17) - 4*wo(i,17,12) - 4*wo(i,17,19)&\n";
  $g[16][16] .= "&- 4*wo(i,12,17) + wo(i,12,12) + wo(i,12,19)&\n";
  $g[16][16] .= "&- 4*wo(i,19,17) + wo(i,19,12) + wo(i,19,19)\n";
    
  $string = "";
  foreach $index2 (1..$orb2) {
    foreach $index1 (1..$orb1) {
      $string .= $g[$index1][$index2];
      $string .= "\n";
    }
  }
  return $string;
}

# This subroutine prints the momentum matrix integral
# solutions in the appropriate order, with l1l2switch
# functionality included.  The usual 16x16 g matrix is 
# actually 16x16x3 for this situation.
sub print_Mom {
  # Declare local variables
  my $string;
  my $momentumCodeFile = "momIntg";

  # Open the output file
  open (PRI,">$momentumCodeFile") || die "Error opening file.\n";

  print PRI "if (l1l2switch.eq.34) then\n\n";
  $string = &Mom_x(4,4);
  $string = &Mom_y(4,4);
  $string = &Mom_z(4,4);
  $string .= &print_g(4,4);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.17) then\n\n";
  $string = &Mom_x(1,1);
  $string = &Mom_y(1,1);
  $string = &Mom_z(1,1);
  $string .= &print_g(1,1);
  print PRI "$string\n";
 
  print PRI "else if (l1l2switch.eq.33) then\n\n";
  $string = &Mom_x(1,4);
  $string = &Mom_y(1,4);
  $string = &Mom_z(1,4);
  $string .= &print_g(1,4);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.65) then\n\n";
  $string = &Mom_x(1,10);
  $string = &Mom_y(1,10);
  $string = &Mom_z(1,10);
  $string .= &print_g(1,9);
  print PRI "$string\n";
 
  print PRI "else if (l1l2switch.eq.129) then\n\n";
  $string = &Mom_x(1,20);
  $string = &Mom_y(1,20);
  $string = &Mom_z(1,20);
  $string .= &print_g(1,16);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.18) then\n\n";
  $string = &Mom_x(4,1);
  $string = &Mom_y(4,1);
  $string = &Mom_z(4,1);
  $string .= &print_g(4,1);
  print PRI "$string\n";
 
  print PRI "else if (l1l2switch.eq.66) then\n\n";
  $string = &Mom_x(4,10);
  $string = &Mom_y(4,10);
  $string = &Mom_z(4,10);
  $string .= &print_g(4,9);
  print PRI "$string\n";
 
  print PRI "else if (l1l2switch.eq.130) then\n\n";
  $string = &Mom_x(4,20);
  $string = &Mom_y(4,20);
  $string = &Mom_z(4,20);
  $string .= &print_g(4,16);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.20) then\n\n";
  $string = &Mom_x(10,1);
  $string = &Mom_y(10,1);
  $string = &Mom_z(10,1);
  $string .= &print_g(9,1);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.36) then\n\n";
  $string = &Mom_x(10,4);
  $string = &Mom_y(10,4);
  $string = &Mom_z(10,4);
  $string .= &print_g(9,4);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.68) then\n\n";
  $string = &Mom_x(10,10);
  $string = &Mom_y(10,10);
  $string = &Mom_z(10,10);
  $string .= &print_g(9,9);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.132) then\n\n";
  $string = &Mom_x(10,20);
  $string = &Mom_y(10,20);
  $string = &Mom_z(10,20);
  $string .= &print_g(9,16);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.24) then\n\n";
  $string = &Mom_x(20,1);
  $string = &Mom_y(20,1);
  $string = &Mom_z(20,1);
  $string .= &print_g(16,1);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.40) then\n\n";
  $string = &Mom_x(20,4);
  $string = &Mom_y(20,4);
  $string = &Mom_z(20,4);
  $string .= &print_g(16,4);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.72) then\n\n";
  $string = &Mom_x(20,10);
  $string = &Mom_y(20,10);
  $string = &Mom_z(20,10);
  $string .= &print_g(16,9);
  print PRI "$string\n";
    
  print PRI "else\n\n";
  $string = &Mom_x(20,20);
  $string = &Mom_y(20,20);
  $string = &Mom_z(20,20);
  $string .= &print_g(16,16);
  print PRI "$string";
  print PRI "end if";
}


