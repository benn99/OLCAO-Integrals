#!use/bin/perl

##################################################################
#
# Author: Ben Walker
#
# This program computes the solution formulas for 3-center nuclear
# attraction integrals.  Positions 1 and 2 involve s-,  p-, d-, 
# or f-type orbitals, while the third (nuclear) position 
# is an s-type Gaussian function.  This formalism is readily 
# extendable to orbitals higher than f-type (g-type orbitals, etc.).
# The azimuthal quantum numbers are arranged according to the 
# "triad" formalism established by Browne & Poshusta in 1962, 
# the reference for which is the following:
# 
# ****************************************************************
# Browne JC, Poshusta RD., "Quantum-Mechanical Integrals 
# over Gaussian Atomic Orbitals.", The Journal of Chemical Physics, 
# Vol. 36, Iss. 7, P. 19337, (1962).
# ****************************************************************
#
# The formalism for the solution to the integrals comes from the
# following source:
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
&printNuc;
&incompleteGamma;
&altGamma;
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

# This subroutine constructs the nuclear potential
# integral solution.
sub nucOverlap {
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
  my $m1;
  my $m2;
  my $n1;
  my $n2;
  my $l;
  my $m;
  my $n;
  my $r;
  my $s;
  my $t;
  my $i;
  my $j;
  my $k;
  my $u;
  my $v;
  my $w;
  my $N;
  my $string = "";

  &triads; # Initialize triads

  foreach $orbital2 (1..$orb2) {
    foreach $orbital1 (1..$orb1) {
      @eTriad1[0..2] = @{$eTriad[$orbital1]};
      @eTriad2[0..2] = @{$eTriad[$orbital2]};
      $string .= "wo($orbital1,$orbital2) = ";
      $string .= "preFactor*("; # Outer enclosing paren
      $l1 = $eTriad1[0];
      $l2 = $eTriad2[0];
      $m1 = $eTriad1[1];
      $m2 = $eTriad2[1];
      $n1 = $eTriad1[2];
      $n2 = $eTriad2[2];
      my $num1;
      my $num2;
      my $num3;
      my $num4;
      my $num5;
      my $num6;
      my $num7;
      my $num8;
      my $num9;
      my $num10;
      my $num11;
      my $num12;
      foreach $l (0..($l1 + $l2)) {
        foreach $r (0..floor($l/2)) {
          foreach $i (0..floor(($l - 2*$r)/2)) {
            foreach $m (0..($m1 + $m2)) {
              foreach $s (0..floor($m/2)) {
                foreach $j (0..floor(($m - 2*$s)/2)) {
                  foreach $n (0..($n1 + $n2)) {
                    foreach $t (0..floor($n/2)) {
                      foreach $k (0..floor(($n - 2*$t)/2)) {
                        foreach $u (0..$l) {
                          foreach $v (0..$m) {
                            foreach $w (0..$n) {
                              if (&choose($l1,$u)*&choose($l2,($l - $u)) == 0
                                  or &choose($m1,$v)*&choose($m2,($m - $v)) == 0
                                  or &choose($n1,$w)*&choose($n2,($n - $w)) == 0) {
                                $string .= "(0)";
                              }
                              else {
                                $num1 = $l1 - $u;
                                $num2 = $l2 - $l + $u;
                                $num3 = $l - 2*$r - 2*$i;
                                $num4 = $r + $i;
                                $num5 = $m1 - $v;
                                $num6 = $m2 - $m + $v;
                                $num7 = $m - 2*$s - 2*$j;
                                $num8 = $s + $j;
                                $num9 = $n1 - $w;
                                $num10 = $n2 - $n + $w;
                                $num11 = $n - 2*$t - 2*$k;
                                $num12 = $t + $k;
                                $N = $l + $m + $n - 2*($r + $s + $t) - ($i + $j + $k);
                                $string .= ((-1)**($l+$m+$n+$i+$j+$k))*&choose($l1,$u)
                                  *&choose($l2,($l - $u))*&choose($m1,$v)*&choose($m2,($m - $v))
                                  *&choose($n1,$w)*&choose($n2,($n - $w))*&factorial($l)
                                  *&factorial($m)*&factorial($n)/(&factorial($i)*&factorial($j)
                                  *&factorial($k)*&factorial($r)*&factorial($s)*&factorial($t)
                                  *&factorial($l-2*$r-2*$i)*&factorial($m-2*$s-2*$j)
                                  *&factorial($n-2*$t-2*$k));
                                $string .= "*(PA(1)**$num1)*(PB(1)**$num2)*(PC(1)**$num3)";
                                $string .= "*(PA(2)**$num5)*(PB(2)**$num6)*(PC(2)**$num7)";
                                $string .= "*(PA(3)**$num9)*(PB(3)**$num10)*(PC(3)**$num11)";
                                $string .= "*(eps**";
                                $string .= $num4 + $num8 + $num12;
                                $string .= ")";
                                $string .= "*F$N";
                              }
                              if ($l != ($l1 + $l2) 
                                  or $r != (floor($l/2)) 
                                  or $i != (floor(($l - 2*$r)/2)) 
                                  or $u != $l 
                                  or $m != ($m1 + $m2)
                                  or $s != (floor($m/2))
                                  or $j != (floor(($m - 2*$s)/2))
                                  or $v != $m
                                  or $n != ($n1 + $n2)
                                  or $t != (floor($n/2))
                                  or $k != (floor(($n - 2*$t)/2))
                                  or $w != $n) {
                                    $string .= " + ";
                                }
                              else {
                                $string .= "";
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      $string .= ")"; # Outer enclosing paren
            
      #**************************************************************
      # Start RegEXP portion here
      #**************************************************************

      # Get rid of terms that are 1
      if ($string =~ m/\(PA\(1\)\*\*0\)/) {
        $string =~ s/\(PA\(1\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(1\)\*\*0\)/) {
        $string =~ s/\(PB\(1\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PC\(1\)\*\*0\)/) {
        $string =~ s/\(PC\(1\)\*\*0\)/1/g;
      }

      if ($string =~ m/\(PA\(2\)\*\*0\)/) {
        $string =~ s/\(PA\(2\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(2\)\*\*0\)/) {
        $string =~ s/\(PB\(2\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PC\(2\)\*\*0\)/) {
        $string =~ s/\(PC\(2\)\*\*0\)/1/g;
      }

      if ($string =~ m/\(PA\(3\)\*\*0\)/) {
        $string =~ s/\(PA\(3\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PB\(3\)\*\*0\)/) {
        $string =~ s/\(PB\(3\)\*\*0\)/1/g;
      }
      if ($string =~ m/\(PC\(3\)\*\*0\)/) {
        $string =~ s/\(PC\(3\)\*\*0\)/1/g;
      }

      if ($string =~ m/\(eps\*\*0\)/) {
        $string =~ s/\(eps\*\*0\)/1/g;
      }
      
      # Simplify terms to the 1 power
      if ($string =~ m/\(PA\(1\)\*\*1\)/) {
        $string =~ s/\(PA\(1\)\*\*1\)/PA\(1\)/g;
      }
      if ($string =~ m/\(PB\(1\)\*\*1\)/) {
        $string =~ s/\(PB\(1\)\*\*1\)/PB\(1\)/g;
      }
      if ($string =~ m/\(PC\(1\)\*\*1\)/) {
        $string =~ s/\(PC\(1\)\*\*1\)/PC\(1\)/g;
      }

      if ($string =~ m/\(PA\(2\)\*\*1\)/) {
        $string =~ s/\(PA\(2\)\*\*1\)/PA\(2\)/g;
      }
      if ($string =~ m/\(PB\(2\)\*\*1\)/) {
        $string =~ s/\(PB\(2\)\*\*1\)/PB\(2\)/g;
      }
      if ($string =~ m/\(PC\(2\)\*\*1\)/) {
        $string =~ s/\(PC\(2\)\*\*1\)/PC\(2\)/g;
      }

      if ($string =~ m/\(PA\(3\)\*\*1\)/) {
        $string =~ s/\(PA\(3\)\*\*1\)/PA\(3\)/g;
      }
      if ($string =~ m/\(PB\(3\)\*\*1\)/) {
        $string =~ s/\(PB\(3\)\*\*1\)/PB\(3\)/g;
      }
      if ($string =~ m/\(PC\(3\)\*\*1\)/) {
        $string =~ s/\(PC\(3\)\*\*1\)/PC\(3\)/g;
      }

      if ($string =~ m/\(eps\*\*1\)/) {
        $string =~ s/\(eps\*\*1\)/eps/g;
      }

      # Get of 1 multiplied by itself repeatedly
      while ($string =~ m/\*1\*/) {
        $string =~s/\*1\*/\*/g;
      }

      while ($string =~ m/ \+ \(0\)/) {
        $string =~ s/ \+ \(0\)//g;
      }

      # wo(1,1)
      if ($string =~ m/preFactor\*\(1\*F0\)/) {
        $string =~ s/preFactor\*\(1\*F0\)/preFactor\*F0/g;
      }

      if ($string =~ m/\+ -/) {
        $string =~ s/\+ -/- /g;
      }

      if ($string =~ m/\(1\*P/) {
        $string =~ s/\(1\*P/\(P/g;
      }
      if ($string =~ m/ 1\*P/) {
        $string =~ s/ 1\*P/ P/g;
      }
      if ($string =~ m/ 1\*\(/) {
        $string =~ s/ 1\*\(/ (/g;
      }
      if ($string =~ m/\(1\*\(/) {
        $string =~ s/\(1\*\(/\(\(/g;
      }
      # eps = e1
      if ($string =~ m/\*eps\*/) {
        $string =~ s/\*eps\*/\*e1\*/g;
      }
      # (eps**2) = e2
      if ($string =~ m/\(eps\*\*2\)/) {
        $string =~ s/\(eps\*\*2\)/e2/g;
      }
      # (eps**3) = e3
      if ($string =~ m/\(eps\*\*3\)/) {
        $string =~ s/\(eps\*\*3\)/e3/g;
      }
      #**************************************************************
      # End of RegEXP portion here
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

# This subroutine computes F_N(T) - which is the Incomplete
# Gamma Function, aka the Boys Function. This comes about as a 
# result of the integral transforms that are required to 
# account for the (1/r) term that is present in the Coulomb 
# integral, to make the integral separable in the Cartesian 
# directions (x,y,z).  The formalism is Eq. (29) on page 41 of the 
# following document:
#
# ****************************************************************
# Petersson T., Hellsing B., "A Detailed Derivation of Gaussian
# Orbital-Based Matrix Elements In Electron Structure Calculations."
# European Journal of Physics, Vol. 31, No. 1, Page 37.  
# ****************************************************************
#
# The Boys Function itself is defined as:
# F_N(T) = integral [from 0 to 1] {(t**(2*N))*exp(-T*(t**2))}dt
#
# where T = (a1 + a2 + a3)*sum(P - C)**2
#
# For our purposes here, N will be an integer >=0 that is a linear
# combination of index values that come from the summations that
# compute the main body of the Coulomb integral, be it Nuclear
# Attraction, or the more complicated Electron Repulsion Integral.
#
# Generally, the solution to this integral is expressed numerically
# in terms of standard error functions (erf), but this closed form
# expression becomes numerically unstable if T = 0 or is very small 
# (~10E-2 or less).  To avoid the large roundoff errors that occur
# as a result of this, a check has to be performed on the value of T.  
# If T <= 10E-2, the following formalism shown by Cook (referenced below) 
# needs to be applied:
#
# F_N(T) = (1/2)*exp(-T)*sum[from i=0 to infinity]{gamma(N+1/2)/gamma(N+i+3/2)}
#
# This infinite series is truncated at i = 6 for our desired accuracy,
# and is implemented in a separate subroutine.  This method also 
# accuratey accounts for the case where T = 0.
#
# More detailed in formation on the Boys Function and numerical 
# approximations of it are detailed in the following sources:
#
# ****************************************************************
# Cook, David B., "Handbook of Computational Quantum Chemistry."
# Dover Publications (August 22, 2005).
# ISBN-10: 0486443078
# ISBN-13: 978-0486443072
# ****************************************************************
#
# as well as:
#
# ****************************************************************
# Helgaker T., Taylor P., "Gaussian Basis Sets and Molecular 
# Integrals", Ch. 12 from "Modern Electronic Structure Theory."
# Pp. 809 - 815.
# ****************************************************************
sub F {
  # Define passed parameters
  my $N = $_[0];

  # Define local variables
  my $string;
  my $N_5;
  my $k;
  my $num;
  my $N4;

  $string = ""; # Initialize the string
  $string .= "F$N = ";
  $string .= &factorial(2*$N)/(2*&factorial($N));
  $string .= "*("; # Outer enclosing paren
  $N_5 = $N + 0.5."d0";
  $N4 = 4**$N;
  $string .= "(sqrt_pi/($N4*XX**($N_5)))*erf_XX";
  $string .= "-";
  $string .= "exp_XX";
  $string .= "*("; # Enclosing paren        
  # Here we can only proceed if N >= 0, otherwise
  # the sum is undefined.
  if (($N - 1) >= 0) {
    foreach $k (0..($N - 1)) {
      $string .= "(";
      $string .= &factorial($N - $k)/((4**$k)*&factorial(2*$N - 2*$k));
      $string .= ")";
      $string .= "*XX**(";
      $string .= -($k + 1);
      $string .= ")";
      if ($k != ($N - 1)) {
        $string .= "+";
      }
      else {
        $string .= "";
      }
    }
  }
  else {
    $string .= "0d0";
  }
  $string .= ")"; # Enclosing paren
  $string .= ")"; # Outer enclosing paren
  return $string;
}

# This subroutine computes Cook's formalism for approximating
# the Boys Function in the case of small values (T <= 10E-2)
# as detailed in the documentation for the previous subroutine.
sub altGamma {
  # Define local variables
  my $i;
  my $string;
  my $approxCodeFile = "approx";

  open (APP,">$approxCodeFile") || die "Error opening file.\n";
  $string = ""; # Initialize string to empty string
  $string .= "S(N+1) = 0.5d0*exp_XX*("; # With enclosing paren
  foreach $i (0..6) {
    $string .= "(gamma(N + 0.5d0)*XX**$i)/gamma(N + $i + 1.5d0)";
    if ($i != 6) {
      $string .= " + ";
    }
    else {
      $string .= "";
    }
  }
  $string .= ")";
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
    print APP $portion;
    $string = ""; # Reset the string before starting the loops again to prevent
                  # the solutions from concatonating onto each other.
}

# This subroutine creates a file of solved Boys Functions using Cook's 
# formalism.  These make up the body of the subroutine that defines 
# F0 through F6 for T > 10E-2.
sub incompleteGamma {
  # Declare local variables
  my $N;
  my $string;
  my $gammaCodeFile = "gamma";

  # Open the output file
  open (GAM,">$gammaCodeFile") || die "Error opening file.\n";
  $string = ""; # Initialize string to empty string
  foreach $N (0..10) {
    $string .= &F($N);
    # Here we get rid of some redundant
    # terms and operations.
    if ($string =~ m/\(1\)\*/) {
      $string =~ s/\(1\)\*//g;
    }
    if ($string =~ m/\*\*\(1\)/) {
      $string =~ s/\*\*\(1\)//g;
    }
    if ($string =~ m/-exp_XX\*\(0d0\)/) {
      $string =~ s/-exp_XX\*\(0d0\)//g;
    }

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
    print GAM $portion;
    $string = ""; # Reset the string before starting the loops again to prevent
                  # the solutions from concatonating onto each other.
  }
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

  $g[1][1] =  "g(1,1) = +wo(1,1)";
  $g[1][2] =  "g(1,2) = +wo(1,2)";
  $g[1][3] =  "g(1,3) = +wo(1,3)";
  $g[1][4] =  "g(1,4) = +wo(1,4)";
  $g[1][5] =  "g(1,5) = +wo(1,8)";
  $g[1][6] =  "g(1,6) = +wo(1,9)";
  $g[1][7] =  "g(1,7) = +wo(1,10)";
  $g[1][8] =  "g(1,8) = +wo(1,5)-wo(1,6)";
  $g[1][9] =  "g(1,9) = +2*wo(1,7) - wo(1,5)-wo(1,6)";
  $g[1][10] = "g(1,10) = +wo(1,11)";
  $g[1][11] = "g(1,11) = +wo(1,13) - wo(1,15)";
  $g[1][12] = "g(1,12) = +wo(1,18) - 3*wo(1,14)";
  $g[1][13] = "g(1,13) = +3*wo(1,12) - wo(1,19)";
  $g[1][14] = "g(1,14) = +2*wo(1,20) - 3*wo(1,13)-3*wo(1,15)";
  $g[1][15] = "g(1,15) = +4*wo(1,16) - wo(1,18)-wo(1,14)";
  $g[1][16] = "g(1,16) = +4*wo(1,17) - wo(1,12)-wo(1,19)";
  $g[2][1] =  "g(2,1) = +wo(2,1)";
  $g[2][2] =  "g(2,2) = +wo(2,2)";
  $g[2][3] =  "g(2,3) = +wo(2,3)";
  $g[2][4] =  "g(2,4) = +wo(2,4)";
  $g[2][5] =  "g(2,5) = +wo(2,8)";
  $g[2][6] =  "g(2,6) = +wo(2,9)";
  $g[2][7] =  "g(2,7) = +wo(2,10)";
  $g[2][8] =  "g(2,8) = +wo(2,5) - wo(2,6)";
  $g[2][9] =  "g(2,9) = +2*wo(2,7) - wo(2,5)-wo(2,6)";
  $g[2][10] = "g(2,10) = +wo(2,11)";
  $g[2][11] = "g(2,11) = +wo(2,13) - wo(2,15)";
  $g[2][12] = "g(2,12) = +wo(2,18) - 3*wo(2,14)";
  $g[2][13] = "g(2,13) = +3*wo(2,12) - wo(2,19)";
  $g[2][14] = "g(2,14) = +2*wo(2,20) - 3*wo(2,13) - 3*wo(2,15)";
  $g[2][15] = "g(2,15) = +4*wo(2,16) - wo(2,18) - wo(2,14)";
  $g[2][16] = "g(2,16) = +4*wo(2,17) - wo(2,12) - wo(2,19)";
  $g[3][1] =  "g(3,1) = +wo(3,1)";
  $g[3][2] =  "g(3,2) = +wo(3,2)";
  $g[3][3] =  "g(3,3) = +wo(3,3)";
  $g[3][4] =  "g(3,4) = +wo(3,4)";
  $g[3][5] =  "g(3,5) = +wo(3,8)";
  $g[3][6] =  "g(3,6) = +wo(3,9)";
  $g[3][7] =  "g(3,7) = +wo(3,10)";
  $g[3][8] =  "g(3,8) = +wo(3,5) - wo(3,6)";
  $g[3][9] =  "g(3,9) = +2*wo(3,7) - wo(3,5)-wo(3,6)";
  $g[3][10] = "g(3,10) = +wo(3,11)";
  $g[3][11] = "g(3,11) = +wo(3,13) - wo(3,15)";
  $g[3][12] = "g(3,12) = +wo(3,18) - 3*wo(3,14)";
  $g[3][13] = "g(3,13) = +3*wo(3,12) - wo(3,19)";
  $g[3][14] = "g(3,14) = +2*wo(3,20) - 3*wo(3,13) - 3*wo(3,15)";
  $g[3][15] = "g(3,15) = +4*wo(3,16) - wo(3,18) - wo(3,14)";
  $g[3][16] = "g(3,16) = +4*wo(3,17) - wo(3,12) - wo(3,19)";
  $g[4][1] =  "g(4,1) = +wo(4,1)";
  $g[4][2] =  "g(4,2) = +wo(4,2)";
  $g[4][3] =  "g(4,3) = +wo(4,3)";
  $g[4][4] =  "g(4,4) = +wo(4,4)";
  $g[4][5] =  "g(4,5) = +wo(4,8)";
  $g[4][6] =  "g(4,6) = +wo(4,9)";
  $g[4][7] =  "g(4,7) = +wo(4,10)";
  $g[4][8] =  "g(4,8) = +wo(4,5) - wo(4,6)";
  $g[4][9] =  "g(4,9) = +2*wo(4,7) - wo(4,5) - wo(4,6)";
  $g[4][10] = "g(4,10) = +wo(4,11)";
  $g[4][11] = "g(4,11) = +wo(4,13) - wo(4,15)";
  $g[4][12] = "g(4,12) = +wo(4,18) - 3*wo(4,14)";
  $g[4][13] = "g(4,13) = +3*wo(4,12) - wo(4,19)";
  $g[4][14] = "g(4,14) = +2*wo(4,20) - 3*wo(4,13) - 3*wo(4,15)";
  $g[4][15] = "g(4,15) = +4*wo(4,16) - wo(4,18) - wo(4,14)";
  $g[4][16] = "g(4,16) = +4*wo(4,17) - wo(4,12) - wo(4,19)";
  $g[5][1] =  "g(5,1) = +wo(8,1)";
  $g[5][2] =  "g(5,2) = +wo(8,2)";
  $g[5][3] =  "g(5,3) = +wo(8,3)";
  $g[5][4] =  "g(5,4) = +wo(8,4)";
  $g[5][5] =  "g(5,5) = +wo(8,8)";
  $g[5][6] =  "g(5,6) = +wo(8,9)";
  $g[5][7] =  "g(5,7) = +wo(8,10)";
  $g[5][8] =  "g(5,8) = +wo(8,5) - wo(8,6)";
  $g[5][9] =  "g(5,9) = +2*wo(8,7) - wo(8,5) - wo(8,6)";
  $g[5][10] = "g(5,10) = +wo(8,11)";
  $g[5][11] = "g(5,11) = +wo(8,13) - wo(8,15)";
  $g[5][12] = "g(5,12) = +wo(8,18) - 3*wo(8,14)";
  $g[5][13] = "g(5,13) = +3*wo(8,12) - wo(8,19)";
  $g[5][14] = "g(5,14) = +2*wo(8,20) - 3*wo(8,13) - 3*wo(8,15)";
  $g[5][15] = "g(5,15) = +4*wo(8,16) - wo(8,18) - wo(8,14)";
  $g[5][16] = "g(5,16) = +4*wo(8,17) - wo(8,12) - wo(8,19)";
  $g[6][1] =  "g(6,1) = +wo(9,1)";
  $g[6][2] =  "g(6,2) = +wo(9,2)";
  $g[6][3] =  "g(6,3) = +wo(9,3)";
  $g[6][4] =  "g(6,4) = +wo(9,4)";
  $g[6][5] =  "g(6,5) = +wo(9,8)";
  $g[6][6] =  "g(6,6) = +wo(9,9)";
  $g[6][7] =  "g(6,7) = +wo(9,10)";
  $g[6][8] =  "g(6,8) = +wo(9,5) - wo(9,6)";
  $g[6][9] =  "g(6,9) = +2*wo(9,7) - wo(9,5) - wo(9,6)";
  $g[6][10] = "g(6,10) = +wo(9,11)";
  $g[6][11] = "g(6,11) = +wo(9,13) - wo(9,15)";
  $g[6][12] = "g(6,12) = +wo(9,18) - 3*wo(9,14)";
  $g[6][13] = "g(6,13) = +3*wo(9,12) - wo(9,19)";
  $g[6][14] = "g(6,14) = +2*wo(9,20) - 3*wo(9,13) - 3*wo(9,15)";
  $g[6][15] = "g(6,15) = +4*wo(9,16) - wo(9,18) - wo(9,14)";
  $g[6][16] = "g(6,16) = +4*wo(9,17) - wo(9,12) - wo(9,19)";
  $g[7][1] =  "g(7,1) = +wo(10,1)";
  $g[7][2] =  "g(7,2) = +wo(10,2)";
  $g[7][3] =  "g(7,3) = +wo(10,3)";
  $g[7][4] =  "g(7,4) = +wo(10,4)";
  $g[7][5] =  "g(7,5) = +wo(10,8)";
  $g[7][6] =  "g(7,6) = +wo(10,9)";
  $g[7][7] =  "g(7,7) = +wo(10,10)";
  $g[7][8] =  "g(7,8) = +wo(10,5) - wo(10,6)";
  $g[7][9] =  "g(7,9) = +2*wo(10,7) - wo(10,5) - wo(10,6)";
  $g[7][10] = "g(7,10) = +wo(10,11)";
  $g[7][11] = "g(7,11) = +wo(10,13) - wo(10,15)";
  $g[7][12] = "g(7,12) = +wo(10,18) - 3*wo(10,14)";
  $g[7][13] = "g(7,13) = +3*wo(10,12) - wo(10,19)";
  $g[7][14] = "g(7,14) = +2*wo(10,20) - 3*wo(10,13) - 3*wo(10,15)";
  $g[7][15] = "g(7,15) = +4*wo(10,16) - wo(10,18) - wo(10,14)";
  $g[7][16] = "g(7,16) = +4*wo(10,17) - wo(10,12) - wo(10,19)";
  $g[8][1] =  "g(8,1) = +wo(5,1)-wo(6,1)";
  $g[8][2] =  "g(8,2) = +wo(5,2)-wo(6,2)";
  $g[8][3] =  "g(8,3) = +wo(5,3)-wo(6,3)";
  $g[8][4] =  "g(8,4) = +wo(5,4)-wo(6,4)";
  $g[8][5] =  "g(8,5) = +wo(5,8)-wo(6,8)";
  $g[8][6] =  "g(8,6) = +wo(5,9)-wo(6,9)";
  $g[8][7] =  "g(8,7) = +wo(5,10) - wo(6,10)";
  $g[8][8] =  "g(8,8) = +wo(5,5)-wo(5,6) - wo(6,5) + wo(6,6)";
  $g[8][9] = "g(8,9) = +2*wo(5,7) - wo(5,5) - wo(5,6)&\n";
  $g[8][9] .= "&- 2*wo(6,7) + wo(6,5) + wo(6,6)\n";
  $g[8][10] = "g(8,10) = +wo(5,11) - wo(6,11)";
  $g[8][11] = "g(8,11) = +wo(5,13) - wo(5,15) - wo(6,13) +wo(6,15)";
  $g[8][12] = "g(8,12) = +wo(5,18) - 3*wo(5,14) - wo(6,18) + 3*wo(6,14)";
  $g[8][13] = "g(8,13) = +3*wo(5,12) - wo(5,19) - 3*wo(6,12) + wo(6,19)";
  $g[8][14] = "g(8,14) = +2*wo(5,20) - 3*wo(5,13) - 3*wo(5,15)&\n";
  $g[8][14] .= "&- 2*wo(6,20) + 3*wo(6,13) +3*wo(6,15)\n";
  $g[8][15] = "g(8,15) = +4*wo(5,16) - wo(5,18) - wo(5,14)&\n";
  $g[8][15] .= "&- 4*wo(6,16) + wo(6,18) + wo(6,14)\n";
  $g[8][16] = "g(8,16) = +4*wo(5,17) - wo(5,12) - wo(5,19)&\n";
  $g[8][16] .= "&- 4*wo(6,17) + wo(6,12) +wo(6,19)\n";
  $g[9][1] = "g(9,1) = +2*wo(7,1) - wo(5,1) - wo(6,1)";
  $g[9][2] = "g(9,2) = +2*wo(7,2) - wo(5,2) - wo(6,2)";
  $g[9][3] = "g(9,3) = +2*wo(7,3) - wo(5,3) - wo(6,3)";
  $g[9][4] = "g(9,4) = +2*wo(7,4) - wo(5,4) - wo(6,4)";
  $g[9][5] = "g(9,5) = +2*wo(7,8) - wo(5,8) - wo(6,8)";
  $g[9][6] = "g(9,6) = +2*wo(7,9) - wo(5,9) - wo(6,9)";
  $g[9][7] = "g(9,7) = +2*wo(7,10) - wo(5,10) - wo(6,10)";
  $g[9][8] = "g(9,8) = +2*wo(7,5) - 2*wo(7,6) - wo(5,5)&\n";
  $g[9][8] .= "&+ wo(5,6) - wo(6,5) + wo(6,6)\n";
  $g[9][9] = "g(9,9) = +4*wo(7,7) - 2*wo(7,5) - 2*wo(7,6)&\n";
  $g[9][9] .= "&- 2*wo(5,7) + wo(5,5) + wo(5,6) - 2*wo(6,7) + wo(6,5) + wo(6,6)\n";
  $g[9][10] = "g(9,10) = +2*wo(7,11) - wo(5,11) - wo(6,11)";
  $g[9][11] = "g(9,11) = +2*wo(7,13) - 2*wo(7,15) - wo(5,13)&\n";
  $g[9][11] .= "&+ wo(5,15) - wo(6,13) + wo(6,15)\n";
  $g[9][12] = "g(9,12) = +2*wo(7,18) - 6*wo(7,14) - wo(5,18)&\n";
  $g[9][12] .= "&+ 3*wo(5,14) - wo(6,18) + 3*wo(6,14)\n";
  $g[9][13] = "g(9,13) = +6*wo(7,12) - 2*wo(7,19) - 3*wo(5,12)&\n";
  $g[9][13] .= "&+ wo(5,19) - 3*wo(6,12) + wo(6,19)\n";
  $g[9][14] = "g(9,14) = +4*wo(7,20) - 6*wo(7,13) - 6*wo(7,15)&\n";
  $g[9][14] .= "&- 2*wo(5,20) + 3*wo(5,13) + 3*wo(5,15) - 2*wo(6,20)&\n";
  $g[9][14] .= "&+ 3*wo(6,13) + 3*wo(6,15)\n";
  $g[9][15] = "g(9,15) = +8*wo(7,16) - 2*wo(7,18) - 2*wo(7,14)&\n";
  $g[9][15] .= "&- 4*wo(5,16) + wo(5,18) + wo(5,14) - 4*wo(6,16)&\n";
  $g[9][15] .= "&+ wo(6,18) + wo(6,14)\n";
  $g[9][16] = "g(9,16) = +8*wo(7,17) - 2*wo(7,12) - 2*wo(7,19)&\n";
  $g[9][16] .= "&- 4*wo(5,17) + wo(5,12) + wo(5,19) - 4*wo(6,17)&\n";
  $g[9][16] .= "&+ wo(6,12) + wo(6,19)\n";
  $g[10][1] =  "g(10,1) = +wo(11,1)";
  $g[10][2] =  "g(10,2) = +wo(11,2)";
  $g[10][3] =  "g(10,3) = +wo(11,3)";
  $g[10][4] =  "g(10,4) = +wo(11,4)";
  $g[10][5] =  "g(10,5) = +wo(11,8)";
  $g[10][6] =  "g(10,6) = +wo(11,9)";
  $g[10][7] =  "g(10,7) = +wo(11,10)";
  $g[10][8] =  "g(10,8) = +wo(11,5) - wo(11,6)";
  $g[10][9] =  "g(10,9) = +2*wo(11,7) - wo(11,5) - wo(11,6)";
  $g[10][10] = "g(10,10) = +wo(11,11)";
  $g[10][11] = "g(10,11) = +wo(11,13) - wo(11,15)";
  $g[10][12] = "g(10,12) = +wo(11,18) - 3*wo(11,14)";
  $g[10][13] = "g(10,13) = +3*wo(11,12) - wo(11,19)";
  $g[10][14] = "g(10,14) = +2*wo(11,20) - 3*wo(11,13) - 3*wo(11,15)";
  $g[10][15] = "g(10,15) = +4*wo(11,16) - wo(11,18) - wo(11,14)";
  $g[10][16] = "g(10,16) = +4*wo(11,17) - wo(11,12) - wo(11,19)";
  $g[11][1] =  "g(11,1) = +wo(13,1) - wo(15,1)";
  $g[11][2] =  "g(11,2) = +wo(13,2) - wo(15,2)";
  $g[11][3] =  "g(11,3) = +wo(13,3) - wo(15,3)";
  $g[11][4] =  "g(11,4) = +wo(13,4) - wo(15,4)";
  $g[11][5] =  "g(11,5) = +wo(13,8) - wo(15,8)";
  $g[11][6] =  "g(11,6) = +wo(13,9) - wo(15,9)";
  $g[11][7] =  "g(11,7) = +wo(13,10) - wo(15,10)";
  $g[11][8] =  "g(11,8) = +wo(13,5) - wo(13,6)-wo(15,5) + wo(15,6)";
  $g[11][9] =  "g(11,9) = +2*wo(13,7) - wo(13,5) - wo(13,6)&\n";
  $g[11][9] .= "&- 2*wo(15,7) + wo(15,5) +wo(15,6)\n";
  $g[11][10] = "g(11,10) = +wo(13,11) - wo(15,11)";
  $g[11][11] = "g(11,11) = +wo(13,13) - wo(13,15) - wo(15,13) + wo(15,15)";
  $g[11][12] = "g(11,12) = +wo(13,18) - 3*wo(13,14) - wo(15,18) + 3*wo(15,14)";
  $g[11][13] = "g(11,13) = +3*wo(13,12) - wo(13,19) - 3*wo(15,12) + wo(15,19)";
  $g[11][14] = "g(11,14) = +2*wo(13,20) - 3*wo(13,13) - 3*wo(13,15)&\n";
  $g[11][14] .= "&- 2*wo(15,20) + 3*wo(15,13) + 3*wo(15,15)\n";
  $g[11][15] = "g(11,15) = +4*wo(13,16) - wo(13,18) - wo(13,14)&\n";
  $g[11][15] .= "&- 4*wo(15,16) + wo(15,18) + wo(15,14)\n";
  $g[11][16] = "g(11,16) = +4*wo(13,17) - wo(13,12) - wo(13,19)&\n";
  $g[11][16] .= "&- 4*wo(15,17) + wo(15,12) + wo(15,19)\n";
  $g[12][1] = "g(12,1) = +wo(18,1) - 3*wo(14,1)";
  $g[12][2] = "g(12,2) = +wo(18,2) - 3*wo(14,2)";
  $g[12][3] = "g(12,3) = +wo(18,3) - 3*wo(14,3)";
  $g[12][4] = "g(12,4) = +wo(18,4) - 3*wo(14,4)";
  $g[12][5] = "g(12,5) = +wo(18,8) - 3*wo(14,8)";
  $g[12][6] = "g(12,6) = +wo(18,9) - 3*wo(14,9)";
  $g[12][7] = "g(12,7) = +wo(18,10) - 3*wo(14,10)";
  $g[12][8] = "g(12,8) = +wo(18,5) - wo(18,6) - 3*wo(14,5) + 3*wo(14,6)";
  $g[12][9] = "g(12,9) = +2*wo(18,7) - wo(18,5)&\n";
  $g[12][9] .= "&- wo(18,6) - 6*wo(14,7) + 3*wo(14,5) + 3*wo(14,6)\n";
  $g[12][10] = "g(12,10) = +wo(18,11) - 3*wo(14,11)";
  $g[12][11] = "g(12,11) = +wo(18,13) - wo(18,15) - 3*wo(14,13) + 3*wo(14,15)  ";
  $g[12][12] = "g(12,12) = +wo(18,18) - 3*wo(18,14) - 3*wo(14,18) + 9*wo(14,14)";
  $g[12][13] = "g(12,13) = +3*wo(18,12) - wo(18,19) - 9*wo(14,12) + 3*wo(14,19)";
  $g[12][14] = "g(12,14) = +2*wo(18,20) - 3*wo(18,13) - 3*wo(18,15)&\n";
  $g[12][14] .= "&- 6*wo(14,20) + 9*wo(14,13) + 9*wo(14,15)\n";
  $g[12][15] = "g(12,15) = +4*wo(18,16) - wo(18,18) - wo(18,14)&\n";
  $g[12][15] .= "&- 12*wo(14,16) + 3*wo(14,18) + 3*wo(14,14)\n";
  $g[12][16] = "g(12,16) = +4*wo(18,17) - wo(18,12) - wo(18,19)&\n";
  $g[12][16] .= "&- 12*wo(14,17) + 3*wo(14,12) + 3*wo(14,19)\n";
  $g[13][1] = "g(13,1) = +3*wo(12,1) - wo(19,1)";
  $g[13][2] = "g(13,2) = +3*wo(12,2) - wo(19,2)";
  $g[13][3] = "g(13,3) = +3*wo(12,3) - wo(19,3)";
  $g[13][4] = "g(13,4) = +3*wo(12,4) - wo(19,4)";
  $g[13][5] = "g(13,5) = +3*wo(12,8) - wo(19,8)";
  $g[13][6] = "g(13,6) = +3*wo(12,9) - wo(19,9)";
  $g[13][7] = "g(13,7) = +3*wo(12,10) - wo(19,10)";
  $g[13][8] = "g(13,8) = +3*wo(12,5) - 3*wo(12,6) - wo(19,5) +wo(19,6)";
  $g[13][9] = "g(13,9) = +6*wo(12,7) - 3*wo(12,5)&\n";
  $g[13][9] .= "&- 3*wo(12,6) - 2*wo(19,7) + wo(19,5) +wo(19,6)\n";
  $g[13][10] = "g(13,10) = +3*wo(12,11) - wo(19,11)";
  $g[13][11] = "g(13,11) = +3*wo(12,13) - 3*wo(12,15) - wo(19,13) + wo(19,15)  ";
  $g[13][12] = "g(13,12) = +3*wo(12,18) - 9*wo(12,14) - wo(19,18) + 3*wo(19,14)";
  $g[13][13] = "g(13,13) = +9*wo(12,12) - 3*wo(12,19) - 3*wo(19,12) + wo(19,19)";
  $g[13][14] = "g(13,14) = +6*wo(12,20) - 9*wo(12,13) - 9*wo(12,15)&\n";
  $g[13][14] .= "&- 2*wo(19,20) + 3*wo(19,13) + 3*wo(19,15)\n";
  $g[13][15] = "g(13,15) = +12*wo(12,16) - 3*wo(12,18) - 3*wo(12,14)&\n";
  $g[13][15] .= "&- 4*wo(19,16) + wo(19,18) + wo(19,14)\n";
  $g[13][16] = "g(13,16) = +12*wo(12,17) - 3*wo(12,12) - 3*wo(12,19)&\n";
  $g[13][16] .= "&- 4*wo(19,17) + wo(19,12) + wo(19,19)\n";
  $g[14][1] = "g(14,1) = +2*wo(20,1) - 3*wo(13,1) - 3*wo(15,1)";
  $g[14][2] = "g(14,2) = +2*wo(20,2) - 3*wo(13,2) - 3*wo(15,2)";
  $g[14][3] = "g(14,3) = +2*wo(20,3) - 3*wo(13,3) - 3*wo(15,3)";
  $g[14][4] = "g(14,4) = +2*wo(20,4) - 3*wo(13,4) - 3*wo(15,4)";
  $g[14][5] = "g(14,5) = +2*wo(20,8) - 3*wo(13,8) - 3*wo(15,8)";
  $g[14][6] = "g(14,6) = +2*wo(20,9) - 3*wo(13,9) - 3*wo(15,9)";
  $g[14][7] = "g(14,7) = +2*wo(20,10) - 3*wo(13,10) - 3*wo(15,10)";
  $g[14][8] = "g(14,8) = +2*wo(20,5) - 2*wo(20,6) - 3*wo(13,5)&\n";
  $g[14][8] .= "&+ 3*wo(13,6) - 3*wo(15,5) +3*wo(15,6)\n";
  $g[14][9] = "g(14,9) = +4*wo(20,7) - 2*wo(20,5) - 2*wo(20,6)&\n";
  $g[14][9] .= "&- 6*wo(13,7) + 3*wo(13,5) + 3*wo(13,6)&\n";
  $g[14][9] .= "&- 6*wo(15,7) + 3*wo(15,5) + 3*wo(15,6)\n";
  $g[14][10] = "g(14,10) = +2*wo(20,11) - 3*wo(13,11) - 3*wo(15,11)";
  $g[14][11] = "g(14,11) = +2*wo(20,13) - 2*wo(20,15) - 3*wo(13,13)&\n";
  $g[14][11] .= "&+ 3*wo(13,15) - 3*wo(15,13) + 3*wo(15,15)\n";
  $g[14][12] = "g(14,12) = +2*wo(20,18) - 6*wo(20,14) - 3*wo(13,18)&\n";
  $g[14][12] .= "&+ 9*wo(13,14) - 3*wo(15,18) + 9*wo(15,14)\n";
  $g[14][13] = "g(14,13) = +6*wo(20,12) - 2*wo(20,19) - 9*wo(13,12)&\n";
  $g[14][13] .= "&+ 3*wo(13,19) - 9*wo(15,12) + 3*wo(15,19)\n";
  $g[14][14] = "g(14,14) = +4*wo(20,20) - 6*wo(20,13) - 6*wo(20,15)&\n";
  $g[14][14] .= "&- 6*wo(13,20) + 9*wo(13,13) + 9*wo(13,15)&\n";
  $g[14][14] .= "&- 6*wo(15,20) + 9*wo(15,13) + 9*wo(15,15)\n";
  $g[14][15] = "g(14,15) = +8*wo(20,16) - 2*wo(20,18) - 2*wo(20,14)&\n";
  $g[14][15] .= "&- 12*wo(13,16) + 3*wo(13,18) + 3*wo(13,14)&\n";
  $g[14][15] .= "&- 12*wo(15,16) + 3*wo(15,18) + 3*wo(15,14)\n";
  $g[14][16] = "g(14,16) = +8*wo(20,17) - 2*wo(20,12) - 2*wo(20,19)&\n";
  $g[14][16] .= "&- 12*wo(13,17) + 3*wo(13,12) + 3*wo(13,19)&\n";
  $g[14][16] .= "&- 12*wo(15,17) + 3*wo(15,12) + 3*wo(15,19)\n";
  $g[15][1] = "g(15,1) = +4*wo(16,1) - wo(18,1) - wo(14,1)";
  $g[15][2] = "g(15,2) = +4*wo(16,2) - wo(18,2) - wo(14,2)";
  $g[15][3] = "g(15,3) = +4*wo(16,3) - wo(18,3) - wo(14,3)";
  $g[15][4] = "g(15,4) = +4*wo(16,4) - wo(18,4) - wo(14,4)";
  $g[15][5] = "g(15,5) = +4*wo(16,8) - wo(18,8) - wo(14,8)";
  $g[15][6] = "g(15,6) = +4*wo(16,9) - wo(18,9) - wo(14,9)";
  $g[15][7] = "g(15,7) = +4*wo(16,10) - wo(18,10) - wo(14,10)";
  $g[15][8] = "g(15,8) = +4*wo(16,5) - 4*wo(16,6) - wo(18,5)&\n";
  $g[15][8] .= "&+ wo(18,6) - wo(14,5) + wo(14,6)\n";
  $g[15][9] = "g(15,9) = +8*wo(16,7) - 4*wo(16,5) - 4*wo(16,6)&\n";
  $g[15][9] .= "&- 2*wo(18,7) + wo(18,5) + wo(18,6)&\n";
  $g[15][9] .= "&- 2*wo(14,7) + wo(14,5)+wo(14,6)\n";
  $g[15][10] = "g(15,10) = +4*wo(16,11) - wo(18,11)-wo(14,11)";
  $g[15][11] = "g(15,11) = +4*wo(16,13) - 4*wo(16,15) - wo(18,13)&\n";
  $g[15][11] .= "&+ wo(18,15) - wo(14,13) + wo(14,15)\n";
  $g[15][12] = "g(15,12) = +4*wo(16,18) - 12*wo(16,14) - wo(18,18)&\n";
  $g[15][12] .= "&+ 3*wo(18,14) - wo(14,18) + 3*wo(14,14)\n";
  $g[15][13] = "g(15,13) = +12*wo(16,12) - 4*wo(16,19) - 3*wo(18,12)&\n";
  $g[15][13] .= "&+ wo(18,19) - 3*wo(14,12) + wo(14,19)\n";
  $g[15][14] = "g(15,14) = +8*wo(16,20) - 12*wo(16,13) - 12*wo(16,15)&\n";
  $g[15][14] .= "&- 2*wo(18,20) + 3*wo(18,13) + 3*wo(18,15)&\n";
  $g[15][14] .= "&- 2*wo(14,20) + 3*wo(14,13) + 3*wo(14,15)\n";
  $g[15][15] = "g(15,15) = +16*wo(16,16) - 4*wo(16,18) - 4*wo(16,14)&\n";
  $g[15][15] .= "&- 4*wo(18,16) + wo(18,18) + wo(18,14)&\n";
  $g[15][15] .= "&- 4*wo(14,16) + wo(14,18) + wo(14,14)";
  $g[15][16] = "g(15,16) = +16*wo(16,17) - 4*wo(16,12) - 4*wo(16,19)&\n";
  $g[15][16] .= "&- 4*wo(18,17) + wo(18,12) + wo(18,19)&\n";
  $g[15][16] .= "&- 4*wo(14,17) + wo(14,12) + wo(14,19)\n";
  $g[16][1] = "g(16,1) = +4*wo(17,1) - wo(12,1) - wo(19,1)";
  $g[16][2] = "g(16,2) = +4*wo(17,2) - wo(12,2) - wo(19,2)";
  $g[16][3] = "g(16,3) = +4*wo(17,3) - wo(12,3) - wo(19,3)";
  $g[16][4] = "g(16,4) = +4*wo(17,4) - wo(12,4) - wo(19,4)";
  $g[16][5] = "g(16,5) = +4*wo(17,8) - wo(12,8) - wo(19,8)";
  $g[16][6] = "g(16,6) = +4*wo(17,9) - wo(12,9) - wo(19,9)";
  $g[16][7] = "g(16,7) = +4*wo(17,10)  -wo(12,10) - wo(19,10)";
  $g[16][8] = "g(16,8) = +4*wo(17,5) - 4*wo(17,6) - wo(12,5)&\n";
  $g[16][8] .= "&+ wo(12,6) - wo(19,5) + wo(19,6)\n";
  $g[16][9] = "g(16,9) = +8*wo(17,7) - 4*wo(17,5) - 4*wo(17,6)&\n";
  $g[16][9] .= "&- 2*wo(12,7) + wo(12,5) + wo(12,6)&\n";
  $g[16][9] .= "&- 2*wo(19,7) + wo(19,5) + wo(19,6)\n";
  $g[16][10] = "g(16,10) = +4*wo(17,11) - wo(12,11) - wo(19,11)";
  $g[16][11] = "g(16,11) = +4*wo(17,13) - 4*wo(17,15) - wo(12,13)&\n";
  $g[16][11] .= "&+ wo(12,15) - wo(19,13) + wo(19,15)\n";
  $g[16][12] = "g(16,12) = +4*wo(17,18) - 12*wo(17,14) - wo(12,18)&\n";
  $g[16][12] .= "&+ 3*wo(12,14) - wo(19,18) + 3*wo(19,14)\n";
  $g[16][13] = "g(16,13) = +12*wo(17,12) - 4*wo(17,19) - 3*wo(12,12)&\n";
  $g[16][13] .= "&+ wo(12,19) - 3*wo(19,12) + wo(19,19)\n";
  $g[16][14] = "g(16,14) = +8*wo(17,20) - 12*wo(17,13) - 12*wo(17,15)&\n";
  $g[16][14] .= "&- 2*wo(12,20) + 3*wo(12,13) + 3*wo(12,15)&\n";
  $g[16][14] .= "&- 2*wo(19,20) + 3*wo(19,13) + 3*wo(19,15)\n";
  $g[16][15] = "g(16,15) = +16*wo(17,16) - 4*wo(17,18) - 4*wo(17,14)&\n";
  $g[16][15] .= "&- 4*wo(12,16) + wo(12,18) + wo(12,14)&\n";
  $g[16][15] .= "&- 4*wo(19,16) + wo(19,18) + wo(19,14)\n";
  $g[16][16] = "g(16,16) = +16*wo(17,17) - 4*wo(17,12) - 4*wo(17,19)&\n";
  $g[16][16] .= "&- 4*wo(12,17) + wo(12,12) + wo(12,19)&\n";
  $g[16][16] .= "&- 4*wo(19,17) + wo(19,12) + wo(19,19)\n";
  
  $string = "";
  foreach $index2 (1..$orb2) {
    foreach $index1 (1..$orb1) {
      $string .= $g[$index1][$index2];
      $string .= "\n";
    }
  }
  return $string;
}

# This subroutine prints the formatted nuclear potential integral
# formulas implementing the l1l2switch formalism.
sub printNuc {
  # Declare local variables
  my $string;
  my $nucCodeFile = "NucIntg";
  
  # Open the output file
  open (PRI,">$nucCodeFile") || die "Error opening file.\n";

  print PRI "if (l1l2switch.eq.34) then\n\n";
  $string = &nucOverlap(4,4);
  $string .= &print_g(4,4);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.17) then\n\n";
  $string = &nucOverlap(1,1);
  $string .= &print_g(1,1);
  print PRI "$string\n";
 
  print PRI "else if (l1l2switch.eq.33) then\n\n";
  $string = &nucOverlap(1,4);
  $string .= &print_g(1,4);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.65) then\n\n";
  $string = &nucOverlap(1,10);
  $string .= &print_g(1,9);
  print PRI "$string\n";
 
  print PRI "else if (l1l2switch.eq.129) then\n\n";
  $string = &nucOverlap(1,20);
  $string .= &print_g(1,16);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.18) then\n\n";
  $string = &nucOverlap(4,1);
  $string .= &print_g(4,1);
  print PRI "$string\n";
 
  print PRI "else if (l1l2switch.eq.66) then\n\n";
  $string = &nucOverlap(4,10);
  $string .= &print_g(4,9);
  print PRI "$string\n";
 
  print PRI "else if (l1l2switch.eq.130) then\n\n";
  $string = &nucOverlap(4,20);
  $string .= &print_g(4,16);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.20) then\n\n";
  $string = &nucOverlap(10,1);
  $string .= &print_g(9,1);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.36) then\n\n";
  $string = &nucOverlap(10,4);
  $string .= &print_g(9,4);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.68) then\n\n";
  $string = &nucOverlap(10,10);
  $string .= &print_g(9,9);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.132) then\n\n";
  $string = &nucOverlap(10,20);
  $string .= &print_g(9,16);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.24) then\n\n";
  $string = &nucOverlap(20,1);
  $string .= &print_g(16,1);
  print PRI "$string\n";
  
  print PRI "else if (l1l2switch.eq.40) then\n\n";
  $string = &nucOverlap(20,4);
  $string .= &print_g(16,4);
  print PRI "$string\n";
    
  print PRI "else if (l1l2switch.eq.72) then\n\n";
  $string = &nucOverlap(20,10);
  $string .= &print_g(16,9);
  print PRI "$string\n";
    
  print PRI "else\n\n";
  $string = &nucOverlap(20,20);
  $string .= &print_g(16,16);
  print PRI "$string";
  print PRI "end if";
}





























