#!/usr/bin/perl
#
# This is a script that will chop a restart file up into individual
# restart files that each contain a single molecule (or atom). It
# handles restart files with only one molecule type for the time being.

# Get the command line arguements. We need the name of the
# file we wish to chop up
@resfile = @ARGV;

# Print some information
print "We are using the file $resfile[0] \n";

# Open the file for reading.
open(FILE, $resfile[0]);

# Molecule name keyword
$moleckey = "MOLECULE";

# Find the molecule key
while ($line = <FILE>) {
    chomp ($line);
    if ( $line =~ /MOLECULE/ ) {
	print "$line \n";
	# Get the molecule name
	@junk = split(/ /,$line);
	# Get everything from the first noncharacter to the end
	$name = $junk[1];
        chomp ($name);
	print "$name\n";
	# Get the number of atoms
	$line = <FILE>;
	@junk = split(/\#/,$line);
	$nmoles = $junk[0];
	$line = <FILE>;
	@junk = split(/\#/,$line);
	$natoms = $junk[0];
	# Check for the beginning of the coordinates
	$line = <FILE>;
	if ( $line =~ /Principal Coordinates/ ) {
	    $n = 1;
	    while ($n <= $nmoles) {

                # Open the new file for writing
		open(OUTF, ">$resfile[0].$n");
		print OUTF "Comment line - Nothing to Say\n";
		print OUTF "        0 \n";
		print OUTF "        0.0E+0 \n";

		print OUTF "_MOLECULE_NAME_: $name \n";
		print OUTF "        1 # Number of molecules \n";
		print OUTF "  $natoms # Number of atoms \n";
		print OUTF "_Principal Coordinates Only_\n";
		$a = 1;
		while ($a <= $natoms) {
		    $line = <FILE>;
		    print OUTF $line;
		    $a = $a + 1;
		}
		close(OUTF);
		$n = $n + 1;
	    }
	}
	#Read the next line - should be velocities
	$line = <FILE>;
	if ( $line =~ /Velocities/ ) {
	    $n = 1;
	    while ($n <= $nmoles) {
		open(OUTF, ">>$resfile[0].$n");

		print OUTF "_Velocities_\n";
		$a = 1;
		while ($a <= $natoms) {
		    $line = <FILE>;
		    print OUTF $line;
		    $a = $a + 1;
		}
		close(OUTF);
		$n = $n + 1;
	    }
	}
	# Read the next line - should be generalized coordinates
	$line = <FILE>;
	if ( $line =~ /Generalized Coordinates/ ) {
	    $n = 1;
	    while ($n <= $nmoles) {
		open(OUTF, ">>$resfile[0].$n");

		print OUTF "_Generalized Coordinates_\n";
		$a = 1;
		while ($a <= $natoms) {
		    $line = <FILE>;
		    print OUTF $line;
		    $line = <FILE>;
		    print OUTF $line;
		    $a = $a + 1;
		}
		close(OUTF);
		$n = $n + 1;
	    }
	}
	# Ok, at this point we need to loop back through the list of
	# files and add the rest of the stuff to the end
	while ($line = <FILE>) {
	    $n = 1;
	    while ($n <= $nmoles) {
		open(OUTF, ">>$resfile[0].$n");
		print OUTF $line;
		close(OUTF);
		$n = $n + 1;
	    }
	}
	# Die in peace
	die "Ok, everything should be finished!\n"
    }
}
# Close the file
close(FILE);


