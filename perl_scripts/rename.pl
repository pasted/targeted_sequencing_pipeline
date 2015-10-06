#!/usr/bin/perl
#rename.pl - rename filenames with Perl regexp substitution
#./rename.pl 's/(.*)_[ACGT]{6}_(.*).fastq/$1_$2.fastq/' ../raw_reads/*.fastq
($regexp = shift @ARGV) || die "Usage:  rename perlexpr [filenames]\n";
if (!@ARGV) {
    @ARGV = <STDIN>;
    chomp(@ARGV);
}
foreach $_ (@ARGV) {
    $old_name = $_;
    eval $regexp;
    die $@ if $@;
    rename($old_name, $_) unless $old_name eq $_;
}
exit(0);
