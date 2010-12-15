#!/bin/perl

use strict; 
use warnings;


my $i;
my @out = ();
my @lines = ();
my $gitrev = `git log | head -n 1`;
chomp($gitrev);
$gitrev =~ s/commit\s+//;

print STDOUT "Editing copyright file\n";
open(FH, "copyright") || die;
@lines = <FH>;
close(FH);

$i=0;
@out = ();
while($i < scalar(@lines)) {
	my $line = $lines[$i];

	if($line =~ m/^### SELECT/) {
		$i += 4;
	}
	elsif($line =~ m/^#/) {
	}
	elsif($line =~ m/<likewise/) {
	}
	else {
		$line =~ s/url:\/\/example.com/http:\/\/www.iontorrent.com/;
		$line =~ s/<put author.*/Nils Homer <nils.homer\@lifetech.com>/;
		$line =~ s/<Put the license.*/GNU GPL V2/;
		$line =~ s/<Copyright \(C\).*/Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved/;
		push(@out, $line);
	}
	$i++;
}

open(FH, ">copyright") || die;
foreach my $line (@out) {
	print FH "$line";
}
close(FH);

print STDOUT "Editing changelog file\n";
open(FH, "changelog");
@lines = <FH>;
close(FH);

$i=0;
@out = ();
while($i < scalar(@lines)) {
	my $line = $lines[$i];
	$line =~ s/Initial release.*/Initial release <$gitrev>/g;
	push(@out, $line);
	$i++;
}

open(FH, ">changelog") || die;
foreach my $line (@out) {
	print FH "$line";
}
close(FH);

print STDOUT "Editing control file\n";
open(FH, "control");
@lines = <FH>;
close(FH);

$i=0;
@out = ();
while($i < scalar(@lines)) {
	my $line = $lines[$i];
	$line =~ s/<insert the upstream URL.*/http:\/\/www.iontorrent.com/;
	$line =~ s/Architecture:.*/Architecture: amd64/;
	$line =~ s/Description:.*/Description: Torrent MAPper/;
	$line =~ s/.*<insert long description.*/ the Torrent MAPping tool/;
	$line =~ s/Section: .*/Section: misc/;
	push(@out, $line);
	$i++;
}

open(FH, ">control") || die;
foreach my $line (@out) {
	print FH "$line";
}
close(FH);

