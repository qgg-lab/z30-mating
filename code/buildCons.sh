# Build consensus sequence to be mapped

# Indels cannot be easily handeled so they are ignored
# SNPs/MNPs are handeled using two differnet approaches and then compare
# 1. consensus sequence, where all lines differ from the reference
# 2. in addition to consensus, put a third base on the position where the lines are different

# extract variant information from the VCF file
# output startLines.variant.bed
# in the format of 
# chr start end ref alt alt_count total_called_number_lines

perl -we 'use List::Util qw(max sum); $header = <>; chomp $header; @line = split /\t/, $header;
  open ID, "<startLines.txt";
  %id = ();
  while (<ID>) {
    chomp $_;
    $id{"DGRP-".$_} = 1;
  }
  close ID;
  @idIndex = ();
  for ($i = 9; $i <= $#line; $i++) {
    if (defined($id{$line[$i]})) {
      push(@idIndex, $i);
    }
  }
  while (<>) { chomp $_; @line = split /\t/, $_;
    @qualcov = split /:/, $line[5];
    $qual = $qualcov[0];
    if ($qual <= 500) { next; }
    if (length($line[3]) != length($line[4])) { 
      next;
    }
    @geno = ();
    for ($i = 0; $i <= $#idIndex; $i++) {
      ($geno, $sr, $or, $gq) = split /:/, $line[$idIndex[$i]];
      $dp = $sr + $or;
      if ($dp > 0 && $gq > 20) {
        if ($geno eq "0/0") {
          push(@geno, 0);
        } elsif ($geno eq "1/1") {
          push(@geno, 2);
        } elsif ($geno eq "1/0" || $geno eq "0/1") {
          push(@geno, "1");
        }
      }
    }
    $total = $#geno + 1;
    $alt = sum(@geno);
    if ($total == 0 || $alt == 0) { next; }
    print $line[0], "\t", $line[1] - 1, "\t", $line[1] - 1 + length($line[3]), "\t", $line[3], "\t", $line[4], "\t", $alt, "\t", $total, "\n";
  }' ~/storage/dgrpFreeze2/andreasGeno2/freeze2.single.jgil.hpp.biallelic.vcf > startLines.variant.bed 2> startLines.variant.log &

# make mended genome sequence
perl -we 'open REF, "</storage/whuang9/flybase/fb-r5.49/dmel.rmsk.fa"; %refseq = ();
  while (<REF>) {
    chomp $_;
    if ($_ =~ m/^>(.*)/) {
      $chr = $1;
      $refseq{$chr} = "";
      next;
    }
    $refseq{$chr} .= uc($_);
  }
  close REF;
  # new reference;
  %newref = %refseq;
  %ntcode = (
             AC => "T",
	     CA => "G",
	     AG => "C",
	     GA => "C",
	     AT => "C",
	     TA => "G",
	     CG => "T",
	     GC => "T",
	     CT => "A",
	     TC => "G",
	     GT => "A",
	     TG => "A",
             AA => "A",
             CC => "C",
             GG => "G",
             TT => "T"
            );
  open CHANGE, "<startLines.variant.bed";
  open TRACK, ">ref/startLines.variant.mend.track";
  while (<CHANGE>) {
    chomp $_;
    @line = split /\t/, $_;
    if ($line[5] < 4 || $line[6] < 20) { next; } # require that at least 20 lines typed and 4 alleles observed (could be two homozyous lines)
    @refallele = split //, $line[3];
    @altallele = split //, $line[4];
    $change = "";
    for (my $i = 0; $i <= $#refallele; $i++) {
      $code = $refallele[$i].$altallele[$i];
      $change .= $ntcode{$code};
    }
    print TRACK join("\t", @line[0..4]), "\t", $change, "\n";
    substr($newref{$line[0]}, $line[1], $line[2] - $line[1], $change);
  }
  foreach my $chr (sort(keys(%newref))) {
    print ">$chr\n";
    for ($i = 0; $i <= int((length($refseq{$chr}) - 1)/50); $i++) {
      print substr($newref{$chr}, $i*50, 50), "\n";
    }
  }' > ref/startLines.mended.fa 2> ref/startLines.mended.log &

# make consensus sequence
perl -we 'open REF, "</storage/whuang9/flybase/fb-r5.49/dmel.rmsk.fa"; %refseq = ();
  while (<REF>) {
    chomp $_;
    if ($_ =~ m/^>(.*)/) {
      $chr = $1;
      $refseq{$chr} = "";
      next;
    }
    $refseq{$chr} .= uc($_);
  }
  close REF;
  # new reference;
  %newref = %refseq;
  open CHANGE, "<startLines.variant.bed";
  open TRACK, ">ref/startLines.variant.cons.track";
  while (<CHANGE>) {
    chomp $_;
    @line = split /\t/, $_;
    if ($line[5] < 4 || $line[6] < 20 || $line[5] != $line[6]*2) { next; } # require that at least 20 lines typed and 4 alleles observed (could be two homozyous lines)
    print TRACK join("\t", @line[0..4]), "\t", $line[4], "\n";
    substr($newref{$line[0]}, $line[1], $line[2] - $line[1], $line[4]);
  }
  foreach my $chr (sort(keys(%newref))) {
    print ">$chr\n";
    for ($i = 0; $i <= int((length($refseq{$chr}) - 1)/50); $i++) {
      print substr($newref{$chr}, $i*50, 50), "\n";
    }
  }' > ref/startLines.cons.fa 2> ref/startLines.cons.log &

# make bwa index
~/software/bwa-0.6.2/bwa index -p ref/startLines.cons ref/startLines.cons.fa > ref/startLines.cons.bwabuild.log 2>&1 &
~/software/bwa-0.6.2/bwa index -p ref/startLines.mended ref/startLines.mended.fa > ref/startLines.mended.bwabuild.log 2>&1 &

~/software/samtools-0.1.18/samtools faidx ref/startLines.cons.fa
~/software/samtools-0.1.18/samtools faidx ref/startLines.mended.fa
