#!/usr/bin/perl

# Adapted from lteng's CalculateProbeQuality_v7.pl
## Generate Probe Quality report using *genome.vcf files, For each probe, calculate
# PercentBaseCoverage>100X: Percentage of bases with coverage >100
# PercentBaseCoverage>250X: Percentage of bases with coverage >250
# MeanBadBases: Average coverage of bases with coverage <100
# Average: Per-base coverage of the probe 
# CV: CV of per-base coverage 
#
# Update aggregate report with BadProbeCount for each sample 
# input: vcf folder, probe information, aggregate report

# First sort the input SampleManifest (index,sample,manifest,probelist,truthpath) file, 
# for each manifest generate probe quality report with corresponding samples and record bad probe count  
#
#
# Previous Analysis: -
# Next Analysis: -
# qli2@illumina.com

# probefile need heading line!
## override the core glob forcing case insensitivity
use File::Glob qw(:globally :nocase);
use File::Basename;
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
## given an array of numbers return the std
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}


$coveragecutoff = 100;
$coveragecutoff2 = 250;
$basepercentcutoff = 0.8;
$Mark = sprintf("%dX",$coveragecutoff);
$Mark2 = sprintf("%dX",$coveragecutoff2);

######
die "Usage $0 alignmentfolder SampleManifest.txt (index,Sample,manifest,probelist,truth) \n"  unless @ARGV == 2;
($dir, $manifest) = @ARGV;
if((substr $dir, -1 ) ne '/'){
   $dir = $dir .'/';
}

%SampleBadProbeCount=();
%Manifests=();
%Samples=();
$fileout = $manifest.".sorted";
`sort -k3,3 -k1,1n -t',' $manifest > $fileout `;

$Total=0;
open(F,"$fileout") or die "cant open $fileout \n";
while(<F>){
   chomp;
   ($index,$sample,$manifest,$probelist) = (split /,/)[0,1,2,3];
   $Total++;
   #print "Read: $manifest\n";
   $Samples{$index}=$sample;
   if(exists $Manifests{$manifest}){
      $ss = $Manifests{$manifest} . ",$sample";
      $Manifests{$manifest} = $ss;	  
   }
   else{  # read a new manifest, generate a subfolder for this manifest
      $Manifests{$manifest} = $probelist.",$sample";	  
	  if(!-e $probelist){
         print "Probe information file $probelist not exist. Probe qualities not calculated.\n";
	     exit;
      }
   }
}
close(F);
# `rm $fileout`;
## Generate probe quality report for each manifest and sample sets 
%BadProbeCount=();
%Bases=();
foreach $manifest (keys %Manifests){
    print "$manifest, $Manifests{$manifest} \n" ;
	$list= $Manifests{$manifest};
	@samples = (split /,/,$list);
	$probelist=$samples[0];
	$numsample = scalar @samples - 1;
	for ($i=1; $i< @samples ; $i++){
        $sample = $samples[$i]; 
        #$sample =~ s/\r|\n//g;
        #$sample =~ s/[^a-zA-Z0-9\-,]+/\-/g;
        #$sample =~ s/[\-]+/-/g;				  		
        $fileg= $dir.$sample.".stitched.genome.vcf";
		# print "Looking for $fileg\n"
		$filep= $dir.$sample.".stitched.phased.genome.vcf";
		if(-e $fileg){
            $filevcf = $fileg;
			print "Calculating probe quality $sample\n";
        }
        elsif(-e $filep){
            $filevcf = $filep;
			print "Calculating probe quality $sample\n";
        }
		else{
            print "vcf file not found in $dir. Probe quality not calculated.\n";
	        exit;
        }
        
        open(F,"$filevcf") or die "can't open $filevcf \n";
        while(<F>){
            chomp;
	        ($chr, $pos, $alternate, $format, $info) = (split/\t/)[0,1,4,7,9];
	        @items = (split /;/,$format);
	        foreach $item (@items){
	            if($item =~/^DP=/){
		            $DP = substr $item,3;
		        }
	        }					    
	        $frequency = (split /:/,$info)[3];		
	        $Bases{$sample}{$chr}{$pos} = $DP;	      
        }
        close(F);
		print "$sample \n";
        $BadProbeCount{$sample}=0;
	}  

	######### Process the probes ####
	$name = basename($manifest);
	$pp = (split /\_/, $name)[0];	  
    $outfile= $dir."ProbeQualityReport_$pp.csv";	
    open(OUT,">$outfile") or die "cant open $outfile\n";
    open(F, "$probelist") or die "cant open $probefile\n";
    $index = 0;
    while(<F>){
        chomp;
	    $_ =~ s/\r|\n//g;	
	    if(/^#/){
	        @items = (split /[,\t]/,$_);
	        $ss = join(',', @items);
	        print OUT $ss;
	        #print OUT $_;
	        for ($i=1; $i<@samples; $i++){
                $sample = $samples[$i];	
	            print OUT ",#BaseinManifest($sample),PerbaseMean($sample),CV($sample),%>$Mark($sample),%>$Mark2($sample),BadProbeFlag(20%<100X)";
	        }
	        print OUT ",TotalBadProbePercent\n";
	        next;
		}
	    $index++;		
        ($chr,$start,$end)= (split /[,\t]/)[1,2,3];
	    if($chr !~ /chr/){
	        $chr = 'chr'.$chr;
	    }
        $badflagcount=0;      
	    print OUT "$_";
        # Each time report performance of one probe 
	    # loop through the samples 
	    # print "$chr, $start, $end\n";
	    for ($i=1; $i<@samples; $i++){
            $sample = $samples[$i];			
	        $num100X = 0;
	        $num250X = 0;
            $numbadbase = 0;
            $meanbadbase = 0;
	        $coveredbase = 0;
	        @PerbaseMean=();        
	        if($dedupmean == 0 ){
	            $dedupCV = 'NA';
            }
            else{
	            $dedupCV = $dedupstd/$dedupmean;
            }

            for($k=$start+1; $k<=$end; $k++){
                if(exists $Bases{$sample}{$chr}{$k}){
                    $DP = $Bases{$sample}{$chr}{$k};
			        $coveredbase++;
			        push @PerbaseMean, $DP;
			        if( $DP > $coveragecutoff){
	                    $num100X++
	                }
  	                if ( $DP > $coveragecutoff2){
	                    $num250X++
	                }				
		            else{
			            $numbadbase++;
			            $meanbadbase = $meanbadbase + $DP;
	                }
                }
                else{
		            #$numbadbase++;
	            }
	        }
	        if($numbadbase == 0){
	            $MeanBadBaseCoverage = $coveragecutoff;		  
	        }
	        else{
	            $MeanBadBaseCoverage = $meanbadbase/$numbadbase;		  
	        }
	        #$BaseCoverage100X=$num100X / ($end-$start+1);
	
	        if ($coveredbase == 0 ){
	            $BaseCoverage100X = -1;
		        $BaseCoverage250X = -1;
		        $perbasemean = -1;
	            $perbasestd = -1;
	            $perbasecv = -1;		
	        }
	        else{
	            $BaseCoverage100X=$num100X / $coveredbase;
		        $BaseCoverage250X=$num250X / $coveredbase;
		        $perbasemean = &average(\@PerbaseMean);
	            $perbasestd = &stdev(\@PerbaseMean);
		        if($perbasemean > 0 ){
		            $perbasecv = $perbasestd/$perbasemean;		
		        }
		        else{
		            $perbasecv = -1;		
		        }		    
	        }		
		
	        $badflag = 0;
	        if($BaseCoverage100X<$basepercentcutoff){
	            $badflag = 1;
	            $badflagcount++;
	            $BadProbeCount{$sample}++;
	        }
	        print OUT sprintf(",%d,%d,%.2f,%.2f%%,%.2f%%,%d",$coveredbase,$perbasemean, $perbasecv, $BaseCoverage100X*100, $BaseCoverage250X*100,$badflag);	
        }	
        print OUT sprintf(",%.0f%%\n",$badflagcount/$numsample*100);	
    }
    close(F);
    print "Generate quality information for $index probes in $outfile.\n";
    close(OUT);
} 