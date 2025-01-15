open IN,"$ARGV[0]";
%tr=();
while(<IN>){
	chomp;
	@a=();@a=split(/\t/,$_);
	@b=();@b=split(/\;/,$a[8]);
	$gi="NA";$gn="NA";$ti="NA";
	if($a[2] eq "exon"){
		foreach $k(@b){
			$gi=$k if($k=~/gene_id/);
			$gn=$k if($k=~/gene_name/);
			$ti=$k if($k=~/transcript_id/);
		}
		$gi=~s/^ //;
		$gi=~s/gene_id //;
		$gi=~s/"//g;
		$gn=~s/^ //;
		$gn=~s/gene_name //;
		$gn=~s/"//g;
		$ti=~s/^ //;
		$ti=~s/transcript_id //;
		$ti=~s/"//g;
		$tr{$ti}[0]="$gi\t$gn";
		$tr{$ti}[1]="$a[0]\t$a[6]";
		push @{$tr{$ti}[2]},$a[3];
		push @{$tr{$ti}[2]},$a[4];

		# print "$gi\t$gn\t$ti\t$a[0]\t$a[6]\n";
	}
}
close IN;

print "TrID\tGeneID\tGeneName\tChr\tStrand\tTrStart\tTrEnd\texonChain\n";
foreach $k1(keys %tr){
	@a=();@a=split(/\t/,$tr{$k1}[1]);
	
	@tmp=();@tmp=sort{$a<=>$b}@{$tr{$k1}[2]};

	$start=$tmp[0];$end=$tmp[-1];
	$tmp[0]="";$tmp[-1]="";
	$str=join("-",@tmp);
	$str=~s/^-//;
	$str=~s/-$//;
	$str="NA" if($str eq "");
		
	print "$k1\t$tr{$k1}[0]\t$tr{$k1}[1]\t$start\t$end\t$str\n";
}
