#! /usr/bin/perl
# perl BAM2exonChain.pl ref_gnome.fa test.bam test.flnc.exonChain test.exonChain.count

open IN,"$ARGV[0]" if($ARGV[0]!~/gz$/);
open IN,"gzip -dc $ARGV[0] | " if($ARGV[0]=~/gz$/);
%fa=();
while(<IN>){
	chomp;
	if($_=~/^>/){
		$_=~s/^>//;
		@a=();@a=split(/ /,$_);
	}else{
		$fa{$a[0]}.=$_;
	}
}
close IN;

open IN,"samtools view $ARGV[1] | ";
open OUT,">$ARGV[2]";
%id=();%ec=();%eci=();
while(<IN>){
	chomp;
	next if($_ =~ /^@/);
	@a=();@a = split(/\t/,$_);
	next if($a[2] eq "*");
	
	# 只保留唯一比对
	next unless ($a[1] == 0 || $a[1] == 16);
	
	$astrand = (($a[1] & 0x10)==0) ? '+' : '-';
	$ts="";$xs="";$tag="";$erro=0;
	foreach $tag (@a) {
		if ($tag=~/XS:A:([\+\-])/) {
			$xs=$1;
		}
		if ($tag=~/ts:A:([\+\-])/) {
			$ts=$1;
		}
		if ($tag=~/^NM:i:(\S+)/){
			$erro=$1;
		}
	}

	if ($ts && !$xs) {
		$xs=$ts;
		if (($a[1] & 0x10)!=0) {
			$xs = ($ts eq '-') ? '+' : '-';
		}
	}
 	$strand=$xs || $astrand;
	
	@b = split(/\D/, $a[5]);
	$a[5] =~ s/^(\d+)//;
	@c = split(/\d+/,$a[5]);
	$end = 0;$gap = 0;@pos = ();
	$clip=0;$cov=0;
	for($i=0;$i<=$#c;$i++){
		if($c[$i] =~ /(S|H|I)/){
			$cov+=$b[$i] if($c[$i] eq "I");
			$clip+=$b[$i] if($c[$i]=~/(S|H)/);
			next;
		}elsif($c[$i] eq "N"){
			$end = $a[3] + $gap -1;
			push @pos,$a[3];
			push @pos,$end;
			$a[3] = $end + $b[$i] + 1;
			$gap = 0;
		}else{
			$cov+=$b[$i] if($c[$i] =~/(D|M)/);
			$gap += $b[$i];
		}
	}
	$end = $a[3] + $gap -1;
	
	push @pos,$a[3];
	push @pos,$end;
	
	$identity=0;$identity=1-($erro/$cov);
	$seqlen=0;$seqlen=length($a[9]);
	$coverage=0;$coverage=($seqlen-$clip)/$seqlen;
	
	$id{$a[0]}++;
	
	$s1=$pos[0];$e1=$pos[-1];
	$pos[0]="";$pos[-1]="";
	$str=join("-",@pos);
	$str=~s/^-//;
	$str=~s/-$//;
	$str="NA" if($str eq "");
	print OUT "$a[0].m$id{$a[0]}\t$a[2]\t$strand\t$s1\t$e1\t$str\t$identity\t$coverage\n";
	$str2="$a[2]\t$strand\t$str";
	$ec{$str2}++ if($str ne "NA");
	$eci{$str2}+=$identity if($str ne "NA");
}
close IN;
close OUT;

open OUT,">$ARGV[3]";
foreach $k(keys %ec){
	$idmean=0;$idmean=$eci{$k}/$ec{$k};
	
	@a=();@a=split(/\t/,$k);
	@b=();@b=split(/\-/,$a[2]);
	$m=0;$ds="";
	if($a[1] eq "+"){
		foreach $k1(@b){
			$m++;
			$str=substr($fa{$a[0]},$k1,2) if ($m % 2 != 0);
			$str=substr($fa{$a[0]},$k1-3,2) if ($m % 2 == 0);
			$ds.="$str-" if ($m % 2 != 0);
			$ds.="$str," if ($m % 2 == 0);
		}
	}elsif($a[1] eq "-"){
		foreach $k1(reverse @b){
			$m++;
			$tmp=substr($fa{$a[0]},$k1-3,2) if ($m % 2 != 0);
			$tmp=substr($fa{$a[0]},$k1,2) if ($m % 2 == 0);
			$str=reverse $tmp;
			$str=~tr/ATCGatcg/TAGCtagc/;
			$ds.="$str-" if ($m % 2 != 0);
			$ds.="$str," if ($m % 2 == 0);
		}
	}
	$ds=~s/\,$//;
	print OUT "$ec{$k}\t$a[0]\t$a[1]\t$a[2]\t$ds\t$idmean\n";
}
close OUT;
