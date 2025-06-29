#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Thread::Queue;
use File::Path qw(make_path);
use File::Temp qw(tempfile);

# perl src/BAM2exonChain.pl ~/work/ISAtools_comparison/0data/Ref/SIRV_isoforms_multi-fasta_200709a.fasta ./PB.sirv.aligned.sorted.bam test.flnc.exonChain exonChain.count 2 0.99 0.97
# 获取命令行参数
my ($ref_gnome, $bam_file, $flnc_exonChain, $exonChain_count, $num_threads, $min_aln_coverage, $min_aln_identity) = @ARGV;
$num_threads //= 5;  # 默认5线程

# bam过滤
my $temp_file = "$bam_file.filtered.sam";
system("samtools view -@ $num_threads -h -q 20 -F 256 -F 2048 -o $temp_file $bam_file") == 0
    or die "Failed to filter BAM file: $!\n";

## 处理BAM
# 线程队列
my $queue = Thread::Queue->new();

# 读取临时BAM 文件并放入队列
open my $TEMP_BAM, "<", $temp_file or die "Failed to open temp file: $!\n";
while (<$TEMP_BAM>) {
    chomp;
    next if /^@/;  # 跳过头部行
    $queue->enqueue($_);  # 将数据加入队列
}
close $TEMP_BAM;

# 终止信号
for my $i (1 .. $num_threads) {
    $queue->enqueue("THREAD_STOP");
}

# 线程子例程
sub process_bam {
    my %id = ();  # 线程内 read ID 计数

    # 创建临时文件
    my $temp_file = $flnc_exonChain . ".temp" . threads->tid();
    # print "Thread " . threads->tid() . " using temp file: $temp_file\n";  # 调试信息
    open my $TEMP, ">", $temp_file or die "Failed to open temp file: $!\n";

    while (my $line = $queue->dequeue()) {
        last if $line eq "THREAD_STOP";
        my @a = split(/\t/, $line);
        next if $a[2] eq "*";
        next unless ($a[1] == 0 || $a[1] == 16);

        my $astrand = ($a[1] & 0x10) ? '-' : '+';
        my ($ts, $xs, $erro) = ("", "", 0);

        foreach my $tag (@a) {
            $xs  = $1 if $tag =~ /XS:A:([\+\-])/;
            $ts  = $1 if $tag =~ /ts:A:([\+\-])/;
            $erro = $1 if $tag =~ /^NM:i:(\d+)/;
        }

        $xs = ($ts && !$xs) ? (($a[1] & 0x10) ? ($ts eq '-' ? '+' : '-') : $ts) : $xs;
        my $strand = $xs || $astrand;

        my $cigar_str = $a[5];
        my @b = split(/\D/, $cigar_str);
        $cigar_str =~ s/^(\d+)//;
        my @c = split(/\d+/, $cigar_str);

        my @pos;
        my ($gap, $clip, $cov, $seqlen) = (0, 0, 0, length($a[9]));
        for (my $i = 0; $i <= $#c; $i++) {
            if ($c[$i] =~ /(S|H|I)/) {
                $cov  += $b[$i] if $c[$i] eq "I";
                $clip += $b[$i] if $c[$i] =~ /(S|H)/;
                next;
            } elsif ($c[$i] eq "N") {
                my $end = $a[3] + $gap - 1;
                push @pos, $a[3], $end;
                $a[3] = $end + $b[$i] + 1;
                $gap = 0;
            } else {
                $cov  += $b[$i] if $c[$i] =~ /(D|M)/;
                $gap  += $b[$i];
            }
        }

        my $end = $a[3] + $gap - 1;
        push @pos, $a[3], $end;

        my $identity = $cov ? 1 - ($erro / $cov) : 0;
        my $coverage = $seqlen ? ($seqlen - $clip) / $seqlen : 0;

        if ($coverage >= $min_aln_coverage && $identity >= $min_aln_identity) {
            $id{$a[0]}++;
            my $s1 = $pos[0];
            my $e1 = $pos[-1];
            shift @pos;
            pop @pos;
            my $str = join("-", @pos) || "NA";
            print $TEMP "$a[0].m$id{$a[0]}\t$a[2]\t$strand\t$s1\t$e1\t$str\t$identity\t$coverage\n";
        }
    }

    close $TEMP;
    return $temp_file;  # 返回临时文件名
}

# 创建线程池
my @threads = map { threads->create(\&process_bam) } (1 .. $num_threads);

# 收集临时文件名
my @temp_files;
for my $thread (@threads) {
    push @temp_files, $thread->join();
}

# 将所有临时文件合并到最终输出文件
my $temp_files_list = join(" ", @temp_files);
system("cat $temp_files_list > $flnc_exonChain") == 0
    or die "Failed to merge temp files: $!\n";

# 删除临时文件
for my $temp_file (@temp_files) {
    unlink $temp_file or warn "Failed to delete temp file '$temp_file': $!\n";
}


### 统计exonChain.count、提取junction碱基对 不计算平均identity
# 唯一exonChain及其数量
my %ec;
open my $IN2, "<", $flnc_exonChain or die "Failed to open input file: $!\n";
while (<$IN2>) {
    chomp;
    my @fields = split(/\t/, $_);
    my $key = "$fields[1]\t$fields[2]\t$fields[5]";
    $ec{$key}++;
}
close $IN2;


# 读取基因组序列
my %fa;
if ($ref_gnome =~ /\.gz$/) {
    open my $IN, "gzip -dc $ref_gnome |" or die "Failed to open genome file: $!\n";
    my $current_seq = "";
    while (<$IN>) {
        chomp;
        if (/^>(\S+)/) {
            $current_seq = $1;
        } else {
            $fa{$current_seq} .= $_;
        }
    }
    close $IN;
} else {
    open my $IN, "<", $ref_gnome or die "Failed to open genome file: $!\n";
    my $current_seq = "";
    while (<$IN>) {
        chomp;
        if (/^>(\S+)/) {
            $current_seq = $1;
        } else {
            $fa{$current_seq} .= $_;
        }
    }
    close $IN;
}


# 写入 `exonChain_count`
open my $OUT_COUNT, ">", $exonChain_count or die "Failed to open output file: $!\n";
for my $k (keys %ec) {
    my @a = split(/\t/, $k);
    my @b = split(/-/, $a[2]);
    my $m = 0;
    my $ds = "";

    # 跳过值为 "NA" 的键
    next if $a[2] eq "NA";

    # 检查 $fa{$a[0]} 是否存在
    unless (exists $fa{$a[0]} && defined $fa{$a[0]}) {
        warn "Warning: Missing sequence for key '$a[0]' in %fa. Skipping.\n";
        next;
    }

    if ($a[1] eq "+") {
        foreach my $k1 (@b) {
            $m++;
            # 检查 $k1 和 $k1 - 3 是否在有效范围内
            my $len = length($fa{$a[0]});
            if ($k1 < 0 || $k1 >= $len || ($k1 - 3) < 0 || ($k1 - 3) >= $len) {
                warn "Warning: Invalid position '$k1' for key '$a[0]'. Skipping.\n";
                next;
            }
            my $str = ($m % 2 != 0) ? substr($fa{$a[0]}, $k1, 2) : substr($fa{$a[0]}, $k1 - 3, 2);
            $ds .= ($m % 2 != 0) ? "$str-" : "$str,";
        }
    } elsif ($a[1] eq "-") {
        foreach my $k1 (reverse @b) {
            $m++;
            # 检查 $k1 和 $k1 - 3 是否在有效范围内
            my $len = length($fa{$a[0]});
            if ($k1 < 0 || $k1 >= $len || ($k1 - 3) < 0 || ($k1 - 3) >= $len) {
                warn "Warning: Invalid position '$k1' for key '$a[0]'. Skipping.\n";
                next;
            }
            my $str = ($m % 2 != 0) ? substr($fa{$a[0]}, $k1 - 3, 2) : substr($fa{$a[0]}, $k1, 2);
            $str = reverse $str;
            $str =~ tr/ATCGatcg/TAGCtagc/;
            $ds .= ($m % 2 != 0) ? "$str-" : "$str,";
        }
    }

    $ds =~ s/,$//;  # 去掉最后的逗号
    print $OUT_COUNT "$ec{$k}\t$a[0]\t$a[1]\t$a[2]\t$ds\n";
}
close $OUT_COUNT;
