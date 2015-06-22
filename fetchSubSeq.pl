#!/opt/local/bin/perl -w
use 5.010;
use strict;
use Bio::DB::Fasta;
use Getopt::Long;       #长参数

#定义变量接受参数
our $upper = 15;
our $down = 15;
our $patter = "K"; #匹配的模式,默认是K
our $display_style = "tab_delimit";
our $delimit = 1; #碱基之间的分割样式，1表示不分割，0表示tab
my  $help = 0;
our $output;  #序列文件
our $positive;#正样本文件
our $input;   #输入文件

Getopt::Long::GetOptions('upper=i' => \$upper, 'down=i' => \$down,
                         'patter=s' => \$patter, 'delimit=i' => \$delimit,
                         'help' => \$help, 'output=s' => \$output,
                        'input=s' => \$input, 'positive=s' => \$positive)
                        or die("Error in command line arguments\n");

if(!defined($input) || !defined($output)) { warn("Please Specify input and output file\n"); exit(-1); }

if($help)
{
    say "--upper    上游氨基酸个数";
    say "--down    下游氨基酸个数";
    say "--patter   匹配的氨基酸模式";
    say "--delimit  结果分割样式";
    say "--input    输入序列文件";
    say "--output   输出样本集";
    say "--positive 正样本文件";
    say "--help     帮助";
    exit(0);
}

if ($delimit eq "1")
{
    $display_style = "no_delimit";
}

#取较大的数来填充星号
my $starNum = ($upper > $down ? $upper : $down) + length($patter); #星号个数

say "填充星号数: ".$starNum;
say "分割样式: ".$display_style;
say "上游: ".$upper;
say "下游: ".$down;

#读入所有的蛋白质序列
my $proteins = Bio::DB::Fasta->new($input);
my @ids = $proteins->ids;

open(POSITIVE,$positive) || die "无法正样本文件\n";
my @positive_samples = <POSITIVE>;
close POSITIVE;

#定义哈希来存放所有的数据
my %allSeq;
foreach my $eachSeq(@positive_samples)
{
    my @id_pos = split(/\t/,$eachSeq);
    
    #由id的名字得到序列
    my $eachSeq = $proteins->seq($id_pos[0]);
    
    #前面和后面都加上星号
    $eachSeq = ('*'x$starNum).$eachSeq.('*'x$starNum);
    #say $eachSeq;

    while($eachSeq =~ /$patter/g)
    {
        #使用pos函数把所有的情况都取出来
        my $pos = pos($eachSeq); #K的位置在 pos-1
        #say substr($eachSeq,$pos-16,31);
        my $mySeq = substr($eachSeq,$pos-$upper-1,$upper+$down+length($patter));
        if( $id_pos[1] == $pos - $starNum ){
            #上游与下游各取多少位
            $allSeq{$mySeq} = "+";
        }else{
            if(!defined($allSeq{$mySeq})){
                $allSeq{$mySeq} = "-";
            }
        }
    }

}

open(SEQ,">$output") || die "打开输出样本失败\n";
while( (my $seq, my $type) = each %allSeq )
{
    if( $display_style eq "tab_delimit" )
    {
        my @AAs = split(//,$seq);
        my $AAs = "";
        foreach(@AAs){ $AAs = $AAs.$_."\t"; }
        print SEQ $AAs.$type."\n";
    }
    if( $display_style eq "no_delimit" )
    {
        print SEQ $seq."\t".$type."\n";
    }
}
close SEQ;

#删除文件
unlink("$input.index");

exit(0);