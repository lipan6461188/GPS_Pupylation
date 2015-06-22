#!/opt/local/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;

#TP : 311
#TN : 3859

my $ubound; #阈值上界
my $dbound; #阈值下界
my $gap;
my $blosum;
my $sample;
my $output; #输出文件
my $silent; #氨基模式

Getopt::Long::GetOptions('ubound=f' => \$ubound, 'dbound=f' => \$dbound, 'gap=f' => \$gap,
                        'blosum=s' => \$blosum, 'sample=s' => \$sample, 'output=s' => \$output, 'silent' => \$silent)
or die("Error in command line arguments\n");

#读入文件
open(HANDLE, "<$sample") || die "打开文件失败\n";
chomp(my @all_seq = <HANDLE>);
close HANDLE;

#正负样本分类
our @positive = ();
our @negative = ();
our @positive_score = ();
our @negative_score = ();
foreach my $eachSeq(@all_seq)
{
    my @seq_prop = split(/\t/,$eachSeq);
    if($seq_prop[1] eq "+")
    { push(@positive, $seq_prop[0]); }
    else
    { push(@negative, $seq_prop[0]); }
}

say '正样本个数：',$#positive+1;
say '负样本个数: ',$#negative+1;
our $pos_num = $#positive+1;
our $neg_num = $#negative+1;

#获取BLOSUM62矩阵
our %BLUSOM = constructBLOSUM62();

#权重向量
#our @weights = (1, 0, 0, 3, 2, 2, 3, 1, 1, 1, 1, 0, 1, 1, 2, 3, 1, 0, 1, -1, 0, 1, 0, 2, 1, 1, 3);

our @weights = (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);

my $pos_s = 0;
#计算正样本
foreach my $one (@positive)
{
    #  say $temp;
    my $sum = 0;#定义一个和，最后要求平均。
    foreach my $each (@positive)
    {
        #不要求自己与自己的相似度，因为这个结果会很高
        if( $one eq  $each) { next; }
        #求两条序列之间的距离
        my $s_grade = similarity($one, $each);
        #如果序列的距离小于0，就设置为0，只把分数大于0的序列分数相加
        if($s_grade > 0) { $sum += $s_grade;}
        #say "正样本$one与$$each比较的分数".$s_grade;
    }
    #say "Ave: ".$average_grade;
    #求一条序列与其他序列距离的平均值，把这个平均值放置到数组中，今后用于统计
    push(@positive_score, $sum/($#positive+1));
    
    $pos_s += $sum;
    # say "正样本$one的得分是: ".$sum;
}

my $neg_s = 0;
#计算负样本
foreach my $one (@negative)
{
    my $sum = 0;#定义一个和，最后要求平均。
    foreach my $each (@positive)
    {
        #求两条序列之间的距离
        my $s_grade = similarity($one, $each);
        #如果序列的距离小于0，就设置为0
        if($s_grade > 0) { $sum += $s_grade;}
        
        #   say "正样本$one与负样本$each比较的分数".$s_grade;
    }
    #say "Ave: ".$average_grade;
    #求一条序列与其他序列距离的平均值，把这个平均值放置到数组中，今后用于统计
    push(@negative_score, $sum/($#positive+1));
    
    $neg_s += $sum;
    # say "负样本$one的得分是: ".$sum;
}

say '$pos_s: ',$pos_s/$#positive;
say '$neg_s: ',$neg_s/$#negative;
say "差值为： ",$pos_s/$#positive - $neg_s/$#negative;


#统计预测情况
open OUTPUT,">$output" || die "Open FIle Failed\n";
for(my $threshold=$dbound; $threshold <= $ubound; $threshold += $gap)
{
    
    #初始化这些变量
    my $TP=0;
    my $FN=0;
    my $TN=0;
    my $FP=0;
    my $Sn=0;
    my $Sp=0;
    my $Pr=0;
    
    foreach(@positive_score)
    {
        if($_ >= $threshold)
        { $TP++; }
        else
        { $FN++; }
    }
    
    foreach(@negative_score)
    {
        if($_ < $threshold)
        { $TN++; }
        else
        { $FP++; }
    }

    $Sp = $TP/($TP+$FN);
    $Sn = $TN/($TN+$FP);
    if($TP+$FP eq 0)
    {
        $Pr = 0;
    }else
    {
        $Pr = $TP/($TP+$FP);
    }
    
    #是否向屏幕输出
    if(!$silent)
    {
        say "\n==>阈值是<==：",$threshold;
        say '$TP: '.$TP;
        say '$FN: '.$FN;
        say '$TN: '.$TN;
        say '$FP: '.$FP;
        
        print "Sn Sp Pr: ";
        print $Sp."\t";
        print $Sn."\t";
        print $Pr."\n";

    }
    
    print OUTPUT 1-$Sp."\t";
    print OUTPUT $Sn."\t";
    print OUTPUT $Pr."\n";

}

close OUTPUT;


#子函数用于构造哈希
sub constructBLOSUM62
{
    #读入哈希
    open HASH,"<$blosum";
    chomp(my @file = <HASH>);
    close HASH;
    
    my %BLOSUM62;
    my @index;
    
    foreach my $eachline(@file)
    {
        #截取信息
        my @info = split(/ +/,$eachline);
        #第一行的信息
        if($info[0] eq "")
        {
            shift @info;
            @index = @info;
            next;
        }
        #其余行的信息
        my $tablet = shift @info;
        for(my $i=0; $i < 24; $i++)
        {
            $BLOSUM62{ $tablet.$index[$i] } = $info[$i];
        }
    }
    return %BLOSUM62;
}

#用于确定两条序列的相似度
sub similarity
{
    my $seq1 = pop @_;
    my $seq2 = pop @_;
    if(length($seq1) ne length($seq2))
    {
        say '$seq1: ',$seq1;
        say '$seq2: ',$seq2;
        warn("==>要比较的序列长度不等<==\n");
        exit(0);
    }
    #  say "==============";
    my $score = 0;
    for (my $var = 0; $var < length($seq1); $var++)
    {
        my $bp1 = substr ($seq1,$var,1);
        my $bp2 = substr ($seq2,$var,1);
        $score += $BLUSOM{ $bp1.$bp2 } * $weights[ $var ];
        #    say $weights[ $var ];
    }
    
    return $score;
}






