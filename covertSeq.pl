#!/opt/local/bin/perl -w
use 5.010;
use strict;
use Bio::DB::Fasta;
use Getopt::Long;        #接受参数的模块

#   这个文件主要是把从Uniprot下载的
#   序列转成更合适的Fasta格式

#定义变量接受参数
my $input;
my $output;
my $help;
Getopt::Long::GetOptions('input=s' => \$input, 'output=s' => \$output, 'help' => \$help)
or die("Error in command line arguments\n");

if($help)
{
    say "--input  输入序列文件";
    say "--output 输出序列文件";
    exit(0);
}

if(!defined($input) || !defined($output)){ warn("Please Specify input and output\n"); exit(-1); }

#读入所有的蛋白质序列
my $old_proteins = Bio::DB::Fasta->new($input);

#say "蛋白质的序列数为：".$old_proteins->seq();

open(CONVERT,">$output") || die "打开文件失败\n";
#首先转一下序列，把序列的>后面转成ID
my $stream = $old_proteins->get_PrimarySeq_stream;
while(my $seq = $stream->next_seq)
{
    #sp|A0QV10|Y2408_MYCS2 截取每个ids的第二位作为头
    my @names = split(/\|/,$seq);
    print CONVERT ">".$names[1]."\n";
    print CONVERT $seq->seq."\n";
}

close CONVERT;

unlink("$input.index") or die "删除文件失败\n";
exit(0);
