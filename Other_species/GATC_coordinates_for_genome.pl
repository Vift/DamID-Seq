print "Enter the whole genome fasta file name (without .fa or fa.gz):\n";
$in=<STDIN>;
chomp $in;
print "uncompressing data\n";
use IO::Uncompress::Gunzip qw($GunzipError);
$input="${in}.fa.gz";
my $gzip=1;
my @genome=();
my $z = IO::Uncompress::Gunzip->new( $input )  or $gzip=0;
my @gatc_list;
if($gzip==1)
{
    @genome=<$z>;
}
else
{
    open (GEN, "${in}.fa") || open (GEN, "${in}.fasta") || die "can't find file ${in}.fa.gz or ${in}.fa";
    @genome=<GEN>;
}
open (OUT, ">All_GATC_list.txt");
my $chr='none';
foreach(@genome)
{
    if ($_=~/>([A-Za-z0-9_]+)/)
    {
        unless($chr eq 'none')
        {
            print "\nprocessing data for ${chr}\n";
        }

        search ($chr,$seq);
        $seq="";
        $chr=$1;
        unless($chr=~/>chr([A-Za-z0-9_]+)/)
        {
            $chr=join("",("chr",$chr))
        }        
        
        
        next
    }
    chomp;
	s/\s//g;
	$seq.= $_;
}
search ($chr,$seq);

sub search
{
	my ($chr, $seq) = @_;
    my $st=0;
    for(my $b=0;$b<=length($seq)-4;$b++)
    {
        my $quard = substr($seq,$b,4);
		if ($quard =~ m/GATC/gi)
        {
            if(($b+2-$st)<10000)
            {
                push @gatc_list, join("\t", ($chr,$st,($b+2)));
            }
            $st=($b+2);
		} 
    }
}
print "\nwriting output file...\n";
@gatc_list = sort	
   {

   (split("\t",$a))[0] cmp (split("\t",$b))[0]
||
   (split("\t",$a))[1] <=> (split("\t",$b))[1]
   
   }@gatc_list;

foreach(@gatc_list)
{
    print OUT $_."\n";
}