print "Enter the input file name (without .txt):\n";
$in=<STDIN>;
chomp $in;
print "\nloading data...\n";
open(DAT, "${in}.txt") || die "no suchafile";
open(GATC, "All_GATC_list.txt") || die "no suchafile";
open(OUT, ">${in}_reads_per_GATC_filtered.txt");

@gatc=<GATC>;
@data=<DAT>;
close (DAT);
close (GATC);

print "\nsorting data...\n";

@data = sort					
{                                               

(split(" "||"\t",$a))[0] cmp (split(" "||"\t",$b))[0]     
||                                                 
(split(" "||"\t",$a))[1] <=> (split(" "||"\t",$b))[1]     
  
}@data;                                         

@gatc = sort					
{                                               

(split("\t",$a))[0] cmp (split("\t",$b))[0]     
||                                                 
(split("\t",$a))[1] <=> (split("\t",$b))[1]     
  
}@gatc;    




print "\nprocessing data...\n";
$lastdot=0;
for($r=0;$r<=$#gatc;$r++)
{
        my ($g_chr,$g_st,$g_end)=split ("\t",$gatc[$r]);
        my @idhash=();
        my $idc=0;
        my $counter=100*sprintf("%.5f", ($r/$#gatc));
        #print "\n$counter% done\t$found GATC reads\t$g_chr";
        my $count=0;
        chomp $g_end;
        my $chrfound=0;
        for($i=$lastdot; $i<=$#data; $i++)
        {
                my ($chr1,$st1,$end1,$id1,$col, $str, @trash1)=split(" "||"\t",$data[$i]);
                chomp $str;
                $id1=~/(.*)(\/)(.)/;
                if($3)
                {
                        $id1=$1;
                }
                if (($g_chr eq $chr1) and ($g_end<$st1))
                {
                        
                        last
                }
                elsif($g_chr lt $chr1)
                {
                        last
                }
                if($g_chr eq $chr1)
                {
                        $chrfound=1;
                        
                        
                        
                        if($str eq '+')
                        {
                                if(abs($g_st-$st1-1)<=3)
                                {
                                        $lastdot=$i+1;
                                        $found++;
                                        $idhash[$idc]=$id1;
                                        $idc++;
                                        $count++;
                                        next;
                                        
                                        
                                }
                        }
                        elsif($str eq '-')
                        {
                                if(abs($g_end-$end1+2)<=3)
                                {
                                        $lastdot=$i+1;
                                        $found++;
                                        $idhash[$idc]=$id1;
                                        $idc++;
                                        $count++;
                                        next;
                                        
                                        
                                }
                        }    
                        elsif($st1>$g_st)
                        {
                                
                                last
                        }
                }elsif($chrfound==1)
                {last}
                
        }
        #убираем повторенные id
        my @s_idhash=sort @idhash;

        for($k=1;$k<=$#s_idhash;$k++)
        {
                if($s_idhash[$k] eq $s_idhash[$k-1])
                {
                        $count--;
                        $found--;
                }
        }
        ###
        print OUT "$g_chr\t$g_st\t$g_end\t$count\n";
      
}

print OUT "\n\n\t$found reads mapped to GATC border";
print "\nDone\n\n$found reads mapped to GATC border\n";
sleep 5; 
