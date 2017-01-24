use Text::NSP::Measures::2D::Fisher::left;
print "Press Enter to use automatic P-value cutoff level. Alternatively type the desired manual P-value cutoff (from 0 to 1):\n";
my $manpval=<STDIN>;
chomp $manpval;
if ($manpval and ($manpval>1 or $manpval<=0))
{
    print "incorrect P-value. Should be from 0 to 1.\nPress Enter to exit";
    my $exit=<STDIN>;
    die 
}
if (!$manpval)
{
    $target_fdr=0.05;
    print "Press Enter to use standard FDR cutoff level (0.05). Alternatively type the desired manual FDR cutoff (from 0 to 1):\n";
    my $tmpfdr=<STDIN>;
    chomp $tmpfdr;
    if ($tmpfdr and ($tmpfdr>1 or $tmpfdr<=0))
    {
        print "incorrect FDR-value. Should be from 0 to 1.\nPress Enter to exit";
        my $exit=<STDIN>;
        die 
    }
    if ($tmpfdr)
    {
        $target_fdr=$tmpfdr;
    }
}
my $damx_over_dam_lim=1;
print "Press Enter to skip minimal Dam-X over Dam excess level. Alternatively type the minimal Dam-X/Dam level (greater than 1) for GATC fragment to become a peak:\n";
my $tmpex=<STDIN>;
chomp $tmpex;
if ($tmpex and $tmpex<1)
{
    print "incorrect Dam-X over Dam excess level. Should be greater or equal to 1.\nPress Enter to exit";
    my $exit=<STDIN>;
    die
}
if ($tmpex)
{
    $damx_over_dam_lim=$tmpex;
}
my $lim_div=sprintf("%.3f", (log($damx_over_dam_lim))/log(2));
print "Enter DamID experiment name (for ex: Protein1_Ovary):\n";
$in5=<STDIN>;
chomp $in5; 

print "Enter the combined Dam datafile name (without _filtered_combined.txt):\n";

$in1=<STDIN>;
chomp $in1;
open(DAM, "${in1}_filtered_combined.txt")  || die "no suchafile";
open(DAMN, "${in1}_filtered_combined.txt")  || die "no suchafile";
print "Enter the averaged Dam-X datafile name (without _filtered_combined.txt):\n";

$in2=<STDIN>;
chomp $in2;
open(DAMX, "${in2}_filtered_combined.txt")  || die "no suchafile";
open(DAMXN, "${in2}_filtered_combined.txt")  || die "no suchafile";

my $found_prob_level=0;
@damn=<DAMN>;
@damxn=<DAMXN>;
@dam=<DAM>;
@damx=<DAMX>;

############################################################################################################auto-pval
if (!$manpval)
{
print "Enter the first Dam-X repeat datafile name (without _filtered.txt):\n";    
$in3=<STDIN>;
chomp $in3;
open(F, "${in3}_filtered.txt")  || die "no suchafile";

print "Enter the second Dam-X repeat datafile name (without _filtered.txt):\n";

$in4=<STDIN>;
chomp $in4;
open(S, "${in4}_filtered.txt")  || die "no suchafile";

print "Wait...\n";

open(OUT, ">${in5}_ROC.txt");
$previous_fh = select(OUT);	# запись в выходной файл (обозванный тут OUT) сразу, без буферизации, вставл€ть после указани€ выходного файла
$| = 1;
select($previous_fh);





@first=<F>;
@second=<S>;

if ($#dam!=$#damx or $#damx!=$#first or $#first!=$#second)
{
    print "Input files have different length! Script will stop...\n";
    my $stop=<STDIN>;
    exit;
}


##########
my $c_dam=0;
my $c_damx=0;
my $c_first=0;
my $c_second=0;
for(my $c=0;$c<=$#dam;$c++)
{
    my ($chrd,$std,$endd,$vald)=split("\t",$dam[$c]);
    my ($chrx,$stx,$endx,$valx)=split("\t",$damx[$c]);
    my ($chrf,$stf,$endf,$valf)=split("\t",$first[$c]);
    my ($chrs,$sts,$ends,$vals)=split("\t",$second[$c]);
    chomp ($vald,$valx,$valf,$vals);
    if ($vald eq "NA" or $valx eq "NA")
    {
        next
    }
    
    $c_dam+=$vald;
    $c_damx+=$valx;
    $c_first+=$valf;
    $c_second+=$vals;    
}

#########

my @true=();
my @false=();
my @true_dost=();
my @false_dost=();

#####


for(my $i=0;$i<=$#dam;$i++)
{
    my ($chrd,$std,$endd,$vald)=split("\t",$dam[$i]);
    chomp $vald;
    my ($chrx,$stx,$endx,$valx)=split("\t",$damx[$i]);
    chomp $valx;
    if ($vald eq "NA" or $valx eq "NA")
    {
        next
    }
    
    if($vald*$c_damx<$valx*$c_dam)
    {
        my $n11 = $vald;
        my $n12 = $valx;
        my $n1p = $vald+$valx;
        my $n21 = sprintf("%.0f", ($n1p*$c_dam/($c_damx+$c_dam)));
        my $n22 = sprintf("%.0f", ($n1p*$c_damx/($c_damx+$c_dam)));
        my $np1 = $n11+$n21;
        my $npp = $n1p+$n21+$n22;
        
        $left_value = calculateStatistic( n11=>$n11,
                                            n1p=>$n1p,
                                            np1=>$np1,
                                            npp=>$npp);
        $p_v=$left_value;
        my $dost=300;
        if($p_v>0)
        {
            $dost=-log($p_v)/log(10);
        }
        push @true_dost, $dost;
        push @true, join("\t", ($chrx,$stx,$endx,$dost));
    }
    
    
}
#FDR
for(my $i=0;$i<=$#first;$i++)
{
    my ($chrd,$std,$endd,$vald)=split("\t",$first[$i]);
    chomp $vald;
    my ($chrx,$stx,$endx,$valx)=split("\t",$second[$i]);
    chomp $valx;
    if ($vald eq "NA" or $valx eq "NA")
    {
        next
    }
    if($vald*$c_second<$valx*$c_first)
    {
        my $n11 = $vald;
        my $n12 = $valx;
        my $n1p = $vald+$valx;
        my $n21 = sprintf("%.0f", ($n1p*$c_first/($c_first+$c_second)));
        my $n22 = sprintf("%.0f", ($n1p*$c_second/($c_first+$c_second)));
        my $np1 = $n11+$n21;
        my $npp = $n1p+$n21+$n22;
        
        $left_value = calculateStatistic( n11=>$n11,
                                            n1p=>$n1p,
                                            np1=>$np1,
                                            npp=>$npp);
        $p_v=$left_value;

        my $dost=300;
        if($p_v>0)
        {
            $dost=-log($p_v)/log(10);
        }
        push @false_dost, $dost;
        push @false, join("\t", ($chrx,$stx,$endx,$dost));
    }
}
#####

if ($#false_dost != $#false or $#true_dost != $#true)
{
    print "ERROR\n";
    my @stop=<STDIN>;
}

my @comb_dost=();
push @comb_dost, @false_dost;
push @comb_dost, @true_dost;

(@comb_dost)=sort {$a<=>$b} @comb_dost;  #  (по возрастанию)

for(my $d=1;$d<=$#comb_dost;$d++)
{
    $comb_dost[$d] = sprintf("%.14f", $comb_dost[$d]);
    if ($comb_dost[$d]==$comb_dost[$d-1])
    {
        splice @comb_dost, $d,1;
        $d--;
    }
    
}
##############################
my $taget_found=0;
my @fdrs=();

for(my $f=0;$f<=$#comb_dost;$f++)
{
    my $pos=0;
    my $neg=0;
    my $lim=$comb_dost[$f];
    if ($lim<1.3010299956639811952137388947245)    #p>0.05
    {
        next;
    }
    my $p=10**(-$lim);
    print "p=$p\n";
    for(my $q=0;$q<=$#true_dost;$q++)
    {
        if ($true_dost[$q]>=$lim)
        {
            $pos++
        }
        if ($false_dost[$q]>=$lim)
        {
            $neg++
        }
    }
    my $fdr=0;
    if($neg>0)
    {
        $fdr=$neg/($neg+$pos);
    }
    else
    {
        print "\n";
        last
    }
    
    if ($fdr<=$target_fdr and $taget_found==0)
    {
        $taget_found=1;
        $found_prob_level=$p;        
        print OUT "target probability:\t$p\nor -log(10)= $lim\nwith FDR:\t$fdr\n\n";
    }
    push @fdrs, join("\t",($p,$pos,$neg,$fdr));
    
}
if ($taget_found==0)
{
    print "Target FDR=${target_fdr} can't be reached. Try to raise FDR cutoff value.
    \nThis problem usually occurs when experimental noise between two replicates of Dam-X is comparable to actual difference between Dam and Dam-X\nPress Enter to exit\n";
    my $exit=<STDIN>;
    die 
}

print OUT "probability\tTrue_peaks\tFalse_peaks\tFDR\n";
for($z=0;$z<=$#fdrs;$z++)
{
    print OUT "$fdrs[$z]\n"
}

}
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
if (!$manpval)
{
    $limit=-log($found_prob_level)/log(10);
}
else
{
    $limit=-log($manpval)/log(10);
}
if ($manpval)
{
    if ($damx_over_dam_lim!=1)
    {
        open(OUT2, ">${in5}_at_Pvalue_${manpval}_quotient_${damx_over_dam_lim}_peaks.bed");
    }
    else
    {
        open(OUT2, ">${in5}_at_Pvalue_${manpval}_peaks.bed");
    }
    open(OUT3, ">${in5}_at_Pvalue_${manpval}_probability_track.wig");
    open(OUT5, ">${in5}_at_Pvalue_${manpval}_divided_track.wig");
    open(OUT6, ">${in5}_at_Pvalue_${manpval}_full_dataset.txt");
}
else
{
    if ($damx_over_dam_lim!=1)
    {
        open(OUT2, ">${in5}_at_FDR_${target_fdr}_quotient_${damx_over_dam_lim}_peaks.bed");
    }
    else
    {
        open(OUT2, ">${in5}_at_FDR_${target_fdr}_peaks.bed");
    }
    open(OUT3, ">${in5}_at_FDR_${target_fdr}_probability_track.wig");
    open(OUT5, ">${in5}_at_FDR_${target_fdr}_divided_track.wig");
    open(OUT6, ">${in5}_at_FDR_${target_fdr}_probability_full_dataset.txt");    
}
my $r1=int(rand(255));
my $g1=int(rand(255));
my $b1=int(rand(255));

print OUT6 "chrom\tstart\tend\tFisher_probability\tlog2(X/Dam)\n";

print OUT3 'track type=wiggle_0 name="'.$in5.'_DamID_Fisher_log10_probability_track" visibility=full autoScale=off viewLimits=-17:17';
print OUT3 " yLineOnOff=on yLineMark=${limit} color=${r1},${g1},${b1}";
print OUT3 "\n";

print OUT5 'track type=wiggle_0 name="'.$in5.'_DamID_divided_track" visibility=full autoScale=off viewLimits=-10:10';
print OUT5 " yLineOnOff=on yLineMark=$lim_div color=${r1},${g1},${b1}";
print OUT5 "\n";

##########
my $c_damn=0;
my $c_damxn=0;
print "Wait...\n";

for(my $c=0;$c<=$#damn;$c++)
{
    my ($chrd,$std,$endd,$vald)=split("\t",$damn[$c]);
    my ($chrx,$stx,$endx,$valx)=split("\t",$damxn[$c]);
    chomp ($vald,$valx);
    if ($vald eq "NA" or $valx eq "NA")
    {
        next
    }
    chomp ($vald,$valx);
    $c_damn+=$vald;
    $c_damxn+=$valx;
 
}
my $coeff=$c_damn/$c_damxn;

#########

if($#damn != $#damxn)
{
    print "\nDam and DamX averaged datasets have different lenth! Script will stop. Press Enter.";
    <STDIN>;
    exit
}

for($i=0;$i<=$#damn;$i++)
{
    my ($chrd,$std,$endd,$vald)=split("\t",$damn[$i]);
    chomp $vald;
    my ($chrx,$stx,$endx,$valx)=split("\t",$damxn[$i]);
    chomp $valx;
    if ($vald eq "NA" or $valx eq "NA")
    {
        print OUT6 "$chrd\t$std\t$endd\tNA\tNA\n";
        next
    }
    if(($valx+$vald)==0)
    {
        print OUT6 "$chrd\t$std\t$endd\t0\t0\n";
        next
    }
    if($vald<$valx*$coeff) 
    {
        #           Dam     Dam-X
        #    obs    n11      n12 | n1p
        #    exp    n21      n22 | n2p
        #           --------------
        #           np1      np2   npp
        my $n11 = $vald;
        my $n12 = $valx;
        my $n1p = $vald+$valx;
        my $n21 = sprintf("%.0f", ($n1p*$c_damn/($c_damxn+$c_damn)));
        my $n22 = sprintf("%.0f", ($n1p*$c_damxn/($c_damxn+$c_damn)));
        my $np1 = $n11+$n21;
        my $npp = $n1p+$n21+$n22;
        $left_value = calculateStatistic( n11=>$n11,
                                          n1p=>$n1p,
                                          np1=>$np1,
                                          npp=>$npp);
        $p=$left_value;
        
        my $dost=300;
        if ($p>0)
        {
            $dost=-log($p)/log(10)
        }
        my $div="NA";
        if ($vald>0 and $valx>0)
        {
            $div=sprintf("%.3f",log(($valx*$coeff)/$vald)/log(2));
        }
        elsif($valx>0 and $vald==0)
        {
            $div="MAX";
        }
        
        if($dost>$limit)
        {
            my $dam_pl1=$vald+1;
            my $damx_pl1=$valx+1;
            my $div_track=sprintf("%.3f",log(($damx_pl1*$coeff)/$dam_pl1)/log(2));
            print OUT5 "$chrd\t$std\t$endd\t$div_track\n";
        }
        
        print OUT6 "$chrd\t$std\t$endd\t$dost\t$div\n";
        print OUT3 "$chrd\t$std\t$endd\t$dost\n";
        if($dost>$limit and ($div>=$lim_div or $div eq "MAX"))
        {
            print OUT2 "$chrd\t$std\t$endd\n";
        }
    }else       #—читаем достоверности "вниз"
    {
        #           Dam     Dam-X
        #    obs    n11      n12 | n1p
        #    exp    n21      n22 | n2p
        #           --------------
        #           np1      np2   npp
        my $n11 = $valx;
        my $n12 = $vald;
        my $n1p = $vald+$valx;
        my $n21 = sprintf("%.0f", ($n1p*$c_damxn/($c_damxn+$c_damn)));
        my $n22 = sprintf("%.0f", ($n1p*$c_damn/($c_damxn+$c_damn)));
        my $np1 = $n11+$n21;
        my $npp = $n1p+$n21+$n22;
        $left_value = calculateStatistic( n11=>$n11,
                                          n1p=>$n1p,
                                          np1=>$np1,
                                          npp=>$npp);
        $p=$left_value;
        my $dost=-300;
        if ($p>0)
        {
            $dost=log($p)/log(10)
        }
        my $div="NA";
        if ($vald>0 and $valx>0)
        {
            $div=sprintf("%.3f",log(($valx*$coeff)/$vald)/log(2));
        }
        elsif($valx==0 and $vald>0)
        {
            $div="MIN";
        }
               
        if($dost<-$limit)
        {
            my $dam_pl1=$vald+1;
            my $damx_pl1=$valx+1;
            my $div_track=sprintf("%.3f",log(($damx_pl1*$coeff)/$dam_pl1)/log(2));
            print OUT5 "$chrd\t$std\t$endd\t$div_track\n";
        }
        print OUT6 "$chrd\t$std\t$endd\t$dost\t$div\n";
        print OUT3 "$chrd\t$std\t$endd\t$dost\n";        
    }
}
print "\nDone\n";