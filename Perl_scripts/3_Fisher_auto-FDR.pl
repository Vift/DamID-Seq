$slow_hardware=1;           #Set 0 if you beleive in your's hardware computing power
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
open(DAM, "${in1}_filtered_combined.txt")  || open(DAM, "${in1}_filtered_averaged.txt")  || die "no suchafile ${in1}_filtered_combined.txt";
open(DAMN, "${in1}_filtered_combined.txt")  || open(DAMN, "${in1}_filtered_averaged.txt")  || die "no suchafile";
print "Enter the combined Dam-X datafile name (without _filtered_combined.txt):\n";

$in2=<STDIN>;
chomp $in2;
open(DAMX, "${in2}_filtered_combined.txt")  || open(DAMX, "${in2}_filtered_averaged.txt")  ||  die "no suchafile ${in2}_filtered_combined.txt";
open(DAMXN, "${in2}_filtered_combined.txt")  || open(DAMXN, "${in2}_filtered_averaged.txt")  || die "no suchafile";

my $found_prob_level=0;
@damn=<DAMN>;
@damxn=<DAMXN>;
@dam=<DAM>;
@damx=<DAMX>;

############################################################################################################auto-pval
if (!$manpval)
{
    
    my $rep=2; 
    if($in2=~/([A-z0-9_\.-]+)(\+)([A-z0-9_\.-]+)(\+)([A-z0-9_\.-]+)/)
    {
        $rep=3;
        $in3=$1;
        open(F, "${in3}_filtered.txt")  || die "no suchafile ${in3}_filtered.txt";
        $in4=$3;
        open(S, "${in4}_filtered.txt")  || die "no suchafile ${in4}_filtered.txt";
        $in44=$5;
        open(T, "${in44}_filtered.txt")  || die "no suchafile ${in44}_filtered.txt";
        @third=<T>;
    }
    elsif($in2=~/([A-z0-9_\.-]+)(\+)([A-z0-9_\.-]+)/)
    {
        $in3=$1;
        open(F, "${in3}_filtered.txt")  || die "no suchafile ${in3}_filtered.txt";
        $in4=$3;
        open(S, "${in4}_filtered.txt")  || die "no suchafile ${in4}_filtered.txt";
    }
    else
    {
        die "Something wrong with input names...";
    }
    print "Calculating FDR:\n";
    open(OUT, ">${in5}_ROC.txt");
    $previous_fh = select(OUT);	
    $| = 1;
    select($previous_fh);   
    
    @first=<F>;
    @second=<S>;
    
    if ($#dam!=$#damx or $#damx!=$#first or $#first!=$#second)
    {
        print "Input files have different lengths! Programm will stop...\n";
        my $stop=<STDIN>;
        exit;
    }
    if($rep==3 and $#first!=$#third)
    {
        print "Input files have different lengths! Programm will stop...\n";
        my $stop=<STDIN>;
        exit;
    }
    ##########  depth of sequencing
    my $c_dam=0;
    my $c_damx=0;
    my $c_first=0;
    my $c_second=0;
    my $c_third=0 if $rep==3;
    for(my $c=0;$c<=$#dam;$c++)
    {
        my ($chrd,$std,$endd,$vald)=split("\t",$dam[$c]);
        my ($chrx,$stx,$endx,$valx)=split("\t",$damx[$c]);
        my ($chrf,$stf,$endf,$valf)=split("\t",$first[$c]);
        my ($chrs,$sts,$ends,$vals)=split("\t",$second[$c]);
        my ($chrt,$stt,$endt,$valt)=split("\t",$third[$c]) if $rep==3;
        chomp ($vald,$valx,$valf,$vals,$valt);
        if ($rep==2 and ($vald eq "NA" or $valx eq "NA"))
        {
            next
        }        
        $c_dam+=$vald if $vald ne "NA";
        $c_damx+=$valx if $valx ne "NA";
        $c_first+=$valf if $valf ne "NA";
        $c_second+=$vals if $vals ne "NA";
        $c_third+=$valt if ($valt ne "NA" and $rep==3);
    }
    #########
    my @true=();
    my @false=();
    my @true_dost=();
    my @false_dost=();
    #####
    my $last_timer=-1;
    # TRUE
    print "\nreal distribution:\n";
    for(my $i=0;$i<=$#dam;$i++)
    {
        my $timer=int(10*$i/($#dam));
        if ($timer>$last_timer)
        {
            print 10*$timer."% done\n";
            $last_timer=$timer;
        }
        my ($chrd,$std,$endd,$vald,$depth_d)=split("\t",$dam[$i]);
        my ($chrx,$stx,$endx,$valx,$depth_x)=split("\t",$damx[$i]);
        chomp ($valx,$vald,$depth_d,$depth_x);
        if ($vald eq "NA" or $valx eq "NA")
        {
            next
        }
        if($rep==2)
        {
            $depth_d=$c_dam;
            $depth_x=$c_damx;
        }
        if($vald*$depth_x<$valx*$depth_d)
        {
            my $n11 = sprintf("%.0f", $vald);
            my $n12 = sprintf("%.0f", $valx);
            my $n1p = sprintf("%.0f", $vald+$valx);
            my $n21 = sprintf("%.0f", ($n1p*$depth_d/($depth_x+$depth_d)));
            my $n22 = sprintf("%.0f", ($n1p*$depth_x/($depth_x+$depth_d)));
            my $p_v=&fishint($n11,$n12,$n21,$n22);
            my $dost=300;
            $p_v=1 if $p_v>1;
            if($p_v>0)
            {
                $dost=-log($p_v)/log(10);
            }
            push @true_dost, $dost;
            push @true, join("\t", ($chrx,$stx,$endx,$dost));
        }
    }
    #FDR
    $last_timer=-1;
    print "\nfalse distribution:\n";
    for(my $i=0;$i<=$#first;$i++)
    {
        my $timer=int(20*$i/($#first));
        if ($timer>$last_timer)
        {
            print 5*$timer."% done\n";
            $last_timer=$timer;
        }
        my $depth_1=$c_first;
        my $depth_2=$c_second;
        my $depth_3=$c_third;
        my ($chrp1,$stp1,$endp1,$valp1)=split("\t",$first[$i]);
        my ($chrp2,$stp2,$endp2,$valp2)=split("\t",$second[$i]);
        my ($chrp3,$stp3,$endp3,$valp3)=split("\t",$third[$i]) if $rep==3;
        chomp ($valp1,$valp2,$valp3);
        my $dost=0;
        if (($valp1 eq "NA" or $valp2 eq "NA" or $valp3 eq "NA") and $rep==3)     # if there is one "NA"
        {
            if ($valp1 ne "NA")
            {
                if($valp3 ne "NA")
                {
                    $valp2=$valp3;
                    $depth_2=$depth_3;
                }
            }
            elsif ($valp1 eq "NA")
            {
                $valp1=$valp3;
                $depth_1=$depth_3;
            }
            my $n11 = sprintf("%.0f", $valp1);
            my $n12 = sprintf("%.0f", $valp2);
            my $n1p = sprintf("%.0f", $valp1+$valp2);
            my $n21 = sprintf("%.0f", ($n1p*$depth_1/($depth_1+$depth_2)));
            my $n22 = sprintf("%.0f", ($n1p*$depth_2/($depth_1+$depth_2)));
            my $p_v=&fishint($n11,$n12,$n21,$n22);
            $dost=300;
            $p_v=1 if $p_v>1;
            if($p_v>0)
            {
                $dost=-log($p_v)/log(10);
            }
        }
        elsif($rep==3)
        {
            my $n11 = sprintf("%.0f", $valp1);
            my $n12 = sprintf("%.0f", $valp2);
            my $n13 = sprintf("%.0f", $valp3);
            my $n1p = sprintf("%.0f", $valp1+$valp2+$valp3);
            my $n21 = sprintf("%.0f", ($n1p*$depth_1/($depth_1+$depth_2+$depth_3)));
            my $n22 = sprintf("%.0f", ($n1p*$depth_2/($depth_1+$depth_2+$depth_3)));
            my $n23 = sprintf("%.0f", ($n1p*$depth_3/($depth_1+$depth_2+$depth_3)));
            my $p_v=&fishint($n11,$n12,$n13,$n21,$n22,$n23);
            $dost=300;
            if($p_v>0)
            {
                $dost=-log($p_v)/log(10);
            }
        }
        else
        {
            my $n11 = sprintf("%.0f", $valp1);
            my $n12 = sprintf("%.0f", $valp2);
            my $n1p = sprintf("%.0f", $valp1+$valp2);
            my $n21 = sprintf("%.0f", ($n1p*$depth_1/($depth_1+$depth_2)));
            my $n22 = sprintf("%.0f", ($n1p*$depth_2/($depth_1+$depth_2)));
            my $p_v=&fishint($n11,$n12,$n21,$n22);
            $dost=300;
            if($p_v>0)
            {
                $dost=-log($p_v)/log(10);
            }
        }
        push @false_dost, $dost;
        push @false, join("\t", ($chrp1,$stp1,$endp1,$dost));
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
    my @fdrs_crash=();
    print "\nCalculating ROC:\n";
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
        for(my $q=0;$q<=$#true_dost;$q++)
        {
            if ($true_dost[$q]>=$lim)
            {
                $pos++
            }
            else
            {
                splice @true_dost, $q,1;
                $q--
            }
        }
        for(my $q=0;$q<=$#false_dost;$q++)
        {
            if ($false_dost[$q]>=$lim)
            {
                $neg++
            }
            else
            {
                splice @false_dost, $q,1;
                $q--
            }            
        }
        my $fdr=0;
        if($neg>0)
        {
            $fdr=$neg/($neg+$pos);
            push @fdrs_crash, join ("\t", ($fdr,$pos));
            print "p=$p\tFDR=$fdr\n";
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
        @fdrs_crash=sort {$a<=>$b} @fdrs_crash;
        my ($fd,$posit)=split("\t",$fdrs_crash[0]);
        print "Target FDR=${target_fdr} can't be reached.\nThe best FDR observed is ${fd} with ${posit} Dam-X peaks.\nYou may set this FDR cutoff value next time.
        \nThis problem usually occurs when an experimental noise between the replicas of Dam-X is comparable to actual difference between Dam and Dam-X.
        \nThe other possible reason is that some samples were mixed up on some preveious stages of the analisys.
        \nPress Enter to exit\n";
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
    open(OUT6, ">${in5}_at_Pvalue_${manpval}_probability_full_dataset.txt");
    open(OUTDIV, ">${in5}_at_Pvalue_${manpval}_full_divided_track.wig"); 
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
    open(OUT5, ">${in5}_at_FDR_${target_fdr}_significant_divided_track.wig");
    open(OUT6, ">${in5}_at_FDR_${target_fdr}_probability_full_dataset.txt");
    open(OUTDIV, ">${in5}_at_FDR_${target_fdr}_full_divided_track.wig"); 
}
my $r1=int(rand(255));
my $g1=int(rand(255));
my $b1=int(rand(255));

print OUT6 "chrom\tstart\tend\tFisher_probability\tlog2(X/Dam)\n";

print OUT3 'track type=wiggle_0 name="'.$in5.'_DamID_Fisher_log10_probability_track" visibility=full autoScale=off viewLimits=-17:17';
print OUT3 " yLineOnOff=on yLineMark=${limit} color=${r1},${g1},${b1}";
print OUT3 "\n";

print OUT5 'track type=wiggle_0 name="'.$in5.'_DamID_significant_divided_track" visibility=full autoScale=off viewLimits=-5:5';
print OUT5 " yLineOnOff=on yLineMark=$lim_div color=${r1},${g1},${b1}";
print OUT5 "\n";

print OUTDIV 'track type=wiggle_0 name="'.$in5.'_DamID_full_divided_track" visibility=full autoScale=off viewLimits=-5:5';
print OUTDIV " yLineOnOff=on yLineMark=$lim_div color=${r1},${g1},${b1}";
print OUTDIV "\n";

##########
my $c_damn=0;
my $c_damxn=0;
my $coeff=1;
print "Looking for peaks and making the profiles...\n";
my $rep=2;
$rep=3 if($in2=~/([A-z0-9_\.-]+)(\+)([A-z0-9_\.-]+)(\+)([A-z0-9_\.-]+)/);
if($rep==2)
{
    for(my $c=0;$c<=$#damn;$c++)
    {
        my ($chrd,$std,$endd,$vald,$depth_d)=split("\t",$damn[$c]);
        my ($chrx,$stx,$endx,$valx,$depth_x)=split("\t",$damxn[$c]);
        chomp ($vald,$valx,$depth_d,$depth_x);
        if ($vald eq "NA" or $valx eq "NA")
        {
            next
        }
        $c_damn+=$vald;
        $c_damxn+=$valx;
    }
    $coeff=$c_damn/$c_damxn;    ##
}
#########
if($#damn != $#damxn)
{
    print "\nDam and DamX averaged datasets have different lenth! Script will stop. Press Enter.";
    <STDIN>;
    exit
}

for($i=0;$i<=$#damn;$i++)
{
    my ($chrd,$std,$endd,$vald,$depth_d)=split("\t",$damn[$i]);
    my ($chrx,$stx,$endx,$valx,$depth_x)=split("\t",$damxn[$i]);
    chomp ($vald,$valx,$depth_d,$depth_x);
    if($rep==3)
    {
        $c_damn=$depth_d;
        $c_damxn=$depth_x;
    }
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
    my $dam_pls1=$vald+1;
    my $damx_pls1=$valx+1;
    my $fdiv_track=sprintf("%.3f",log(($damx_pls1*$c_damn)/($dam_pls1*$c_damxn))/log(2));
    print OUTDIV "$chrd\t$std\t$endd\t$fdiv_track\n";
    if($vald*$c_damxn<$valx*$c_damn)    # Dam-X > Dam
    {
        #           Dam     Dam-X
        #    obs    n11      n12 | n1p
        #    exp    n21      n22 | n2p
        #           --------------
        #           np1      np2   npp
        my $n11 = sprintf("%.0f", $vald);
        my $n12 = sprintf("%.0f", $valx);
        my $n1p = sprintf("%.0f", $vald+$valx);
        my $n21 = sprintf("%.0f", ($n1p*$c_damn/($c_damxn+$c_damn)));
        my $n22 = sprintf("%.0f", ($n1p*$c_damxn/($c_damxn+$c_damn)));
        my $p=&fishint($n11,$n12,$n21,$n22);
        my $dost=300;
        $p=1 if $p>1;
        if ($p>0)
        {
            $dost=-log($p)/log(10)
        }
        my $div="NA";
        if ($vald>0 and $valx>0)
        {
            $div=sprintf("%.3f",log(($valx*$c_damn)/($vald*$c_damxn))/log(2));
        }
        elsif($valx>0 and $vald==0)
        {
            $div="MAX";
        }
        
        if($dost>$limit)
        {
            my $dam_pl1=$vald+1;
            my $damx_pl1=$valx+1;
            my $div_track=sprintf("%.3f",log(($damx_pl1*$c_damn)/($dam_pl1*$c_damxn))/log(2));
            print OUT5 "$chrd\t$std\t$endd\t$div_track\n";
        }
        
        print OUT6 "$chrd\t$std\t$endd\t$dost\t$div\n";
        print OUT3 "$chrd\t$std\t$endd\t$dost\n";
        if($dost>$limit and ($div>=$lim_div or $div eq "MAX"))
        {
            print OUT2 "$chrd\t$std\t$endd\n";
        }
    }
    else       # Dam > Dam-X
    {
        #           Dam     Dam-X
        #    obs    n11      n12 | n1p
        #    exp    n21      n22 | n2p
        #           --------------
        #           np1      np2   npp
        my $n11 = sprintf("%.0f", $vald);
        my $n12 = sprintf("%.0f", $valx);
        my $n1p = sprintf("%.0f", $vald+$valx);
        my $n21 = sprintf("%.0f", ($n1p*$c_damn/($c_damxn+$c_damn)));
        my $n22 = sprintf("%.0f", ($n1p*$c_damxn/($c_damxn+$c_damn)));
        my $p=&fishint($n11,$n12,$n21,$n22);
        $p=1 if $p>1;
        my $dost=-300;
        if ($p>0)
        {
            $dost=log($p)/log(10)
        }
        my $div="NA";
        if ($vald>0 and $valx>0)
        {
            $div=sprintf("%.3f",log(($valx*$c_damn)/($vald*$c_damxn))/log(2));
        }
        elsif($valx==0 and $vald>0)
        {
            $div="MIN";
        }
               
        if($dost<-$limit)
        {
            my $dam_pl1=$vald+1;
            my $damx_pl1=$valx+1;
            my $div_track=sprintf("%.3f",log(($damx_pl1*$c_damn)/($dam_pl1*$c_damxn))/log(2));
            print OUT5 "$chrd\t$std\t$endd\t$div_track\n";
        }
        print OUT6 "$chrd\t$std\t$endd\t$dost\t$div\n";
        print OUT3 "$chrd\t$std\t$endd\t$dost\n";        
    }
}
print "\nDone\n";

sub fisher
{
    my @inp=@_;
    if($#inp==11)
    {
        my @num=sort{$b<=>$a} (@inp[7..11]);
        my @denom=sort{$b<=>$a} (@inp[0..6]);
        my $logf=0;
        for(my $n=($num[0]+1);$n<=$denom[0];$n++)
        {
            $logf-=log($n)
        }
        if ($num[1]>0)
        {
            if ($num[1]>$denom[1])
            {
                for(my $n=($denom[1]+1);$n<=$num[1];$n++)
                {
                    $logf+=log($n)
                }
            }
            if ($num[2]>0)
            {
                if ($num[2]>$denom[2])
                {
                    for(my $n=($denom[2]+1);$n<=$num[2];$n++)
                    {
                        $logf+=log($n)
                    }
                }
                if ($num[3]>0)
                {
                    if ($num[3]>$denom[3])
                    {
                        for(my $n=($denom[3]+1);$n<=$num[3];$n++)
                        {
                            $logf+=log($n)
                        }
                    }
                    if ($num[4]>0)
                    {
                        if ($num[4]>$denom[4])
                        {
                            for(my $n=($denom[4]+1);$n<=$num[4];$n++)
                            {
                                $logf+=log($n)
                            }
                        }
                    }
                }
            }
        }
        if($denom[5]>0)
        {
            for(my $n=1;$n<=$denom[5];$n++)
            {
                $logf-=log($n)
            }
        }
        if($denom[6]>0)
        {
            for(my $n=1;$n<=$denom[6];$n++)
            {
                $logf-=log($n)
            }
        }    
        return ($logf)
    }
    elsif($#inp==8)
    {
        my @num=sort{$b<=>$a} (@inp[5..8]);
        my @denom=sort{$b<=>$a} (@inp[0..4]);
        my $logf=0;
        for(my $n=($num[0]+1);$n<=$denom[0];$n++)
        {
            $logf-=log($n)
        }
        if ($num[1]>0)
        {
            if ($num[1]>$denom[1])
            {
                for(my $n=($denom[1]+1);$n<=$num[1];$n++)
                {
                    $logf+=log($n)
                }
            }
            if ($num[2]>0)
            {
                if ($num[2]>$denom[2])
                {
                    for(my $n=($denom[2]+1);$n<=$num[2];$n++)
                    {
                        $logf+=log($n)
                    }
                }
                if ($num[3]>0)
                {
                    if ($num[3]>$denom[3])
                    {
                        for(my $n=($denom[3]+1);$n<=$num[3];$n++)
                        {
                            $logf+=log($n)
                        }
                    }

                }
            }
        }
        if($denom[4]>0)
        {
            for(my $n=1;$n<=$denom[4];$n++)
            {
                $logf-=log($n)
            }
        }
        return ($logf)
    }
    else
    {
        return "ERROR"
    }
}
sub fishint
{
    my @inp=@_;
    if($#inp==5)
    {
        my ($n11,$n12,$n13,$n21,$n22,$n23)=@inp;
        chomp $n23;
        unless($n11>-1){next}
        $n11=sprintf("%.0f",$n11); #a
        $n12=sprintf("%.0f",$n12); #b
        $n13=sprintf("%.0f",$n13); #c
        $n21=sprintf("%.0f",$n21); #d
        $n22=sprintf("%.0f",$n22); #e
        $n23=sprintf("%.0f",$n23); #f
        my $n1p = $n11+$n12+$n13;  #R1
        my $n2p = $n21+$n22+$n23;  #R2
        my $np1 = $n11+$n21;       #C1
        my $np2 = $n12+$n22;       #C2
        my $np3 = $n13+$n23;       #C3
        my $npp = $n1p+$n2p;       #N
        if($npp>3000 and $slow_hardware==1)
        {
            my $optimiser=10;
            $optimiser=100 if $npp>30000;
            $n11=sprintf("%.0f",$n11/$optimiser); #a
            $n12=sprintf("%.0f",$n12/$optimiser); #b
            $n13=sprintf("%.0f",$n13/$optimiser); #c
            $n21=sprintf("%.0f",$n21/$optimiser); #d
            $n22=sprintf("%.0f",$n22/$optimiser); #e
            $n23=sprintf("%.0f",$n23/$optimiser); #f
            $n1p = $n11+$n12+$n13;  #R1
            $n2p = $n21+$n22+$n23;  #R2
            $np1 = $n11+$n21;       #C1
            $np2 = $n12+$n22;       #C2
            $np3 = $n13+$n23;       #C3
            $npp = $n1p+$n2p;       #N
        }
        my @num=sort{$b<=>$a} ($n1p,$n2p,$np1,$np2,$np3);
        my @denom=sort{$b<=>$a} ($n11,$n12,$n13,$n21,$n22,$n23,$npp);
        my @integralpval=();
        my $l0pval=&fisher(@denom,@num);
        my $last_row_p=0;
        my $firstfound=0;
        my $last_col_p=0;
        my $ee0=0;
        my $cc0=0;
        my $bb00=0;
        my $cc00=0;
        my $dd00=0;
        for(my $aa=0;$aa<=$n1p and $aa<=$np1;$aa++)
        {
            my $row_first_table_found=0;
            for(my $bb=(0>($n1p-$aa-$np3) ? 0 : ($n1p-$aa-$np3)) ;$bb<=($n1p-$aa) and $bb<=$np2;$bb++)
            {
                my $cc=$n1p-$aa-$bb;
                my $dd=$np1-$aa;
                my $ee=$np2-$bb;
                my $ff=$np3-$cc;
                if($row_first_table_found==0)
                {
                    $row_first_table_found=1;
                    if($firstfound==0)
                    {
                        $last_row_p=&fisher($aa,$bb,$cc,$dd,$ee,$ff,$npp,$n1p,$n2p,$np1,$np2,$np3);
                        $firstfound=1;
                        $bb00=$bb;
                        $cc00=$cc;
                        $dd00=$dd;
                    }
                    else
                    {
                        if($bb00!=$bb)
                        {
                            $last_row_p=$last_row_p+log($bb00)+log($dd00)-log($aa)-log($ee);
                        }
                        else
                        {
                            $last_row_p=$last_row_p+log($cc00)+log($dd00)-log($aa)-log($ff);
                        }
                        $bb00=$bb;
                        $cc00=$cc;
                        $dd00=$dd;
                    }
                    push @integralpval, $last_row_p if sprintf("%.13f",$last_row_p)<=sprintf("%.13f",$l0pval);
                    $last_col_p=$last_row_p;
                }
                else
                {
                    $last_col_p=$last_col_p+log($ee0)+log($cc0)-log($bb)-log($ff);
                    push @integralpval,$last_col_p if sprintf("%.13f",$last_col_p)<=sprintf("%.13f",$l0pval);
                }
                $ee0=$ee;
                $cc0=$cc;
            }
        }
        my $intp=0;
        foreach(@integralpval)
        {
            $intp+=exp($_);
        }
        return $intp;
    }
    elsif($#inp==3)
    {
        my ($n11,$n12,$n21,$n22)=@inp;
        chomp $n22;
        unless($n11>-1){next}
        $n11=sprintf("%.0f",$n11); #a
        $n12=sprintf("%.0f",$n12); #b
        $n21=sprintf("%.0f",$n21); #d
        $n22=sprintf("%.0f",$n22); #e
        my $n1p = $n11+$n12;       #R1
        my $n2p = $n21+$n22;       #R2
        my $np1 = $n11+$n21;       #C1
        my $np2 = $n12+$n22;       #C2
        my $npp = $n1p+$n2p;       #N
        my @num=sort{$b<=>$a} ($n1p,$n2p,$np1,$np2);
        my @denom=sort{$b<=>$a} ($n11,$n12,$n21,$n22,$npp);
        my @integralpval=();
        my $l0pval=&fisher(@denom,@num);
        my $last_p=0;
        my $firstfound=0;
        my $bb0=0;
        my $dd0=0;
        my $tmp=0;
        for(my $aa=(0>($n11-$n22) ? 0 : ($n11-$n22));$aa<=$n1p and $aa<=$np1;$aa++)
        {
            my $bb=$n1p-$aa;
            my $dd=$np1-$aa;
            my $ee=$np2-$bb;
            if($firstfound==0)
            {
                $last_p=&fisher($aa,$bb,$dd,$ee,$npp,$n1p,$n2p,$np1,$np2);
                $firstfound=1;
                $bb0=$bb;
                $dd0=$dd;
                push @integralpval, $last_p if sprintf("%.13f",$last_p)<=sprintf("%.13f",$l0pval);
            }
            else
            {
                $last_p=$last_p+log($bb0)+log($dd0)-log($aa)-log($ee);
                $bb0=$bb;
                $dd0=$dd;
                push @integralpval, $last_p if sprintf("%.13f",$last_p)<=sprintf("%.13f",$l0pval);
            }
        }
        my $intp=0;
        foreach(@integralpval)
        {
            $intp+=exp($_);
        }
        return $intp;
    }
    else
    {
        return "ERROR_in_Fisher_inputs"
    }
}