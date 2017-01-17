$sd_border=2;
$min_dots_in_win=200;
print 'Enter the first Dam repeat file name (without "_reads_per_GATC_filtered.txt"):'."\n";
$in1=<STDIN>;
chomp $in1;

print 'Enter the second Dam repeat file name (without "_reads_per_GATC_filtered.txt"):'."\n";
$in2=<STDIN>;
chomp $in2;

print 'Enter the first Dam-X repeat file name (without "_reads_per_GATC_filtered.txt"):'."\n";
$in3=<STDIN>;
chomp $in3;

print 'Enter the second Dam-X repeat file name (without "_reads_per_GATC_filtered.txt"):'."\n";
$in4=<STDIN>;
chomp $in4;

open(DAT1, "${in1}_reads_per_GATC_filtered.txt") || die "no suchafile";
open(DAT2, "${in2}_reads_per_GATC_filtered.txt") || die "no suchafile";
open(DAT3, "${in3}_reads_per_GATC_filtered.txt") || die "no suchafile";
open(DAT4, "${in4}_reads_per_GATC_filtered.txt") || die "no suchafile";

open(OUT, ">${in1}_filtered.txt");
open(OUT2, ">${in2}_filtered.txt");
open(OUT3, ">${in3}_filtered.txt");
open(OUT4, ">${in4}_filtered.txt");
open(SUM, ">${in1}+${in2}_filtered_combined.txt");
open(STAT, ">${in3}+${in4}_vs_${in1}+${in2}_stat.txt");
open(SUM2, ">${in3}+${in4}_filtered_combined.txt");



@data1=<DAT1>;
@data2=<DAT2>;

if($#data1 != $#data2)
{
    print "\ndatasets have different length! Script will stop. Press Enter.";
    my $st=<STDIN>;
    exit
}

my @pmass1=();
my @pmass2=();
my $pcount=0;
my $averval1=0;
my $averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$data1[$i]);
        chomp $val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
        chomp $val2;
        if ($val1==0 and $val2==0)
        {
                next
        }
        push @pmass1, $val1;
        push @pmass2, $val2;
        $pcount++;
        $averval1+=$val1;
        $averval2+=$val2;
        
}
$averval1/=$pcount;
$averval2/=$pcount;

$cov=0;
$sigma1=0;
$sigma2=0;

for($q=0;$q<=$#pmass1;$q++)
{
        $cov+=($pmass1[$q]-$averval1)*($pmass2[$q]-$averval2);
        $sigma1+=($pmass1[$q]-$averval1)**2;
        $sigma2+=($pmass2[$q]-$averval2)**2;
        
}

my $pearson=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);
my $psq=sprintf("%.3f", $pearson**2);
print STAT "${in1}+${in2} statistics:\nCorrelation before filtering:\tR=$pearson\n";

my @mass0=();
my @mass1=();
my @mass2=();
my @mass3=();
my @mass4=();

close (DAT1);
close (DAT2);

my $c1=0;
my $c2=0;
my $not_null=0;
my @all_dots=();
my @pre_norm=();
my @norm=();

my @povt1=();
my @povt2=();
my @null=();
my $inside=0;
my $outside=0;
my $passed_summ=0;
my $filt_summ=0;
print "${in1}+${in2}:\n";
for($i=0;$i<=$#data1-3;$i++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$data1[$i]);
        chomp $val1;
        $c1+=$val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
        chomp $val2;
        $c2+=$val2;
        if($val1+$val2>0)
        {
                push @pre_norm, join("\t",($chr1,$st1,$end1,$val1,$val2));
                $not_null++;
        }else
        {
                $inside++;
                push @null, join("\t",($chr1,$st1,$end1,$val1));
        }
}

my $mult1=1000000/$c1;
my $mult2=1000000/$c2;
print "not null pairs $not_null\nWait...\n";

foreach (@pre_norm)
{
        my ($chr,$st,$end,$val1,$val2)=split("\t",$_);
        $val1*=$mult1;
        $val2*=$mult2;
        push @all_dots, $val1;
        push @all_dots, $val2;
        push @norm, join("\t",($chr,$st,$end,$val1,$val2,($val1-$val2)));
}
(@dummy)=sort {$a<=>$b} @all_dots;
my $max=$dummy[$#dummy];
my $norm_size=$#norm;

  
for($lim_up=-1;$lim_up<=$max*2;)
{
        my @dots_inside=();
        my $di_count=0;
        #print "___border___\n";
        #print LINE "$lim_up\t";       
        for($nabor_tochek=0;$nabor_tochek<=$norm_size;$nabor_tochek++)
        {
                $lim_up++;
                $lim_up=sprintf("%.0f", $lim_up);
                $perc=sprintf("%.2f",(100*$lim_up/($max*2)));
                
                #print "$perc%\n";
                for (my $d=0;$d<=$#norm;$d++)
                {
                        my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$norm[$d]);
                        chomp $dif;
                        if(sprintf("%.0f", ($val1+$val2))==$lim_up)
                        {
                                push @dots_inside, join("\t",($chr,$st,$end,$val1,$val2,$dif));
                                splice @norm, $d,1;
                                $d--;
                                $di_count++;
                        }
                }
                if($di_count>=$min_dots_in_win or $lim_up>$max*2)   
                {
                        last
                }
        }
        #print "dots_inside\t$di_count\nnorm_size\t$norm_size\n\n";

        my @line1=();
        for($s=0;$s<=$#dots_inside;$s++)
        {
                my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$dots_inside[$s]);
                chomp $dif;
                push @line1, $dif;
        }
        #################### #SD counter, needs @line1 with data
        my $averdata=0;
        foreach(@line1)
        {
            $averdata+=$_/($#line1+1)
        }
        my $summdata=0;
        foreach(@line1)
        {
            $summdata+=($_-$averdata)*($_-$averdata);
        }
        if($#line1>0)
        {$sd=sqrt($summdata/($#line1));}
        else{$sd=0}
        ####################
        $three_sd=$sd*$sd_border;
        #print LINE "$three_sd\n";
        ####ищем выпаденцев
        for($s=0;$s<=$#dots_inside;$s++)
        {
                my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$dots_inside[$s]);
                chomp $dif;
                if($dif>($averdata+$three_sd) or $dif<($averdata-$three_sd))
                {
                        $outside++;
                        push @null, join("\t",($chr,$st,$end,"NA"));
                        
                }else
                {
                        $inside++;
                        push @povt1, join("\t",($chr,$st,$end,$val1));
                        push @povt2, join("\t",($chr,$st,$end,$val2));
                }
        }
        
}

my $out_perc= sprintf("%.2f", (100*$outside/($inside+$outside)));
print STAT "Passed filter:\t$inside\nFiltered out:\t$outside\nFiltered portion:\t$out_perc%\n";
print "${in1}+${in2} statistics:\nPassed filter:\t$inside\nFiltered out:\t$outside\nFiltered portion:\t$out_perc%\n\n";
print "Wait...\n";

my $cc1=0;
my $cc2=0;


push @povt1,@null;
push @povt2,@null;


@povt1 = sort
   {

   (split("\t",$a))[0] cmp (split("\t",$b))[0]
||
   (split("\t",$a))[1] <=> (split("\t",$b))[1]
   
   }@povt1;
@povt2 = sort
   {

   (split("\t",$a))[0] cmp (split("\t",$b))[0]
||
   (split("\t",$a))[1] <=> (split("\t",$b))[1]
   
   }@povt2;
if($#povt1 != $#povt2)
{
        print "\nRepeats lengths appear different. First is $#povt1, the second one $#povt2. Closing.";
        <STDIN>;
        exit
}
for($k=0;$k<=$#povt1;$k++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$k]);
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$k]);
        $cc1+=$val1;
        $cc2+=$val2;
        
}


for($v=0;$v<=$#povt1;$v++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$v]);
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$v]);
        
        if ($val1 ne "NA" and $val2 ne "NA")
        {
                $val1=sprintf("%.0f", $val1/$mult1);
                $val2=sprintf("%.0f", $val2/$mult2);
        }
        else
        {
                $val1="NA";
                $val2="NA";
        }
        push @mass0, join("\t",($chr1,$st1,$end1));
        push @mass1, $val1;
        push @mass2, $val2;
}
@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#povt1;$i++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$i]);
        chomp $val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$i]);
        chomp $val2;
        if ($val1==0 and $val2==0)
        {
                next
        }
        push @pmass1, $val1;
        push @pmass2, $val2;
        $pcount++;
        $averval1+=$val1;
        $averval2+=$val2;
        
}
$averval1/=$pcount;
$averval2/=$pcount;

$cov=0;
$sigma1=0;
$sigma2=0;

for($q=0;$q<=$#pmass1;$q++)
{
        $cov+=($pmass1[$q]-$averval1)*($pmass2[$q]-$averval2);
        $sigma1+=($pmass1[$q]-$averval1)**2;
        $sigma2+=($pmass2[$q]-$averval2)**2;
        
}

$pearson=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);
$psq=sprintf("%.3f", $pearson**2);
print STAT "Correlation after filtering:\tR=$pearson\n\n";
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
print "${in3}+${in4}:\n";
@data1=();
@data2=();

@data1=<DAT3>;
@data2=<DAT4>;

if($#data1 != $#data2)
{
    print "\ndatasets have different length! Script will stop. Press Enter.";
    my $st=<STDIN>;
    exit
}

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$data1[$i]);
        chomp $val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
        chomp $val2;
        if ($val1==0 and $val2==0)
        {
                next
        }
        push @pmass1, $val1;
        push @pmass2, $val2;
        $pcount++;
        $averval1+=$val1;
        $averval2+=$val2;
        
}
$averval1/=$pcount;
$averval2/=$pcount;

$cov=0;
$sigma1=0;
$sigma2=0;

for($q=0;$q<=$#pmass1;$q++)
{
        $cov+=($pmass1[$q]-$averval1)*($pmass2[$q]-$averval2);
        $sigma1+=($pmass1[$q]-$averval1)**2;
        $sigma2+=($pmass2[$q]-$averval2)**2;
        
}

$pearson=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);
$psq=sprintf("%.3f", $pearson**2);
print STAT "${in3}+${in4} statistics:\nCorrelation before filtering:\tR=$pearson\n";


$c1=0;
$c2=0;
$not_null=0;
@all_dots=();
@pre_norm=();
@norm=();

@povt1=();
@povt2=();
@null=();
$inside=0;
$outside=0;

for($i=0;$i<=$#data1-3;$i++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$data1[$i]);
        chomp $val1;
        $c1+=$val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
        chomp $val2;
        $c2+=$val2;
        if($val1+$val2>0)
        {
                push @pre_norm, join("\t",($chr1,$st1,$end1,$val1,$val2));
                $not_null++;
        }else
        {
                push @null, join("\t",($chr1,$st1,$end1,$val1));
                $inside++;
        }
}

$mult1=1000000/$c1;
$mult2=1000000/$c2;
print "not null pairs $not_null\n";

foreach (@pre_norm)
{
        my ($chr,$st,$end,$val1,$val2)=split("\t",$_);
        $val1*=$mult1;
        $val2*=$mult2;
        push @all_dots, $val1;
        push @all_dots, $val2;
        push @norm, join("\t",($chr,$st,$end,$val1,$val2,($val1-$val2)));
}
(@dummy)=sort {$a<=>$b} @all_dots;
$max=$dummy[$#dummy];
$norm_size=$#norm;

  
for($lim_up=-1;$lim_up<=$max*2;)
{
        my @dots_inside=();
        my $di_count=0;
        #print "___border___\n";
        #print LINE "$lim_up\t";       
        for($nabor_tochek=0;$nabor_tochek<=$norm_size;$nabor_tochek++)
        {
                $lim_up++;
                $lim_up=sprintf("%.0f", $lim_up);
                $perc=sprintf("%.2f",(100*$lim_up/($max*2)));
                
                #print "$perc%\n";
                for (my $d=0;$d<=$#norm;$d++)
                {
                        my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$norm[$d]);
                        chomp $dif;
                        if(sprintf("%.0f", ($val1+$val2))==$lim_up)
                        {
                                push @dots_inside, join("\t",($chr,$st,$end,$val1,$val2,$dif));
                                splice @norm, $d,1;
                                $d--;
                                $di_count++;
                        }
                }
                if($di_count>=$min_dots_in_win or $lim_up>$max*2)   
                {
                        last
                }
        }
        #print "dots_inside\t$di_count\nnorm_size\t$norm_size\n\n";

        my @line1=();
        for($s=0;$s<=$#dots_inside;$s++)
        {
                my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$dots_inside[$s]);
                chomp $dif;
                push @line1, $dif;
        }
        #################### #SD counter, needs @line1 with data
        my $averdata=0;
        foreach(@line1)
        {
            $averdata+=$_/($#line1+1)
        }
        my $summdata=0;
        foreach(@line1)
        {
            $summdata+=($_-$averdata)*($_-$averdata);
        }
        if($#line1>0)
        {$sd=sqrt($summdata/($#line1));}
        else{$sd=0}
        ####################
        $three_sd=$sd*$sd_border;
        #print LINE "$three_sd\n";
        ####ищем выпаденцев
        for($s=0;$s<=$#dots_inside;$s++)
        {
                my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$dots_inside[$s]);
                chomp $dif;
                if($dif>($averdata+$three_sd) or $dif<($averdata-$three_sd))
                {
                        $outside++;
                        push @null, join("\t",($chr,$st,$end,"NA"));
                        
                }else
                {
                        $inside++;
                        push @povt1, join("\t",($chr,$st,$end,$val1));
                        push @povt2, join("\t",($chr,$st,$end,$val2));
                }
        }
        
}

$out_perc= sprintf("%.2f", (100*$outside/($inside+$outside)));
print STAT "Passed filter:\t$inside\nFiltered out:\t$outside\nFiltered portion:\t$out_perc%\n";
print "${in3}+${in4} statistics:\nPassed filter:\t$inside\nFiltered out:\t$outside\nFiltered portion:\t$out_perc%\n\n";
print "Wait...\n";

$cc1=0;
$cc2=0;


push @povt1,@null;
push @povt2,@null;


@povt1 = sort
   {

   (split("\t",$a))[0] cmp (split("\t",$b))[0]
||
   (split("\t",$a))[1] <=> (split("\t",$b))[1]
   
   }@povt1;
@povt2 = sort
   {

   (split("\t",$a))[0] cmp (split("\t",$b))[0]
||
   (split("\t",$a))[1] <=> (split("\t",$b))[1]
   
   }@povt2;
if($#povt1 != $#povt2)
{
        print "\nRepeats lengths appear different. First is $#povt1, the second one $#povt2. Closing.";
        <STDIN>;
        exit
}
for($k=0;$k<=$#povt1;$k++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$k]);
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$k]);
        $cc1+=$val1;
        $cc2+=$val2;
        
}


for($v=0;$v<=$#povt1;$v++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$v]);
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$v]);
        
        if ($val1 ne "NA" and $val2 ne "NA")
        {
                $val1=sprintf("%.0f", $val1/$mult1);
                $val2=sprintf("%.0f", $val2/$mult2);
        }
        else
        {
                $val1="NA";
                $val2="NA";
        }
        push @mass3, $val1;
        push @mass4, $val2;
}
@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#povt1;$i++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$i]);
        chomp $val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$i]);
        chomp $val2;
        if ($val1==0 and $val2==0)
        {
                next
        }
        push @pmass1, $val1;
        push @pmass2, $val2;
        $pcount++;
        $averval1+=$val1;
        $averval2+=$val2;
        
}
$averval1/=$pcount;
$averval2/=$pcount;

$cov=0;
$sigma1=0;
$sigma2=0;

for($q=0;$q<=$#pmass1;$q++)
{
        $cov+=($pmass1[$q]-$averval1)*($pmass2[$q]-$averval2);
        $sigma1+=($pmass1[$q]-$averval1)**2;
        $sigma2+=($pmass2[$q]-$averval2)**2;
        
}

$pearson=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);
$psq=sprintf("%.3f", $pearson**2);
print STAT "Correlation after filtering:\tR=$pearson\n\n";

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

if ($#mass3 != $#mass2 or $#mass1 != $#mass4 or $#mass3 != $#mass1 or $#mass1 != $#mass0)
{
        print "Something wrong!\n";
        my $qqq=<STDIN>;
}

for (my $z=0;$z<=$#mass1;$z++)
{
        if ($mass1[$z] eq "NA" or $mass2[$z] eq "NA" or $mass3[$z] eq "NA" or $mass4[$z] eq "NA")
        {
                $filt_summ++;
                print OUT "$mass0[$z]\tNA\n";
                print OUT2 "$mass0[$z]\tNA\n";
                print OUT3 "$mass0[$z]\tNA\n";
                print OUT4 "$mass0[$z]\tNA\n";
                print SUM "$mass0[$z]\tNA\n";
                print SUM2 "$mass0[$z]\tNA\n";
        }
        else
        {
                $passed_summ++;
                print OUT "$mass0[$z]\t$mass1[$z]\n";
                print OUT2 "$mass0[$z]\t$mass2[$z]\n";
                print OUT3 "$mass0[$z]\t$mass3[$z]\n";
                print OUT4 "$mass0[$z]\t$mass4[$z]\n";
                my $av1=($mass1[$z]+$mass2[$z]);
                my $av2=($mass3[$z]+$mass4[$z]);
                print SUM "$mass0[$z]\t$av1\n";
                print SUM2 "$mass0[$z]\t$av2\n";                
        }
        
}
$out_perc= sprintf("%.2f", (100*$filt_summ/($passed_summ+$filt_summ)));
print STAT "Final statistics:\nPassed filter:\t$passed_summ\nFiltered out:\t$filt_summ\nFiltered portion:\t$out_perc%\n\n";





























