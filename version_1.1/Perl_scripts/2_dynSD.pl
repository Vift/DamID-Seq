$sd_border=2;
$min_dots_in_win=200;
print 'Enter the first Dam repeat file name (without "_reads_per_GATC_filtered.txt"):'."\n";
$in1=<STDIN>;
chomp $in1;
print "\n";
open(DAT1, "${in1}_reads_per_GATC_filtered.txt") || (print 'no such a file: "'."${in1}_reads_per_GATC_filtered.txt".'"'."\nScript will stop. Press Enter" and my $exit=<STDIN> and die);
print 'Enter the second Dam repeat file name (without "_reads_per_GATC_filtered.txt"):'."\n";
$in2=<STDIN>;
chomp $in2;
print "\n";
open(DAT2, "${in2}_reads_per_GATC_filtered.txt") || (print 'no such a file: "'."${in2}_reads_per_GATC_filtered.txt".'"'."\nScript will stop. Press Enter" and $exit=<STDIN> and die);
print 'If there is the third Dam repeat, please enter it\'s file name (without "_reads_per_GATC_filtered.txt")'."\nOR".' just press "Enter" if there are only two Dam replicas:'."\n";
$in22=<STDIN>;
chomp $in22;
print "\n";
my $repeats=2;
if ($in22)
{
    $repeats=3;
    open(DAT22, "${in22}_reads_per_GATC_filtered.txt") || (print 'no such a file: "'."${in22}_reads_per_GATC_filtered.txt".'"'."\nScript will stop. Press Enter" and $exit=<STDIN> and die);
}

print 'Enter the first Dam-X repeat file name (without "_reads_per_GATC_filtered.txt"):'."\n";
$in3=<STDIN>;
chomp $in3;
print "\n";
open(DAT3, "${in3}_reads_per_GATC_filtered.txt") || (print 'no such a file: "'."${in3}_reads_per_GATC_filtered.txt".'"'."\nScript will stop. Press Enter" and $exit=<STDIN> and die);
print 'Enter the second Dam-X repeat file name (without "_reads_per_GATC_filtered.txt"):'."\n";
$in4=<STDIN>;
chomp $in4;
print "\n";
open(DAT4, "${in4}_reads_per_GATC_filtered.txt") || (print 'no such a file: "'."${in4}_reads_per_GATC_filtered.txt".'"'."\nScript will stop. Press Enter" and $exit=<STDIN> and die);
if ($repeats==3)
{
    
    print 'Enter the third Dam-X repeat file name (without "_reads_per_GATC_filtered.txt"):'."\n";
    $in44=<STDIN>;
    chomp $in44;
    open(DAT44, "${in44}_reads_per_GATC_filtered.txt") || (print 'no such a file: "'."${in44}_reads_per_GATC_filtered.txt".'"'."\nScript will stop. Press Enter" and $exit=<STDIN> and die);
}
print "\n\n";

#################################################################### 2 repeats ################################################################################################


if ($repeats==2)
{

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
        $val1=~s/\015?\012/\n/g;
        $c1+=$val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
        $val2=~s/\015?\012/\n/g;
        $c2+=$val2;
        if($val1+$val2>0)
        {
                push @pre_norm, join("\t",($chr1,$st1,$end1,$val1,$val2));
                $not_null++;
        }else
        {
                
                push @null, join("\t",($chr1,$st1,$end1,$val1));
        }
}

my $mult1=1000000/$c1;
my $mult2=1000000/$c2;
print "not null GATC: $not_null\nWait...\n";

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
my @tummy=sort {$a<=>$b} ($mult1,$mult2);
my $mindiff=$tummy[0];
  
for($lim_up=-$mindiff;$lim_up<=$max*2;)
{
        my @dots_inside=();
        my $di_count=0;
        for($nabor_tochek=0;$nabor_tochek<=$max*2;$nabor_tochek++)
        {
                $lim_up+=$mindiff;
                for (my $d=0;$d<=$#norm;$d++)
                {
                        my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$norm[$d]);
                        $dif=~s/\015?\012/\n/g;
                        if(($val1+$val2)<=$lim_up)
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
        my @line1=();
        for($s=0;$s<=$#dots_inside;$s++)
        {
                my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$dots_inside[$s]);
                $dif=~s/\015?\012/\n/g;
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
                $dif=~s/\015?\012/\n/g;
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
        $val1=~s/\015?\012/\n/g;
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$i]);
        $val2=~s/\015?\012/\n/g;
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
        $val1=~s/\015?\012/\n/g;
        my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
        $val2=~s/\015?\012/\n/g;
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
        $val1=~s/\015?\012/\n/g;
        $c1+=$val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
        $val2=~s/\015?\012/\n/g;
        $c2+=$val2;
        if($val1+$val2>0)
        {
                push @pre_norm, join("\t",($chr1,$st1,$end1,$val1,$val2));
                $not_null++;
        }else
        {
                push @null, join("\t",($chr1,$st1,$end1,$val1));
                
        }
}

$mult1=1000000/$c1;
$mult2=1000000/$c2;
print "not null GATC: $not_null\n";

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
@tummy=sort {$a<=>$b} ($mult1,$mult2);
$mindiff=$tummy[0];
  
for($lim_up=-$mindiff;$lim_up<=$max*2;)
{
        my @dots_inside=();
        my $di_count=0; 
        for($nabor_tochek=0;$nabor_tochek<=$max*2;$nabor_tochek++)
        {
                $lim_up+=$mindiff;
                for (my $d=0;$d<=$#norm;$d++)
                {
                        my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$norm[$d]);
                        $dif=~s/\015?\012/\n/g;
                        if(($val1+$val2)<=$lim_up)
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
        my @line1=();
        for($s=0;$s<=$#dots_inside;$s++)
        {
                my ($chr,$st,$end,$val1,$val2,$dif)=split("\t",$dots_inside[$s]);
                $dif=~s/\015?\012/\n/g;
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
                $dif=~s/\015?\012/\n/g;
                if($dif>($averdata+$three_sd) or $dif<($averdata-$three_sd))
                {
                        $outside++;
                        push @null, join("\t",($chr,$st,$end,"NA"));
                        
                }else
                {
                        if($val1+$val2>0)
                        {
                                $inside++;
                        }
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
        $val1=~s/\015?\012/\n/g;
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$i]);
        $val2=~s/\015?\012/\n/g;
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

}
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
#################################################################### 3 repeats ###########################
else
{


open(OUT, ">${in1}_filtered.txt");
open(OUT2, ">${in2}_filtered.txt");
open(OUT3, ">${in3}_filtered.txt");
open(OUT4, ">${in4}_filtered.txt");
open(OUT22, ">${in22}_filtered.txt");
open(OUT44, ">${in44}_filtered.txt");
open(SUM, ">${in1}+${in2}+${in22}_filtered_combined.txt");
open(STAT, ">${in3}+${in4}+${in44}_vs_${in1}+${in2}+${in22}_stat.txt");
open(SUM2, ">${in3}+${in4}+${in44}_filtered_combined.txt");


@data1=<DAT1>;
@data2=<DAT2>;
@data22=<DAT22>;

if($#data1 != $#data2 or $#data22 != $#data2)
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
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
    $val2=~s/\015?\012/\n/g;
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
my $pearson12=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$data2[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$data22[$i]);
    $val2=~s/\015?\012/\n/g;
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
my $pearson23=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$data1[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$data22[$i]);
    $val2=~s/\015?\012/\n/g;
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
my $pearson13=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

print STAT "Dam (${in1}+${in2}+${in22}) statistics:\n\n";
print STAT "${in1} vs ${in2}:\nCorrelation before filtering:\tR=$pearson12\n";
print STAT "${in2} vs ${in22}:\nCorrelation before filtering:\tR=$pearson23\n";
print STAT "${in1} vs ${in22}:\nCorrelation before filtering:\tR=$pearson13\n\n";
######################## here
my @mass0=();                   
my @mass1=();
my @mass2=();
my @mass22=();
my @mass3=();
my @mass4=();
my @mass44=();

close (DAT1);
close (DAT2);
close (DAT22);

my $c1=0;
my $c2=0;
my $c3=0;
my $not_null=0;
my @all_dots=();
my @pre_norm=();
my @norm=();

my @povt1=();
my @povt2=();
my @povt3=();
my @null=();
my $inside=0;
my $outside=0;
my $outside1=0;
my $outside2=0;
my $outside3=0;
my $passed_summ=0;
my $filt_summ=0;
print "Dam (${in1}+${in2}+${in22}):\n";

for($i=0;$i<=$#data1-3;$i++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$data1[$i]);
        $val1=~s/\015?\012/\n/g;
        $c1+=$val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
        $val2=~s/\015?\012/\n/g;
        $c2+=$val2;
        my ($chr3,$st3,$end3,$val3)=split("\t",$data22[$i]);
        $val3=~s/\015?\012/\n/g;
        $c3+=$val3;
        if($val1+$val2+$val3>0)
        {
                push @pre_norm, join("\t",($chr1,$st1,$end1,$val1,$val2,$val3));
                $not_null++;
        }else
        {
                push @null, join("\t",($chr1,$st1,$end1,0));
        }
}

my $mult1=1000000/$c1;
my $mult2=1000000/$c2;
my $mult3=1000000/$c3;
print "not null GATC: $not_null\nWait...\n";

foreach (@pre_norm)
{
        my ($chr,$st,$end,$val1,$val2,$val3)=split("\t",$_);
        $val1*=$mult1;
        $val2*=$mult2;
        $val3*=$mult3;
        push @all_dots, $val1;
        push @all_dots, $val2;
        push @all_dots, $val3;
        push @norm, join("\t",($chr,$st,$end,$val1,$val2,$val3,($val1-$val2),($val2-$val3),($val3-$val1)));
}

my @dummy=sort {$a<=>$b} @all_dots;
my $max=$dummy[$#dummy];
my @tummy=sort {$a<=>$b} ($mult1,$mult2,$mult3);
my $mindiff=$tummy[0];
for($lim_up=-$mindiff;$lim_up<=$max*3;)
{
        my @dots_inside=();
        my $di_count=0;
        for($nabor_tochek=0;$nabor_tochek<=$max*3;$nabor_tochek++)
        {
            $lim_up+=$mindiff;
            for (my $d=0;$d<=$#norm;$d++)
            {
                my ($chr,$st,$end,$val1,$val2,$val3,$dif12,$dif23,$dif31)=split("\t",$norm[$d]);
                if(($val1+$val2+$val3)<=$lim_up )
                {
                        push @dots_inside, join("\t",($chr,$st,$end,$val1,$val2,$val3,$dif12,$dif23,$dif31));
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
        my @line1=();
        for($s=0;$s<=$#dots_inside;$s++)
        {
            my ($chr,$st,$end,$val1,$val2,$val3,$dif12,$dif23,$dif31)=split("\t",$dots_inside[$s]);
            push @line1, ($dif12,$dif23,$dif31);
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
        ####Outsiders_lookout
        for($s=0;$s<=$#dots_inside;$s++)
        {
            my ($chr,$st,$end,$val1,$val2,$val3,$dif12,$dif23,$dif31)=split("\t",$dots_inside[$s]);
            my $d12=0;
            my $d23=0;
            my $d31=0;
            if  ($dif12>($averdata+$three_sd) or $dif12<($averdata-$three_sd))
            {
                $d12=1;
            }
            if  ($dif23>($averdata+$three_sd) or $dif23<($averdata-$three_sd))
            {
                $d23=1;
            }
            if  ($dif31>($averdata+$three_sd) or $dif31<($averdata-$three_sd))
            {
                $d31=1;
            }
            if ($d12+$d23+$d31>2)
            {
                $outside++;
                $outside1++;
                $outside2++;
                $outside3++;
                push @null, join("\t",($chr,$st,$end,"NA"));
            }else
            {
                $inside++;
                if ($d12+$d31<2)
                {
                    push @povt1, join("\t",($chr,$st,$end,$val1));
                }
                else
                {
                    push @povt1, join("\t",($chr,$st,$end,"NA"));
                    $outside1++;
                }
                if ($d23+$d12<2)
                {
                    push @povt2, join("\t",($chr,$st,$end,$val2));
                }
                else
                {
                    push @povt2, join("\t",($chr,$st,$end,"NA"));
                    $outside2++;
                }
                if ($d23+$d31<2)
                {
                    push @povt3, join("\t",($chr,$st,$end,$val3));
                }
                else
                {
                    push @povt3, join("\t",($chr,$st,$end,"NA"));
                    $outside3++;
                }
            }
        }
}
my $out_perc= sprintf("%.2f", (100*$outside/($inside+$outside)));
my $out1_perc= sprintf("%.2f", (100*$outside1/($inside+$outside)));
my $out2_perc= sprintf("%.2f", (100*$outside2/($inside+$outside)));
my $out3_perc= sprintf("%.2f", (100*$outside3/($inside+$outside)));
my $part_filt_out=$outside1+$outside2+$outside3;
print STAT "Passed filter:\t$inside\nFiltered out:\t$outside\nFiltered portion:\t$out_perc%\n\n";
print "\nPassed filter:\t$inside\nFiltered out in all replicas:\t$outside ($out_perc%)\nFiltered out in at least one replica:\t$part_filt_out\n\n";
print "Wait...\n";



push @povt1,@null;
push @povt2,@null;
push @povt3,@null;

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
@povt3 = sort
   {

   (split("\t",$a))[0] cmp (split("\t",$b))[0]
||
   (split("\t",$a))[1] <=> (split("\t",$b))[1]
   
   }@povt3;
   
if(($#povt1 != $#povt2) or ($#povt1 != $#povt3))
{
        print "\nRepeats lengths appear different. First is $#povt1, the second one $#povt2, the third one $#povt3. Closing.";
        <STDIN>;
        exit
}

for($v=0;$v<=$#povt1;$v++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$v]);
        if ($val1 ne "NA")
        {
             $val1=sprintf("%.0f", $val1/$mult1);
        }
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$v]);
        if ($val2 ne "NA")
        {
             $val2=sprintf("%.0f", $val2/$mult2);
        }        
        my ($chr3,$st3,$end3,$val3)=split("\t",$povt3[$v]);
        if ($val3 ne "NA")
        {
             $val3=sprintf("%.0f", $val3/$mult3);
        }        
        push @mass0, join("\t",($chr1,$st1,$end1));
        push @mass1, $val1;
        push @mass2, $val2;
        push @mass22, $val3;
}

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$i]);
    $val2=~s/\015?\012/\n/g;
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
$pearson12=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$povt2[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$povt3[$i]);
    $val2=~s/\015?\012/\n/g;
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
$pearson23=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$povt3[$i]);
    $val2=~s/\015?\012/\n/g;
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
$pearson13=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

print STAT "GATC filtered in ${in1}:\t$outside1 ($out1_perc%)\nGATC filtered in ${in2}:\t$outside2 ($out2_perc%)\nGATC filtered in ${in22}:\t$outside3 ($out3_perc%)\n";
print STAT "\n${in1} vs ${in2}:\nCorrelation after filtering:\tR=$pearson12\n";
print STAT "${in2} vs ${in22}:\nCorrelation after filtering:\tR=$pearson23\n";
print STAT "${in1} vs ${in22}:\nCorrelation after filtering:\tR=$pearson13\n\n";

#########################################################################################################
#########################################################################################################

@data1=();
@data2=();
@data22=();
@data1=<DAT3>;
@data2=<DAT4>;
@data22=<DAT44>;

if($#data1 != $#data2 or $#data22 != $#data2)
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
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
    $val2=~s/\015?\012/\n/g;
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
$pearson12=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$data2[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$data22[$i]);
    $val2=~s/\015?\012/\n/g;
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
$pearson23=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$data1[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$data22[$i]);
    $val2=~s/\015?\012/\n/g;
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
$pearson13=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);
print STAT "----------------------------------------------------\n\nDam-X (${in3}+${in4}+${in44}) statistics:\n\n";
print STAT "${in3} vs ${in4}:\nCorrelation before filtering:\tR=$pearson12\n";
print STAT "${in4} vs ${in44}:\nCorrelation before filtering:\tR=$pearson23\n";
print STAT "${in3} vs ${in44}:\nCorrelation before filtering:\tR=$pearson13\n\n";
######################## here
$c1=0;
$c2=0;
$c3=0;
$not_null=0;
@all_dots=();
@pre_norm=();
@norm=();

@povt1=();
@povt2=();
@povt3=();
@null=();
$inside=0;
$outside=0;
$outside1=0;
$outside2=0;
$outside3=0;
$passed_summ=0;
$filt_summ=0;
print "Dam-X (${in3}+${in4}+${in44}):\n";
for($i=0;$i<=$#data1-3;$i++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$data1[$i]);
        $val1=~s/\015?\012/\n/g;
        $c1+=$val1;
        my ($chr2,$st2,$end2,$val2)=split("\t",$data2[$i]);
        $val2=~s/\015?\012/\n/g;
        $c2+=$val2;
        my ($chr3,$st3,$end3,$val3)=split("\t",$data22[$i]);
        $val3=~s/\015?\012/\n/g;
        $c3+=$val3;
        if($val1+$val2+$val3>0)
        {
                push @pre_norm, join("\t",($chr1,$st1,$end1,$val1,$val2,$val3));
                $not_null++;
        }else
        {
                push @null, join("\t",($chr1,$st1,$end1,$val1));
        }
}

$mult1=1000000/$c1;
$mult2=1000000/$c2;
$mult3=1000000/$c3;
print "not null GATC: $not_null\nWait...\n";

foreach (@pre_norm)
{
        my ($chr,$st,$end,$val1,$val2,$val3)=split("\t",$_);
        $val1*=$mult1;
        $val2*=$mult2;
        $val3*=$mult3;
        push @all_dots, $val1;
        push @all_dots, $val2;
        push @all_dots, $val3;
        push @norm, join("\t",($chr,$st,$end,$val1,$val2,$val3,($val1-$val2),($val2-$val3),($val3-$val1)));
        
}
@dummy=sort {$a<=>$b} @all_dots;
$max=$dummy[$#dummy];
@tummy=sort {$a<=>$b} ($mult1,$mult2,$mult3);
$mindiff=$tummy[0];
for($lim_up=-$mindiff;$lim_up<=$max*3;)
{
        my @dots_inside=();
        my $di_count=0;
        for($nabor_tochek=0;$nabor_tochek<=$max*3;$nabor_tochek++)
        {
            $lim_up+=$mindiff;
            for (my $d=0;$d<=$#norm;$d++)
            {
                my ($chr,$st,$end,$val1,$val2,$val3,$dif12,$dif23,$dif31)=split("\t",$norm[$d]);
                $dif31=~s/\015?\012/\n/g;
                if(($val1+$val2+$val3)<=$lim_up )
                {
                        push @dots_inside, join("\t",($chr,$st,$end,$val1,$val2,$val3,$dif12,$dif23,$dif31));
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
        my @line1=();
        for($s=0;$s<=$#dots_inside;$s++)
        {
            my ($chr,$st,$end,$val1,$val2,$val3,$dif12,$dif23,$dif31)=split("\t",$dots_inside[$s]);
            push @line1, ($dif12,$dif23,$dif31);
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
        ####Outsiders_lookout
        for($s=0;$s<=$#dots_inside;$s++)
        {
            my ($chr,$st,$end,$val1,$val2,$val3,$dif12,$dif23,$dif31)=split("\t",$dots_inside[$s]);
            my $d12=0;
            my $d23=0;
            my $d31=0;
            if  ($dif12>($averdata+$three_sd) or $dif12<($averdata-$three_sd))
            {
                $d12=1;
            }
            if  ($dif23>($averdata+$three_sd) or $dif23<($averdata-$three_sd))
            {
                $d23=1;
            }
            if  ($dif31>($averdata+$three_sd) or $dif31<($averdata-$three_sd))
            {
                $d31=1;
            }
            if ($d12+$d23+$d31>2)
            {
                $outside++;
                $outside1++;
                $outside2++;
                $outside3++;
                push @null, join("\t",($chr,$st,$end,"NA"));
            }else
            {
                $inside++;
                if ($d12+$d31<2)
                {
                    push @povt1, join("\t",($chr,$st,$end,$val1));
                }
                else
                {
                    push @povt1, join("\t",($chr,$st,$end,"NA"));
                    $outside1++;
                }
                if ($d23+$d12<2)
                {
                    push @povt2, join("\t",($chr,$st,$end,$val2));
                }
                else
                {
                    push @povt2, join("\t",($chr,$st,$end,"NA"));
                    $outside2++;
                }
                if ($d23+$d31<2)
                {
                    push @povt3, join("\t",($chr,$st,$end,$val3));
                }
                else
                {
                    push @povt3, join("\t",($chr,$st,$end,"NA"));
                    $outside3++;
                }
            }
        }
}

$out_perc= sprintf("%.2f", (100*$outside/($inside+$outside)));
$out1_perc= sprintf("%.2f", (100*$outside1/($inside+$outside)));
$out2_perc= sprintf("%.2f", (100*$outside2/($inside+$outside)));
$out3_perc= sprintf("%.2f", (100*$outside3/($inside+$outside)));
$part_filt_out=$outside1+$outside2+$outside3;
print STAT "Passed filter:\t$inside\nFiltered out:\t$outside\nFiltered portion:\t$out_perc%\n\n";
print "\nPassed filter:\t$inside\nFiltered out in all replicas:\t$outside ($out_perc%)\nFiltered out in at least one replica:\t$part_filt_out\n\n";
print "Wait...\n";



push @povt1,@null;
push @povt2,@null;
push @povt3,@null;

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
@povt3 = sort
   {

   (split("\t",$a))[0] cmp (split("\t",$b))[0]
||
   (split("\t",$a))[1] <=> (split("\t",$b))[1]
   
   }@povt3;
   
if(($#povt1 != $#povt2) or ($#povt1 != $#povt3))
{
        print "\nRepeats lengths appear different. First is $#povt1, the second one $#povt2, the third one $#povt3. Closing.";
        <STDIN>;
        exit
}

for($v=0;$v<=$#povt1;$v++)
{
        my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$v]);
        if ($val1 ne "NA")
        {
             $val1=sprintf("%.0f", $val1/$mult1);
        }
        my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$v]);
        if ($val2 ne "NA")
        {
             $val2=sprintf("%.0f", $val2/$mult2);
        }        
        my ($chr3,$st3,$end3,$val3)=split("\t",$povt3[$v]);
        if ($val3 ne "NA")
        {
             $val3=sprintf("%.0f", $val3/$mult3);
        }        

        push @mass3, $val1;
        push @mass4, $val2;
        push @mass44, $val3;
}
@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$povt2[$i]);
    $val2=~s/\015?\012/\n/g;
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
$pearson12=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$povt2[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$povt3[$i]);
    $val2=~s/\015?\012/\n/g;
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
$pearson23=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);

@pmass1=();
@pmass2=();
$pcount=0;
$averval1=0;
$averval2=0;
for(my $i=0;$i<=$#data1-3;$i++)
{
    my ($chr1,$st1,$end1,$val1)=split("\t",$povt1[$i]);
    $val1=~s/\015?\012/\n/g;
    my ($chr2,$st2,$end2,$val2)=split("\t",$povt3[$i]);
    $val2=~s/\015?\012/\n/g;
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
$pearson13=sprintf("%.3f", $cov/($sigma1*$sigma2)**0.5);
print STAT "GATC filtered in ${in3}:\t$outside1 ($out1_perc%)\nGATC filtered in ${in4}:\t$outside2 ($out2_perc%)\nGATC filtered in ${in44}:\t$outside3 ($out3_perc%)\n";
print STAT "\n${in3} vs ${in4}:\nCorrelation after filtering:\tR=$pearson12\n";
print STAT "${in4} vs ${in44}:\nCorrelation after filtering:\tR=$pearson23\n";
print STAT "${in3} vs ${in44}:\nCorrelation after filtering:\tR=$pearson13\n\n";
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

if ($#mass3 != $#mass2 or $#mass1 != $#mass4 or $#mass3 != $#mass1 or $#mass1 != $#mass0 or $#mass1 != $#mass22 or $#mass1 != $#mass44)
{
        open (M0, ">M0.txt");
        foreach(@mass0){print M0 "$_\n"};
        open (M1, ">M1.txt");
        foreach(@mass1){print M1 "$_\n"};
        open (M3, ">M3.txt");
        foreach(@mass3){print M3 "$_\n"};        
        print "$#mass0\t$#mass1\t$#mass2\t$#mass22\t$#mass3\t$#mass4\t$#mass44\n";
        print 'Something wrong! Report Vift@mcb.nsc.ru'."\n";
        my $qqq=<STDIN>;
        exit
}
my $depth_f=0;
my $depth_s=0;
my $depth1=0;
my $depth2=0;
my $depth22=0;
my $depth3=0;
my $depth4=0;
my $depth44=0;
for(my $d=0;$d<=$#mass1;$d++)
{
    $depth1+=$mass1[$d] if $mass1[$d] ne "NA";
    $depth2+=$mass2[$d] if $mass2[$d] ne "NA";
    $depth22+=$mass22[$d] if $mass22[$d] ne "NA";
    $depth3+=$mass3[$d] if $mass3[$d] ne "NA";
    $depth4+=$mass4[$d] if $mass4[$d] ne "NA";
    $depth44+=$mass44[$d] if $mass44[$d] ne "NA";
}

for (my $z=0;$z<=$#mass1;$z++)
{
    my ($av1,$av2)=0;
    my ($rep_passed1,$rep_passed2)=0;
    if(
           ($mass1[$z] eq "NA" and $mass2[$z] eq "NA")
        or ($mass2[$z] eq "NA" and $mass22[$z] eq "NA")
        or ($mass22[$z] eq "NA" and $mass1[$z] eq "NA")
        or ($mass3[$z] eq "NA" and $mass4[$z] eq "NA")
        or ($mass4[$z] eq "NA" and $mass44[$z] eq "NA")
        or ($mass44[$z] eq "NA" and $mass3[$z] eq "NA")
       )
    {
            $filt_summ++;
            print OUT "$mass0[$z]\tNA\t$rep_passed1\n";
            print OUT2 "$mass0[$z]\tNA\t$rep_passed1\n";
            print OUT22 "$mass0[$z]\tNA\t$rep_passed1\n";
            print OUT3 "$mass0[$z]\tNA\t$rep_passed2\n";
            print OUT4 "$mass0[$z]\tNA\t$rep_passed2\n";
            print OUT44 "$mass0[$z]\tNA\t$rep_passed2\n";
            print SUM "$mass0[$z]\tNA\t$rep_passed1\n";
            print SUM2 "$mass0[$z]\tNA\t$rep_passed2\n";
    }
    else
    {
            $passed_summ++;
            print OUT "$mass0[$z]\t$mass1[$z]\n";
            print OUT2 "$mass0[$z]\t$mass2[$z]\n";
            print OUT22 "$mass0[$z]\t$mass22[$z]\n";
            print OUT3 "$mass0[$z]\t$mass3[$z]\n";
            print OUT4 "$mass0[$z]\t$mass4[$z]\n";
            print OUT44 "$mass0[$z]\t$mass44[$z]\n";
            if
            (($mass1[$z] ne "NA") and ($mass2[$z] ne "NA") and ($mass22[$z] ne "NA") )
            {
                $av1=($mass1[$z]+$mass2[$z]+$mass22[$z]);
                $depth_f=$depth1+$depth2+$depth22;
            }elsif
            (($mass1[$z] ne "NA") and ($mass2[$z] ne "NA") and ($mass22[$z] eq "NA") )
            {
                $av1=($mass1[$z]+$mass2[$z]);
                $depth_f=$depth1+$depth2;
            }elsif
            (($mass1[$z] ne "NA") and ($mass2[$z] eq "NA") and ($mass22[$z] ne "NA") )
            {
                $av1=($mass1[$z]+$mass22[$z]);
                $depth_f=$depth1+$depth22;
            }elsif
            (($mass1[$z] eq "NA") and ($mass2[$z] ne "NA") and ($mass22[$z] ne "NA") )
            {
                $av1=($mass2[$z]+$mass22[$z]);
                $depth_f=$depth2+$depth22;
            }
            if
            (($mass3[$z] ne "NA") and ($mass4[$z] ne "NA") and ($mass44[$z] ne "NA") )
            {
                $av2=($mass3[$z]+$mass4[$z]+$mass44[$z]);
                $depth_s=$depth3+$depth4+$depth44;
            }elsif
            (($mass3[$z] ne "NA") and ($mass4[$z] ne "NA") and ($mass44[$z] eq "NA") )
            {
                $av2=($mass3[$z]+$mass4[$z]);
                $depth_s=$depth3+$depth4;
            }elsif
            (($mass3[$z] ne "NA") and ($mass4[$z] eq "NA") and ($mass44[$z] ne "NA") )
            {
                $av2=($mass3[$z]+$mass44[$z]);
                $depth_s=$depth3+$depth44;
            }elsif
            (($mass3[$z] eq "NA") and ($mass4[$z] ne "NA") and ($mass44[$z] ne "NA") )
            {
                $av2=($mass4[$z]+$mass44[$z]);
                $depth_s=$depth4+$depth44;
            }
            print SUM "$mass0[$z]\t$av1\t$depth_f\n";
            print SUM2 "$mass0[$z]\t$av2\t$depth_s\n";                
    }
    
}
$out_perc= sprintf("%.2f", (100*$filt_summ/($passed_summ+$filt_summ)));
print STAT "\n-------------------------------------\n\n";
print STAT "Final statistics:\nPassed filter:\t$passed_summ\nFiltered out:\t$filt_summ\nFiltered portion:\t$out_perc%\n\n";

}
print "Done\n";
