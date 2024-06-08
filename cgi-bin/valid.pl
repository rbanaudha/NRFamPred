#!/usr/bin/perl

############################### Header Information ##############################
require 'cgi.perl';
use CGI;;
$query = new CGI;
&ReadParse;
print &PrintHeader;

################################ Reads Inoput Data ##############################
$atom = $query->param('atom');
$file = $query->param('file');
$svm_th = $query->param('svm_th');

#################Validation Of Input Sequence Data (file upload) ###################################
if($file ne '' && $atom eq '')
{
    $file=~m/^.*(\\|\/)(.*)/; 
    while(<$file>) 
    {
	$seqfi .= $_;
    }
}
elsif($atom ne '' && $file eq ''){

    $seqfi="$atom";
}

##############ACTUAL PROCESS BEGINS FROM HERE#######################
$infut_file = "/webservers/cgi-bin/nrfampred";
$ran= int(rand 10000);
#$dir = "/webservers/cgidocs/mkumar/temp/nrfam$ran";
$dir = "/webservers/cgidocs/mkumar/temp/Ravindra/NRfamPred/nrfam$ran";
system "mkdir $dir";
system "chmod 777 $dir";
#$nam = 'input.'.'fasta';
open(FP1,">$dir/input.fasta");
print FP1 "$seqfi\n";
#print "$seqfi\n";
close FP1;

system "/usr/bin/perl $infut_file/fasta.pl $dir/input.fasta >$dir/twoline";
system "/bin/grep -c '>' $dir/twoline >$dir/total_seq";
system "/bin/grep '>' $dir/twoline |/usr/bin/cut -d '|' -f3 |/usr/bin/cut -d ' ' -f1 >$dir/protein_id";
$total_seq=`head -1 $dir/total_seq`;chomp($total_seq);
if($total_seq <= 25)
{
    system "/usr/bin/perl $infut_file/dipep.pl $dir/twoline >$dir/dipep_comp";
    system "/bin/sed -e 's/^/+1 /' $dir/dipep_comp >$dir/final";
    system "/usr/local/bin/svm_classify $dir/final $infut_file/Models/level1_model $dir/svm_score_leve1 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/final $infut_file/Models/level2_model_NR0 $dir/svm_score_leve2_NR0 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/final $infut_file/Models/level2_model_NR1 $dir/svm_score_leve2_NR1 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/final $infut_file/Models/level2_model_NR2 $dir/svm_score_leve2_NR2 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/final $infut_file/Models/level2_model_NR3 $dir/svm_score_leve2_NR3 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/final $infut_file/Models/level2_model_NR4 $dir/svm_score_leve2_NR4 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/final $infut_file/Models/level2_model_NR5 $dir/svm_score_leve2_NR5 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/final $infut_file/Models/level2_model_NR6 $dir/svm_score_leve2_NR6 >/dev/null";
    system "/usr/bin/paste $dir/protein_id $dir/final $dir/svm_score_leve1 $dir/svm_score_leve2_NR0 $dir/svm_score_leve2_NR1 $dir/svm_score_leve2_NR2 $dir/svm_score_leve2_NR3 $dir/svm_score_leve2_NR4 $dir/svm_score_leve2_NR5 $dir/svm_score_leve2_NR6 >$dir/final_svm";
    system "/usr/bin/tr '\t' '#' <$dir/final_svm >$dir/final_pred";
    #open(RESULT,">>$dir/Result") or die "$!";
    #print RESULT "Protein ID\tPrediction\n";
    #close RESULT;
    
    open(FINAL_PRED,"$dir/final_pred") or die "$!";
    while($pred=<FINAL_PRED>)
    {
	chomp($pred);
	@svm=split(/#/,$pred);
	#print "$svm[2]\n";
	if($svm[2] < $svm_th)
	{
	    open(RESULT,">>$dir/Result") or die "$!";
	    print RESULT "$svm[0]\tNon-NR family\n";
	    close RESULT;
	}
	else
	{
	    if(($svm[3] > $svm[4])&&($svm[3] > $svm[5])&&($svm[3] > $svm[6])&&($svm[3] > $svm[7])&&($svm[3] > $svm[8])&&($svm[3] > $svm[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$svm[0]\tNR0 sub-family\n";
		close RESULT;
	    }
	    if(($svm[4] > $svm[3])&&($svm[4] > $svm[5])&&($svm[4] > $svm[6])&&($svm[4] > $svm[7])&&($svm[4] > $svm[8])&&($svm[4] > $svm[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$svm[0]\tNR1 sub-family\n";
		close RESULT;
	    }
	    if(($svm[5] > $svm[3])&&($svm[5] > $svm[4])&&($svm[5] > $svm[6])&&($svm[5] > $svm[7])&&($svm[5] > $svm[8])&&($svm[5] > $svm[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$svm[0]\tNR2 sub-family\n";
		close RESULT;
	    }
	    if(($svm[6] > $svm[3])&&($svm[6] > $svm[4])&&($svm[6] > $svm[5])&&($svm[6] > $svm[7])&&($svm[6] > $svm[8])&&($svm[6] > $svm[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$svm[0]\tNR3 sub-family\n";
		close RESULT;
	    }
	    if(($svm[7] > $svm[3])&&($svm[7] > $svm[4])&&($svm[7] > $svm[5])&&($svm[7] > $svm[6])&&($svm[7] > $svm[8])&&($svm[7] > $svm[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$svm[0]\tNR4 sub-family\n";
		close RESULT;
	    }
	    if(($svm[8] > $svm[3])&&($svm[8] > $svm[4])&&($svm[8] > $svm[5])&&($svm[8] > $svm[6])&&($svm[8] > $svm[7])&&($svm[8] > $svm[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$svm[0]\tNR5 sub-family\n";
		close RESULT;
	    }
	    if(($svm[9] > $svm[3])&&($svm[9] > $svm[4])&&($svm[9] > $svm[5])&&($svm[9] > $svm[6])&&($svm[9] > $svm[7])&&($svm[9] > $svm[8]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$svm[0]\tNR6 sub-family\n";
		close RESULT;
	    }
	}
    }
}

if($total_seq > 25)
{
    system "/usr/bin/perl $infut_file/2line_to_singleline.pl $dir/twoline |grep '>' >$dir/singleline";
    #system "/bin/grep '>' $dir/singleline > $dir/ttt";
    system "/usr/bin/head -25 $dir/singleline >$dir/top_25_seq";
    system "/usr/bin/perl $infut_file/singleline_to_2line.pl $dir/top_25_seq >$dir/top25_seq_twoline";
    system "/bin/grep '>' $dir/top25_seq_twoline |/usr/bin/cut -d '|' -f3 |/usr/bin/cut -d ' ' -f1 >$dir/top_25_protein_id";
    system "/usr/bin/perl $infut_file/dipep.pl $dir/top25_seq_twoline >$dir/top_25_dipep";
    system "/bin/sed -e 's/^/+1 /' $dir/top_25_dipep >$dir/top_25_final";
    system "/usr/local/bin/svm_classify $dir/top_25_final $infut_file/Models/level1_model $dir/top_25_svm_score_leve1 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/top_25_final $infut_file/Models/level2_model_NR0 $dir/top_25_svm_score_leve2_NR0 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/top_25_final $infut_file/Models/level2_model_NR1 $dir/top_25_svm_score_leve2_NR1 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/top_25_final $infut_file/Models/level2_model_NR2 $dir/top_25_svm_score_leve2_NR2 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/top_25_final $infut_file/Models/level2_model_NR3 $dir/top_25_svm_score_leve2_NR3 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/top_25_final $infut_file/Models/level2_model_NR4 $dir/top_25_svm_score_leve2_NR4 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/top_25_final $infut_file/Models/level2_model_NR5 $dir/top_25_svm_score_leve2_NR5 >/dev/null";
    system "/usr/local/bin/svm_classify $dir/top_25_final $infut_file/Models/level2_model_NR6 $dir/top_25_svm_score_leve2_NR6 >/dev/null";
    system "/usr/bin/paste $dir/top_25_protein_id $dir/top_25_final $dir/top_25_svm_score_leve1 $dir/top_25_svm_score_leve2_NR0 $dir/top_25_svm_score_leve2_NR1 $dir/top_25_svm_score_leve2_NR2 $dir/top_25_svm_score_leve2_NR3 $dir/top_25_svm_score_leve2_NR4 $dir/top_25_svm_score_leve2_NR5 $dir/top_25_svm_score_leve2_NR6 >$dir/top_25_final_svm";
    system "/usr/bin/tr '\t' '#' <$dir/top_25_final_svm >$dir/top_25_final_pred";

    open(TOP_25_FINAL_PRED,"$dir/top_25_final_pred") or die "$!";
    while($toppred=<TOP_25_FINAL_PRED>)
    {
	chomp($toppred);
	@val1=split(/#/,$toppred);
	if($val1[2] < $svm_th)
	{
	    open(RESULT,">>$dir/Result") or die "$!";
	    print RESULT "$val1[0]\tNon-NR family\n";
	    close RESULT;
	}
	else
	{
	    if(($val1[3] > $val1[4])&&($val1[3] > $val1[5])&&($val1[3] > $val1[6])&&($val1[3] > $val1[7])&&($val1[3] > $val1[8])&&($val1[3] > $val1[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$val1[0]\tNR0 sub-family\n";
		close RESULT;
	    }
	    if(($val1[4] > $val1[3])&&($val1[4] > $val1[5])&&($val1[4] > $val1[6])&&($val1[4] > $val1[7])&&($val1[4] > $val1[8])&&($val1[4] > $val1[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$val1[0]\tNR1 sub-family\n";
		close RESULT;
	    }
	    if(($val1[5] > $val1[3])&&($val1[5] > $val1[4])&&($val1[5] > $val1[6])&&($val1[5] > $val1[7])&&($val1[5] > $val1[8])&&($val1[5] > $val1[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$val1[0]\tNR2 sub-family\n";
		close RESULT;
	    }
	    if(($val1[6] > $val1[3])&&($val1[6] > $val1[4])&&($val1[6] > $val1[5])&&($val1[6] > $val1[7])&&($val1[6] > $val1[8])&&($val1[6] > $val1[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$val1[0]\tNR3 sub-family\n";
		close RESULT;
	    }
	    if(($val1[7] > $val1[3])&&($val1[7] > $val1[4])&&($val1[7] > $val1[5])&&($val1[7] > $val1[6])&&($val1[7] > $val1[8])&&($val1[7] > $val1[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$val1[0]\tNR4 sub-family\n";
		close RESULT;
	    }
	    if(($val1[8] > $val1[3])&&($val1[8] > $val1[4])&&($val1[8] > $val1[5])&&($val1[8] > $val1[6])&&($val1[8] > $val1[7])&&($val1[8] > $val1[9]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$val1[0]\tNR5 sub-family\n";
		close RESULT;
	    }
	    if(($val1[9] > $val1[3])&&($val1[9] > $val1[4])&&($val1[9] > $val1[5])&&($val1[9] > $val1[6])&&($val1[9] > $val1[7])&&($val1[9] > $val1[8]))
	    {
		open(RESULT,">>$dir/Result") or die "$!";
		print RESULT "$val1[0]\tNR6 sub-family\n";
		close RESULT;
	    }
	}
    }

}
print  "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
print  "<html><HEAD>\n";
print  "<TITLE>NRfamPred::Prediction Result</TITLE>\n";
print  "<META NAME=\"description\" CONTENT=\"NRfamPred, University of Delhi South Campus, INDIA\">\n";
print  "</HEAD><body bgcolor=\"\#FFFFE0\">\n";
print  "<h2 ALIGN = \"CENTER\"> NRfamPred Prediction Result</h2>\n";
print  "<HR ALIGN =\"CENTER\"> </HR>\n";
print  "<p align=\"center\"><font size=4 color=black><b>The submitted protein/proteins belongs to <font color='red'></p>";
print "<table border='1' width='400' align='center'><tr><th>Protein ID</th><th>Prediction</th></tr>";
open(PREDICTION,"$dir/Result") or die "$!";
while($pre=<PREDICTION>)
{
    chomp($pre);
    @pred=split(/\t/,$pre);
    print "<tr align='center'><td>$pred[0]</td><td>$pred[1]</td></tr>";
}
print "</table>";
print "</font></b></font></p>\n";
print  "<p align=\"center\"><font size=3 color=black><b>Thanks for using NRfamPred Prediction Server</b></font></p>\n";
print  "<p align=\"center\"><font size=3 color=black><b>If you have any problem or suggestions please contact <a href='mailto:manish@south.du.ac.in'>Dr. Manish Kumar</a></b></font>. Please mention your job number in any communication.</p></br>\n";
print  "<p ALIGN=\"CENTER\"><b>Your job number is <font color=\"red\">$ran</b></font></p>\n";
print  "</body>\n";
print  "</html>\n";
system "chmod 000 $dir";
