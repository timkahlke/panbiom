#!/usr/bin/perl

use strict;
use Data::Dumper;
use Getopt::Std;
use warnings;

eval{
    _main();
};

if($@){
    print $@;
}



sub _main{
    our %opts;
    getopts("t:g:i:o:x:am:",\%opts);

    my $in = $opts{'i'};
    my $out = $opts{'o'};
    my $rel_ab = $opts{'a'};
    my $horizontal = $opts{'x'};
    my $min = $opts{'m'}?$opts{'m'}:0;
    my ($samples,$reps) = $opts{'g'}?_get_map($opts{'g'}):(0,0);
    my $rep_out = $opts{'t'};

    if(!$in||!$out){
        _usage();
    }

    my ($data,$colh,$rowh) = _get_data($in,$horizontal,$samples);

    if(!($rel_ab)){
        $data = _get_rel($data);
    }

    my $core;
    if($reps){
        if(!($rep_out)){
            print STDOUT "\n[ERROR] No replicate outlier count specified with option -t!!!\n\n";
            _usage();
        }
        $core = _get_rep_core($data,$min,$colh,$reps,$rep_out);
    }
    else{
        $core = _get_core($data,$min,$colh);
    }

    print STDOUT "\n\nCore of given data and parameters: ".scalar(@$core)."\n\n";
}


sub _get_rep_core{
    my ($data,$min,$colh,$reps,$ol) = @_;
    my $core = [];
    for(my $i = 0;$i<scalar(@$data);$i++){
        my $row = $data->[$i];
        my $cgc = 0;
        my $ca = [];
        foreach my $rg(keys(%$reps)){
            my $ro = 0;
            foreach my $s(keys(%{$reps->{$rg}})){
                if($row->[$colh->{$s}]<$min){
                    $ro++;
                }
                else{
                    push @$ca,$row->[$colh->{$s}];
                }
            }
            if($ro<=$ol){
                warn Dumper $ca;
                $cgc++;
            }
        }
        if($cgc == scalar(keys(%$reps))){
            push @$core,$i;
        }
    }
    die Dumper $core;
    return $core;
}


sub _get_map{
    my $file = shift;
    open(my $ih,"<",$file) or die "Failed to open $file";
    my $samples = {};
    my $reps = {};
    while(my $line=<$ih>){
        $line=~s/[\r\n]//g;
        my @s = split(/\t/,$line);
        $samples->{$s[0]} = 1;
        if($s[1]){
            $samples->{$s[0]} = $s[1];
            $reps->{$s[1]}->{$s[0]} = 0;
        }
    }
    close($ih);
    return ($samples,$reps);
}

sub _get_core{
    my ($data,$min,$colh) = @_;
    my $core = [];
    for(my $i = 0;$i<scalar(@$data);$i++){
        if($colh){
            my $nd = [];
            foreach(values(%$colh)){
                push @$nd, $data->[$i]->[$_];
            }
            next unless _check($nd,$min);
        }
        else{
            next unless _check($data->[$i],$min);
        }
        push @$core,$i;
    }
    return $core;
}


sub _check{
    my ($r,$min) = @_;
    foreach(@$r){
        if($_<$min){ 
            return 0;
        }
    }
    return 1;
}


sub _get_col_sums{
    my $d = shift;
    my $cs = [];
    my $mi = scalar(@{$d->[1]});
    foreach my $r(@$d){
        for(my $i=0;$i<$mi;$i++){
            if($cs->[$i]){
                $cs->[$i]+=$r->[$i];
            }
            else{
                $cs->[$i] = $r->[$i];
            }
        }
    }
    return $cs;
}


sub _get_rel{
    my $d = shift;
    my @rel_data = ();
    my $col_sums = _get_col_sums($d);
    foreach my $r(@$d){
        my $nl = [];
        for(my $i=0;$i<scalar(@$r);$i++){
            $nl->[$i] = (($r->[$i]*100)/$col_sums->[$i]);     
        }
        push @rel_data,$nl;
    }
    return \@rel_data;
}


sub _rel{
    my ($v,$s) = @_;
    return (($v*100)/$s);
}


sub _get_sum{
    my $a = shift;
    my $s = 0;
    foreach(@$a){
        $s+=$_;
    }
    return $s;
}




# create data structure with samples as rows
sub _get_data{
    my ($file,$hor,$samples) = @_;

    open(my $ih,"<",$file) or die "Failed to open input file $file";
    my $rh = {};
    my $ch = {};
    my $data = [];
    
    while(my $line=<$ih>){
        if($.==1){
            $ch = _get_rheader($line,$samples);
            next;
        }
        my $s=_split_line($line);
        $rh->{shift(@$s)} = scalar(keys(%$rh));

        if(!($hor)){
            push @$data, $s;
        }
        else{
            for(my $i=0;$i<scalar(@$s);$i++){
                if(ref($data->[$i])){
                    push @{$data->[$i]},$s->[$i];
                }
                else{
                    $data->[$i] = [$s->[$i]];
                }
            }
        }
    }
    close($ih);
    return ($data,$ch,$rh);
}


sub _split_line{
    my $line = shift;
    $line=~s/[\r\n]//g;
    my @s = split(/\t/,$line);
    return \@s;
}


sub _get_rheader{
    my ($line,$samples) = @_;
    my $rheader = {};
    my $s=_split_line($line);
    for(my $i=1;$i<scalar(@$s);$i++){
        if($samples){
            next unless $samples->{$s->[$i]};
        }
        $rheader->{$s->[$i]} = $i-1;
    }   
    return $rheader;
}

sub _usage{
    print STDOUT "\nScript to estimate the core microbiome of a set of samples in a tab-separated tetx file\n";
    print STDOUT "Parameters:\n";
    print STDOUT "i : input text file (tab-separated) with first column and first row headers\n";
    print STDOUT "o : output file\n";
    print STDOUT "a : (optional) if set count data is assumed to be relative abundance\n";
    print STDOUT "x : (optional) if set samples are assumed to be rows not columns\n";
    print STDOUT "m : minimum an OTU/Species has to have in percent per sample\n";
    print STDOUT "g : mapping file with groups of samples defining which samples to take into account (first column) and which are replicates (second column, optional)\n";
    print STDOUT "t : (optional) number of outliers allowed per replicate group\n";
    print STDOUT "\n\n";
    exit;
}





