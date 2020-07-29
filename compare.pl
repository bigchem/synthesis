use strict;
use Cwd qw(cwd);

my $fileIn = shift;
my $fileRes = shift;
my $top = shift;
my $isomeric = shift;
my $largest = shift;
my $multi = shift;
my $augmentations = shift; 
my $beams = shift; 
my $skipcanonical = shift;
my $donotconvert = shift;
my $single = shift;

$single = defined $single ? $single: 0;

$isomeric = defined $isomeric ? $isomeric: 0;

$multi = defined $multi?$multi:0;

my $decode = 0;

$donotconvert += 0;

$largest = defined $largest?$largest:0;

my $fix = 0;
my $usemass=0;
#my $fix = 20;
my $ID ="id";
my $INITIAL="initialreaction";

print "version 2.03 23/07/2020\n";

#$fileIn =~ /\.can/ && die("$fileIn should not contain .can");
#$fileRes =~ /\.can/ && die("$fileRes should not contain .can");

#(-e $fileIn) or die("File is absent $fileIn\n");

if(!defined $top){
    if(length($fileRes)>0 && int($fileRes) >0){
        $top = $fileRes;
        $fileRes = undef;
    }
}

if(!defined $fileIn){
    open(I,"config.cfg") or die ("use as: $0 test_set_file.csv <results.csv> <ntops> <isomeric = 0/all, 1/iso, 2/non > <largest> <multi> <augmentations> <beams> <skipcanonical: all, 1-yes 2-canonical> <donotconvert: 1> <single: 1>\n\n");
    while(<I>){
        chomp;
        /^\#/ and next;
        if(/apply_data_file\s+=\s+(\S+)/){$fileIn =$1;}
        if(/result_file\s+=\s+(\S+)/){$fileRes =$1;}
    }
    close(I);
}

if(!(defined $fileRes)){
    my @a = split(/\//,$fileIn);
    $fileRes = 	"result_$a[-1]";
} 

if (!(-e $fileRes)) {
    if(-e "result_patents_test20.csv.can"){
        $fileRes =~ /result_patents_test(\d+).csv/ or die("File $fileRes is absent");
        print "attempting conversion to N=$1\n";
        system("perl /transformer/marios/split.pl  result_patents_test20.csv $1 >/dev/null");
        system("perl /transformer/marios/split.pl  result_patents_test20.csv.can $1 >/dev/null");
    }
}

my (%AConvUnique,@toUnique,$number,$read,$all,$corr,%errors) =({},0,0,0,0,0);

my $initialCGRCanonicalFile = "";

my $globalerror = 0;


#$number = 1000;

$top = $top ?$top:1;

$augmentations += 0;
$beams += 0;

print "\nusing as:\nperl $0 $fileIn $fileRes $top $isomeric $largest $multi $augmentations   $beams $skipcanonical $donotconvert $single\n";
print "perl $0 $fileIn $fileRes best_top=$top isomeric=",$isomeric==0?"all":$isomeric==1?"iso":"no"," largest=",$largest?"yes":"no"," multi=",$multi?"yes":"no"," augmentations=",
$augmentations?$augmentations:"all"," beams=",$beams?$beams:"all"," canonical=", $skipcanonical==0?"all":$skipcanonical==1?"skip":"only",($donotconvert?" donotconvert":""),($single?" single":""),"\n\n";

my (@errors,$errors,$lines,$non,$empty,@ids); # will be used to skip errors in files
my ($allprocessed,$inchies);


my ($smilesI,$startsI,$hashesI) = compress(readSmiles($fileIn,1));
my ($smilesR,$startsT,$hashesR) = compress(readSmiles($fileRes,0));

defined $smilesI or die("$fileIn  does not contain data or absent");
defined $smilesR or die("$fileRes does not contain data or absent");

print "\n$fileIn\tmols: ",scalar @{$smilesI}," total lines:",`wc -l < $fileIn`;
print "$fileRes\tmols: ",scalar @{$smilesR}," total lines:",`wc -l < $fileRes`,"\n";

$decode && warn("N.B.! Files with \"I\" inside of names were back converted to original SMILES\n.");

system("rm inp >/dev/null 2>/dev/null");

scalar @errors && scalar @{$smilesI} != scalar @{$smilesR} && die("Input data contains error and sizes of input and target files are different\n");

open(R,">all.csv");
print R "mol,id,top,target,predicted,prob,result,inputdata\n";

my (@emptyi,@emptyr);

for(my $i = 0; $i<scalar @{$smilesI}; $i++){
    my $res = compareS($smilesI->[$i],$smilesR->[$i],$i);
    if(empty($smilesI->[$i])){
		push @emptyi,$i;
		#warn("empy initial: $i\n");
        $non++;
        next;
    }
    if(empty($smilesR->[$i])){
		#warn("empy results: $i\n");
		push @emptyr,$i;
        $empty++;
        next;
    }
    $all++;
    $corr += $res;
}

@emptyi && warn("empty initial ",scalar @emptyi,": @emptyi\n");
@emptyr && warn("empty results ",scalar @emptyr,": @emptyr\n");

my $dir = cwd;
warn("$dir   ", $isomeric==0?"all":$isomeric==1?"isomeric":"non-isomeric","\n");

warn $errors?"skipped: $errors processed: $lines ":"","ERRORS=$non NO-RES=$empty TP=$corr ALL=$all accuracy=",
	int(1000*$corr/$all+0.49999)/10,"(",int(1000*$corr/($all + $empty)+0.49999)/10,")% for TOP=$top predictions for ",scalar @{$smilesI}," processed molecules\n\n";

$fix && warn("CAUTION: $fix SMILES were forcedly merged for each molecule!\n");
$usemass && warn("CAUTION: CGRS with different atoms were not considered for analysis!\n\n");
#$allprocessed && warn("ALL=$allprocessed  ‎InChI=$inchies ",int($inchies*100/$allprocessed+0.499999),"%\n");

sub compareOne{
    my ($a,$e,$n,$top) = @_;
    
    my $res = 1;
    
    if(!defined $a || !defined $e || scalar @{$e} != scalar @{$a} || scalar @{$e} == 0){
        $res = 0; # not found!
    }
    
    for(my $j=0;$j<scalar @{$a};$j++){
        $a->[$j] == -1 && do{$res =0; last;}; # used for errors!
        $a->[$j] != $e->[$j] && do{$res =0; last;};
    }
   
   	print R "$n,$startsI->[$n],$top,";
    printMol($a);
    print R ",";
    printMol($e);
    print R ",$hashesR->[$n]->[$top],$res,$ids[$startsT->[$n]]->[2],$errors{$startsI->[$n]}->[0],$errors{$startsT->[$n]}->[1]\n";
    return $res;
}

sub compareMass{
    my ($a,$e) = @_;
    $usemass or return 1;
    (scalar @{$e} != scalar @{$a} || scalar @{$e} == 0) && return 0;
    
    my ($aa,$bb) = ("","");
    for(my $j=0;$j<scalar @{$a};$j++){
        $aa .= $toUnique[$a->[$j]];
        $bb .= $toUnique[$e->[$j]];
    };
    $aa =~ s/[^a-zA-Z0-9,]//g;
    $bb =~ s/[^a-zA-Z0-9,]//g;
    my @a= sort split(//,$aa);
    my @b= sort split(//,$bb);
    return "@a" eq "@b";
}

sub printMol{
    my ($a) = @_;
    for(my $j=0;$j<scalar @{$a};$j++){
        $j > 0 && print R ".";
        print R $a->[$j] == -1?"error": $toUnique[$a->[$j]];
    }
}

sub empty{
    my $smiles = shift;
    defined $smiles or return 1;
    scalar @{$smiles} > 0 or return 1;
    return $smiles->[0]->[0] == -1;
}

sub compareS{
    my ($a,$b,$nn) = @_;
    #print "top ",scalar @{$a},"\n";
    defined $b or return 0;
    for(my $l =0; $l< scalar @{$a}; $l++){ # comparing with any of possible variants of the training set
        my @a = @{$a->[$l]}; # always the first one
        my $k = 0;
        my $found = 0;
        for(my $i =0; $i< scalar @{$b} && $k< $top; $i++){
            #print "TOP $top $i -> ", scalar @{$b},"\n";
            my $e = $b->[$i];
            compareMass(\@a,$e) or next;
            if(compareOne(\@a,$e,$nn,$k)){
                $found = 1; # any match is good!
            };
           $k++;
        }
       $found && return 1;
    }
    return 0;
}

sub decode{
    $decode = 1;
    #print("decode\n");
    my $s = shift;
    $s =~ s/A/\[nH]/g;
    $s =~ s/T/Cl/g;
    $s =~ s/b/Br/g;
    $s =~ s/D/@@/g;
    $s =~ s/U/@/g;
    return $s;
}

sub isIsomeric{
    my $f = shift;
    
    $isomeric == 0 and return 1;
    
    if($isomeric == 1){
        return  ($f =~ /\\/ or $f =~ /\// or $f =~ /\@/)?1:0;
    }
    
    return  ($f =~ /\\/ or $f =~ /\// or $f =~ /\@/)?0:1;
}


sub readSmiles{
    my ($file, $experimental) = @_;
    my ($res,$smiles,$remove);
    
    my $can;
    
    if($experimental){
        open(I,$file) or die();
        my ($name,$initial,$input);
        while(<I>){
            chomp;
            my @a = split(/,/,$_);
            if(!defined $name){
                for(my $i=0;$i<scalar @a;$i++){
                    $a[$i] eq $ID && do{$name=$i};
                    $a[$i] eq $INITIAL && do{$initial=$i};
                    $a[$i] eq "input" && do{$input=$i};
                };
                next;
            }
            $ids[$.-2]->[0] = $a[$name];
            $ids[$.-2]->[1] = $a[$initial];
            $ids[$.-2]->[2] = $a[$input];
        }
        close(I);
    }
    
    if($file =~ /\.can/  || $donotconvert){
        $can = $file;
    }else{
        my $inititalFile= $file;
        $can = "$file.can";
        
        if( !open(I,$can)){
            print "converting $file to canonical $can\n";
            
            if($file =~ /I/){
                open(I,$file);
                open(W,">$file".".i");
                while(<I>){
                    length($_)>1 or next;
                    print W decode($_);
                }
                close(I);
                close(W);
                $remove = $file = "$file".".i";
            }elsif($file =~ /cgr/){
                print "converting CGR lines: ".`wc -l < $fileIn`;
                $file = convertCGR($file);
            }
            system("python canonical.py $initialCGRCanonicalFile $file > $can 2>/dev/null");
            my $size = -s $can;
            if($size < 10){
                system("rm $can");
                die("$inititalFile is empty\n");
            }
        }
    }
    
    open(I,$can) or die();
    
    $initialCGRCanonicalFile = $file =~ /cgr/ ? $can : "";
    
    my (@targert);
    my ($seq, $count,$err,$nn) = ("",0,0,0);
    while(<I>){
        chomp;
        $number && $. > $number && last;
        if(/input/){
            @targert = split(/,/,$_);
            next;
        }
        
        if(/error/){
            if($experimental){
                $errors[$nn]=1;$err++;
            }
            $errors{$nn}->[$experimental] = $_;
        }
        
        ($read && $read % 100 == 0)  && warn("processed $nn records from $file\n");
        chomp;
        s/\r//;
        my $line = $_;
        my @a =split(/,/,$line);
        for(my $i=0;$i<scalar @a;$i++){
            $a[$i] = standartiseReaction($a[$i]);
        }
        
        isIsomeric($a[0]) or do{
            #print "$a[0]\n";
            next;
        };
        
        if($seq eq $a[0]){
            $count++;
        }else{
            $seq = $a[0];
            $count = 1;
        }
        
        if(!$experimental && $skipcanonical > 0){
            $skipcanonical == 1 && $count == 1 && next; # skip canonical
            $skipcanonical != 1  && $count != 1 && next; # keep only canonical
        }
        
        $augmentations && $count > $augmentations && next; # skip those with more than > augmentations
        
        @targert or die("\n\ninput,target is absent in the file. Did you use correct input file: $can ?\n");
        
        
        my %a;
        for(my $n=0;$n<scalar @a;$n++){
            $beams && $n > $beams && next; # skip more than number of beams
            
            if($targert[$n] =~ /target/){ # store results
                my $val = convert($a[$n],$largest);
                if(!$multi){
                    $a{hash($val)} && next; # only one per sequence if not multi!
                }
                push @{$res->[$nn]},$val;
                $a{hash($val)}=1;
                #print "$nn  @{$res->[$nn]}\n";
            }elsif($targert[$n] =~ /input/i){
                push @{$smiles->[$nn]},convert($a[$n],0); # store also smiles
            }
        }
        
        
        $nn++;
    }
    close(I);
    if(defined $remove){
        system("rm $remove");
    }
    
    print "$can - processed: $nn errors: $err\n";
    
    return ($smiles,$res);
}

sub convertCGR{
    my $f = shift;
    
    my $cgr = "cgr.csv";
    
    open(I,$f) or die("Cannot open $f");
    
    my (%cgr, @names,@data,@ids,$id,%converted);
    while(<I>){
        chomp;
        $. == 1 && do{
            @names = split(/,/,$_);
            for($id=0;$id<scalar @names;$id++){
                $names[$id] eq $ID and last;	
            }
            next;
        };
        
        my @a = split(/,/,$_);
        
        for(my $i=1; $i< scalar @a;$i++){
            $names[$i] =~ /target/ or next;
            length($a[$i])>2 or next;
            $cgr{$a[$i]} and next;
            $cgr{$a[$i]} = 1; # will be substituted if corrected
            push @data,$a[$i];
            push @ids,$a[$id];
        }
    }
    close(I);
    
    my $script =  "canonizeCGR.py";
    
    open(C,">$cgr") or die("Cannot create $cgr");
    print C "cgr,$ID\n";
    
    for(my $i=0;$i<scalar @data;$i++){
        print C trim($data[$i]),",",(defined($ids[$i])?$ids[$i]:$i),"\n"; # for conversion
    }
    close(C);
    
    my $rescgr = "$cgr.converted.csv";
    
    system("python $script $cgr 2>/dev/null > $rescgr");
    open(R, $rescgr);
    
    my $n =0;
    while(<R>){
        /converted/ && next;
        chomp;
        if( !( /\{/ and /\}/ ) or /error:/ ){
            $cgr{$data[$n]} = "error:$data[$n]";
        }else{
            my @a=split(/,/,$_);
            $cgr{$data[$n]} = $a[0];
            $converted{$data[$n]}=$a[3];
        };
        $n++;
    }
    close(R);
    
    $n == scalar @data or die("conversion using $script failed $n !=".(scalar @data));
    
    #unlink($cgr);
    #unlink($rescgr);
    
    open(I,$f) or die;
    $f = "$f.cgr"; # contains cgr but not yet canonical input data
    open(O,">$f") or die;
    
    while(<I>){
        chomp;
        $. == 1 && do { print O; next;};
        my @a = split(/,/,$_);
        
        my $s = "\n".standartiseReaction($a[0]);
        my $first = 1;
        for(my $i=1; $i< scalar @a;$i++){
            $names[$i] =~ /target/ or next; # only those with target are used
            my $v = $a[$i];
            my $a = $cgr{$v};
            $s .= ",$a";
        }
        
        print O $s;
    }
    close(O);
    close(I);
    return ($f);
}


sub standartiseReaction{
    my $aa = shift;
    
    $aa =~ "error" && return "error";
    $aa =~ />>/ or return $aa;
    $aa = trim($aa);
    my @bb = split(/>>/,$aa);
    for(my $i=0;$i<scalar @bb; $i++){
        my @aa = split(/\./,$bb[$i]);
        @aa = sort { length $a <=> length $b || $a cmp $b } @aa;
        $bb[$i]="@aa";
        $bb[$i] =~ s/ /./g;
    }
    return length($bb[1])?"$bb[0]>>$bb[1]":$bb[0];
}

sub  trim{
    my $s = shift;
    $s =~ s/^\s+|\s+$//g;
    return $s
};

sub hash{
    my $n = shift;
    my $s ="hash";
    foreach my $e (@{$n}){
        $s .= $e."_";
    }
    return $s;
}

sub compress{
    my ($smiles,$res)=@_;
    
    my ($newres,$newsmiles,$starts,$hashes);
    my $nn=0;
    for(my $i = 0; $i<scalar @{$smiles};){
        my (%rr,$j,%val);
        
        while(!defined $smiles->[$i]->[0] || scalar @{$smiles->[$i]->[0]} == 0){ # failed molecule in the training set
            $errors[$i] != 1 && warn("Inconsistent error for smiles $i\n");
            $errors++;
            $i++;
        };
        
        for($j = $i; $j<$i+$fix; $j++){
            $smiles->[$j]->[0] = $smiles->[$i]->[0];
        }
        
        for($j = $i; $i == $j || ($j<scalar @{$smiles} && compareOne($smiles->[$i]->[0],$smiles->[$j]->[0],0,0));$j++){ # identical SMILES, predictions will be accumulated
            
            my $num = defined $res->[$j]? scalar @{$res->[$j]} : 0;
            
            for(my ($n,$c)=(0,1); $n<$num;$n++){
                my $hash = hash($res->[$j]->[$n]);
                $hash eq "hash" && next;
                if(!defined   $rr{$hash}){
                    $rr{$hash} = 0;
                    $val{$hash} = $res->[$j]->[$n];
                }
                $rr{$hash} += 1./(1. + 0.00001*$c++); # to provide always the same accuracy number; the SMILES at the top are getting the higher score
                #$rr{$hash} += 1./( ++$c + $j/9999999.); # to provide always the same accuracy number; the SMILES at the top are getting the higher score
                #$rr{$hash}++; # just counts
            }
            if($single){ $j++; last;};
            $lines++;
        }
        
        $res && $fix && (($j -$i) != $fix) && warn("lines $i -- $j were merged for ", ($j-$i)/$fix, " identical molecules\n");
        
        if(%rr){
            #print "found: ", scalar keys %rr,"\n";
            my @a = sort{$rr{$b} <=> $rr{$a}} keys %rr;
            foreach my $a (@a){
                my $aa = $val{$a};
                #print "$i $a @{$val{$a}} -- $rr{$a} @{$aa}\n";
                push @{$newres->[$nn]},$aa;
                push @{$hashes->[$nn]},int($rr{$a}); # the count itself
            }
        }else{
            my @a;$a[0]=-1; # for errors only
            push @{$newres->[$nn]},\@a;
            push @{$hashes->[$nn]},0; # the error
        }
        
        $starts->[$nn]=$i;
        $nn++;
        $i=$j;
    }
    
    return ($newres,$starts,$hashes);
}

sub convert{
    my ($s,$large) = @_;
    my @res;
    
    if($globalerror && $s =~ /error/){
    	return \@res;
    }
        
    my @a =split(/\./,$s);
    
    if($large){
        @a = sort { length $b <=> length $a || $a cmp $b } @a;
        length($a[1]) > length($a[0]) && die("error in sorting @a");
        #print "@a --> $a[0]\n";
        $#a=0; # only the first is left
    }
    
    foreach my $e (@a){
        
        length($e) < 1 && next;
        $e =~ /error/ && next; # we are not interested in the failed by any reason predictions!
        
        my $full = $e; # FULL SMILES or CGR
        
        if( /\{/ and /\}/){ # CGRsssss
			$allprocessed++;

    	    if($e =~ /\|/){ # special case to use InChies together with CGRs
        		$inchies++;
        		my @a = split(/\|/,$e);
        		$e = $a[1];
       	 	}
       	 };
        
        if(!defined $AConvUnique{$e}){
            $AConvUnique{$e} = scalar keys %AConvUnique;
            $toUnique[$AConvUnique{$e}] = $full;
        }
        push @res,$AConvUnique{$e};
    }
    @res = sort @res;
    return \@res;
}
