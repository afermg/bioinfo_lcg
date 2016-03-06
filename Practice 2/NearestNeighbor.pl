#Completar el código con el diagrama y hacer una predicción con los sitios de e.coli (formato consensus, secuencia+locus tag, -400 +50, no sabemos el inicio de
#la transcripción)

#!/usr/bin/perl -w

# prog1.1 
# Bruno Contreras-Moreira
# Nearest Neighbor dG calculator

use strict;

# global variables
my $T           = 37; # temperature(C)
my $windowL     = 15;  # window length, http://www.biomedcentral.com/1471-2105/6/1
my ($i,$j,$wind,$temp_deltas,$temp_sum,$cutoff1,$cutoff2,$prom_num,$cont);
my (@deltas,@E1,@E2,@D,@prom_pos,@prom_seq, @prom_all,@prom_all_pos);
my %NNparams    = ( 
	# SantaLucia J (1998) PNAS 95(4): 1460-1465.
	# [NaCl] 1M, 37C & pH=7 
	# H(enthalpy): kcal/mol	, S(entropy): cal/kmol
	# stacking dinucleotides
	'AA/TT' , {'H',-7.9, 'S',-22.2},
	'AT/TA' , {'H',-7.2, 'S',-20.4},
	'TA/AT' , {'H',-7.2, 'S',-21.3},
	'CA/GT' , {'H',-8.5, 'S',-22.7},
	'GT/CA' , {'H',-8.4, 'S',-22.4},
	'CT/GA' , {'H',-7.8, 'S',-21.0},
	'GA/CT' , {'H',-8.2, 'S',-22.2},
	'CG/GC' , {'H',-10.6,'S',-27.2},
	'GC/CG' , {'H',-9.8, 'S',-24.4},
	'GG/CC' , {'H',-8.0, 'S',-19.9},
	# initiation costs
	'G'     , {'H', 0.1, 'S',-2.8 },
	'A'     , {'H', 2.3, 'S',4.1  },
	# symmetry correction
	'sym'   , {'H',   0, 'S',-1.4 } );

my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";

print "# parameters: Temperature=$T\C Window=$windowL\n\n";

open(SEQ, $infile) || die "# cannot open input $infile : $!\n";

my $cutoff1= 3.4;
my $cutoff2= -15.99;

$cont = 0;
open(my $FILE, '>',"Data/E1.txt");
open(my $FILE2, '>',"Data/E2.txt");
open(my $FILE3, '>',"Data/D.txt");



while(<SEQ>)
{
	$cont++;
	if(/^(b\d{4}) \\ ([ATGC]+)/)
	{
		my ($name,$seq) = ($1,$2); 
		#Borramos los datos en nuestros vectores para hacer iteraciones.
		undef(@deltas);
		undef(@E1);
		undef(@E2);
		undef(@D);
		undef(@prom_pos);
		undef(@prom_seq);
		#hacer la ventana -> ir moviendola de n en n.
		#la n, siempre se refiere a posicion
		
		#Creamos las ventanas a lo largo de la secuencia, a excepción de las orillas del tamaño de mitad de la ventana.
		for ($i=int($windowL/2); $i < length($seq)-int($windowL/2); $i++){
				$wind = substr($seq,$i-$windowL/2,$windowL);
				#llamas a la subrutina duplex_deltaG e introducimos el valor a nuestro vector de deltas G
				push @deltas,duplex_deltaG($wind,$T);	
		}
		#Comenzamos en nucleótido 51. 
		#Hacemos el cálculo para cada n a partir de ahí hasta 150 antes de acabar 
		#(para poder obtener valores E2 y D también).
		#Luego introducimos los valores E1 y E2 a los vectores @E1 y @E2
		
		#E1=SUM (n a n+49)*$total_dG/50
		for ($i=50; $i < length($seq)-150;$i++){
			$temp_sum = 0;
			for ($j=$i-50;$j<$i;$j++){
			$temp_sum += $deltas[$j];
			}
			$temp_deltas = $temp_sum/50;
			push(@E1,$temp_deltas);
		}
		
		#E2=SUM (n+99 a n+199)*$total_dG/100
		for ($i=50; $i < length($seq)-150;$i++){
			$temp_sum = 0;
			for ($j=$i+50;$j<$i+150;$j++){
			$temp_sum += $deltas[$j];
			}
			$temp_deltas = $temp_sum/100;
			push(@E2,$temp_deltas);
		}
		
		#Obtenemos los valores D y los introducimos a su propio vector.
		#D(n)=E1(n)-E2(n)
		$prom_num = 0;
		for ($i=0;$i < scalar(@E1);$i++){
			push(@D,$E1[$i]-$E2[$i]);
			
			#Despúes, usando E1 y D rescatamos los motivos de acuerdo a los valores de cutoff
			if ($D[$i]>$cutoff1 && $E1[$i]>$cutoff2) { 			
				#Obtenemos secuencia y valor real (sumando los windowL/2 que se omitieron al inicio)
				push @prom_seq,substr($seq,$i+50,$windowL);
				push @prom_pos,($i+50+int($windowL/2));

				$prom_num++;
				
				#if (n>) Si n está a 25n a -25n
				if ($prom_num>1){
					#Si están demasiado cerca los consideramos el mismo motivo.
					if ($prom_pos[$prom_num-1]-$prom_pos[$prom_num-2]<25){
						pop @prom_pos;
						pop @prom_seq;
						$prom_num--;
					}
				}

			}
		  
			
		}
		
	


	}
	push @prom_all,@prom_seq;
	push @prom_all_pos,@prom_pos;
	if ($cont==1){
		for my $i (1..$#E1+1){
			printf $FILE "%d\t", $i;
			printf $FILE2 "%d\t", $i;
			printf $FILE3 "%d\t", $i;
		}
	
		printf $FILE "\n";
		printf $FILE2 "\n";
		printf $FILE3 "\n";
	}

	for my $i (0..$#E1) {
	printf $FILE "%f\t", $E1[$i];
	printf $FILE2 "%f\t", $E2[$i];
	printf $FILE3 "%f\t", $D[$i];
	}
	
	printf $FILE "\n";
	printf $FILE2 "\n";
	printf $FILE3 "\n";
}
close(SEQ);
close($FILE);
close($FILE2);
close($FILE3);
#print "@prom_all \n";

open(my $FILE4, '>',"predicted_promoters.txt");
printf $FILE4 "Promoter\tPos\tSequence\n";
for my $i (0..$#prom_all) {
	printf $FILE4 "%d\t%d\t%s\n",$i+1, $prom_all_pos[$i], $prom_all[$i];
}	
close($FILE4);

# calculate NN free energy of a DNA duplex , dG(t) = (1000*dH - t*dS) / 1000
# parameters: 1) DNA sequence string; 2) Celsius temperature
# returns; 1) free energy scalar
# uses global hash %NNparams
sub duplex_deltaG
{
   	my ($seq,$tCelsius) = @_; 
	
	my ($DNAstep,$nt,$dG,$total_dG) = ('','',0,0);
	my @sequence = split(//,uc($seq));
	my $tK = 273.15 + $tCelsius;
	my $i;
	
	sub complement{ $_[0] =~ tr/ATGC/TACG/; return $_[0] }
	
	# add dG for overlapping dinculeotides
	for(my $n=0;$n<$#sequence;$n++) 
	{
			$DNAstep = $sequence[$n].$sequence[$n+1].'/'.
				complement($sequence[$n].$sequence[$n+1]);
			
			if(!defined($NNparams{$DNAstep}))
			{
				$DNAstep = reverse($DNAstep);
			}
			

			$dG = ((1000*$NNparams{$DNAstep}{'H'})-
					($tK*$NNparams{$DNAstep}{'S'}))
					/ 1000 ;
			
			$total_dG += $dG; 
	}
	
	# add correction for helix initiation
	$nt = $sequence[0]; # first pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) } 
	$total_dG += ((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'}))
					/ 1000; 
	
	$nt = $sequence[$#sequence]; # last pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) }
	$total_dG += ((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'}))
					/ 1000;
					
	# please complete for symmetry correction
	if($seq == reverse($seq)){
		$total_dG += ((1000*$NNparams{'sym'}{'H'})-
					($tK*$NNparams{'sym'}{'S'}))
					/ 1000;
	}
	
	
	
	return $total_dG;
}

#------------------------------------------------Diccionario de Variables--------------------------------------------------------------#
=pod
$T           = 37; Temperatura
$windowL     = 15;  Tamaño de la ventana
$cutoff1= 3.4; Valor de corte para E1
$cutoff2= -15.99; Valor de corte para D
$i	Contador 1, usado para recorrer las posiciones de las que se obtendrá dG y éstos mismos valores en su vector.
$j	Contador 2, usado para hacer sumatorias de dG en cálculos de E1 y E2
$wind Valor scalar temporal que contiene las secuencias de las ventanas para introducir a la función duplex_deltaG
$temp_deltas	Valor escalar temporal para guardar el valor final al obtener E1 y E2
$temp_sum	Valor escalar temporal para guardar la sumatoria de dGs y mandarlas a temp_deltas
$prom_num	Valor escalar, contador usado para evaluar posiciones mediante el vector prom_pos
$count Valor escalar, contador usado para insertar los headers a los archivos exportados hacia R
@deltas Vector unidimensional para guardar los valores dG
@E1	Vector unidimensional para guardar los valores E1
@E2 Vector unidimensional para guardar los valores E2
D	Vector unidimensional para guardar los valores D
@prom_pos Vector unidimensional para guardar las posiciones de los promotores seleccionados en cada secuencia.
@prom_seq Vector unidimensional para guardar las secuencias de los promotores seleccionados en cada secuencias.
@prom_all Vector unidimensional para guardar las secuencias de todos los promotores a lo largo del análisis.

=cut
