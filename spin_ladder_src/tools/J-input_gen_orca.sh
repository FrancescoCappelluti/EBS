#!/bin/sh

DATAPATH="/home/narzi/SinglePoint_ModelB/"

rm M_tab.dat
rm Energies.dat

#Make a Loop over all the directory with distinct BS state calculations
#for i in Up-Up-Up-Up Up-Down-Up-Down Up-Down-Down-Up Up-Up-Down-Down Down-Up-Up-Up Up-Down-Up-Up Up-Up-Down-Up Up-Up-Up-Down
for i in 0 1 2 3 4 5 6 7
do
    # Extract from the output the MULLIKEN SPIN POPULATION for each Metal and the DFT Energy
#    echo $i >> M_tab.dat
    grep "  0 Fe: " $DATAPATH/$i/out | awk '{printf("%10.6lf ",$7) }' >> M_tab.dat
    grep "  1 Fe: " $DATAPATH/$i/out | awk '{printf("%10.6lf ",$7) }' >> M_tab.dat
    grep "  2 Fe: " $DATAPATH/$i/out | awk '{printf("%10.6lf ",$7) }' >> M_tab.dat
    grep "  3 Fe: " $DATAPATH/$i/out | awk '{printf("%10.6lf \n",$7) }' >> M_tab.dat

    grep "FINAL SINGLE POINT ENERGY " $DATAPATH/$i/out | awk '{printf("%20.12lf \n",$9) }' >> Energies.dat

done

    awk 'function abs(value){return (value<0?-value:value) }
    {a1+=abs($1);s1+=$1*$1;a2+=abs($2);s2+=$2*$2;a3+=abs($3);s3+=$3*$3;a4+=abs($4);s4+=$4*$4}
    END{a1/=8;s1/=8;a2/=8;s2/=8;a3/=8;s3/=8;a4/=8;s4/=8;
    s1-=(a1*a1);s2-=(a2*a2);s3-=(a3*a3);s4-=(a4*a4);
    printf(" %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf \n",a1,sqrt(s1),a2,sqrt(s2),a3,sqrt(s3),a4,sqrt(s4))}' M_tab.dat > M_ave.dat

    cat M_ave.dat M_tab.dat | awk  '{if(NR==1){m1=$1;m2=$3;m3=$5;m4=$7}else{
    printf("%10.6lf %10.6lf %10.6lf %10.6lf\n",($1<0?-m1:m1),($2<0?-m2:m2),($3<0?-m3:m3),($4<0?-m4:m4))}}' > M_tab-new.dat

    awk '{printf(" %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf \n",$1*$2/4,$1*$3/4,$1*$4/4,$2*$3/4,$2*$4/4,$3*$4/4)}' M_tab-new.dat > M_matrix.dat
