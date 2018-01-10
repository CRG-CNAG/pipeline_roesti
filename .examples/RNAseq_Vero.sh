export PATH=$PATH:/software/ls/MAQ/:/software/ls/

if [ "$#" -ne 2 ]
then
        echo "Performs the pre-processing and mapping of paired-end to obtain the pile-ups." 
        echo "ARGS: (both mandatory)"
        echo "1. Read_1 file"
        echo "2. Read_2 file"
        exit 1
else
        ##source ~/.myfunctions
        r1=`basename $1`
        r2=`basename $2`

        ###fastq to bfq
        outbfq=`echo $r1 | awk '{ sub(/fastq/, "bfq"); print}'`
        outbfq2=`echo $r2 | awk '{ sub(/fastq/, "bfq"); print}'`
        maq fastq2bfq $r1 $outbfq
        maq fastq2bfq $r2 $outbfq2
        echo "Fastq converted to bfq. Getting indexed genome sequence..."

        ###get genome.bfa (copy from elsewhere)
        cp ~/NC_000912.bfa .
        echo "Reference genome index is correct. Mapping with MAQ..."

        ###maq map
        outmap="$outbfq.map"
        outmap2="$outbfq2.map"
        maq map -1 40 -n 1 $outmap NC_000912.bfa $outbfq
        maq map -1 40 -n 1 $outmap2 NC_000912.bfa $outbfq2

        echo "Mapping finished successfully! Converting to text..."

        ###map to txt
        outtxt="$outmap.txt"
        outtxt2="$outmap2.txt"

        maq mapview $outmap > $outtxt
        maq mapview $outmap2 > $outtxt2

        outfiltered="$outtxt.filt"
        outfiltered2="$outtxt2.filt"

        cat $outtxt | awk '{if ($12==1 || ($12==0 && $13==1)) print $0}' > $outfiltered  ##takes only reads that uniquely map to 1 position, either with no mismatches or only with one
        cat $outtxt2 | awk '{if ($12==1 || ($12==0 && $13==1)) print $0}' > $outfiltered2  ##takes only reads that uniquely map to 1 position, either with no mismatches or only with one

        echo "Map converted to text, calculating pileup..."

        ###pileup

        result=`echo $r1 | awk '{ sub(/read1.fastq/, "read.pile"); print}'`
        pileupExcp $outfiltered $outfiltered2 40 816394 $result 0 MAQPESS
        echo "Pileup calculated. You are done!"
fi 
