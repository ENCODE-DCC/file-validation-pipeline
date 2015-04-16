#!/bin/bash
# comp-bams 1.0.0

#set -x
#set +e

main() {
    #echo "Value of test dataset: '$test_dir'"
    #echo "Value of standard dataset: '$data_dir'"
    #dx download project-BKjyBz00ZZ0PZF5V7Gv00zqG:/"$test_dir"/*
    #dx download project-BKjyBz00ZZ0PZF5V7Gv00zqG:/test/"$data_dir"/*
    #find

    echo "*****"
    echo "* Running: comp-bams.sh [1.0.0]"
    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "* also using: diff, ls, md5sum, tee, diss and sammy"
    echo "*****"

    echo "* md5sum of all reads: '$all'"
    echo "* md5sum of all reads sorted: '$allsort'"
    echo "* md5sum of unique reads: '$uniqs'"
    echo "* md5sum of mutli-mapped reads: '$multis'"
    echo "* md5sum of unmapped reads: '$unmapped'"
    echo "* Diff of bam file headers: '$headers'"

    bamA_fn=`dx describe "$bamA" --name`
    bamB_fn=`dx describe "$bamB" --name`
    bamA_root=${bamA_fn%.bam}
    bamB_root=${bamB_fn%.bam}
    log_diff_fn=compBams_${bamA_root}_${bamB_root}.txt

    if [ "${bamA_fn}" == "${bamB_fn}" ]; then
        # avoid problems with identically named files
        bamA_root="bamA"
        bamB_root="bamB"
    fi
    echo "* Download ${bamA_fn} to '${bamA_root}.bam'..."
    dx download "$bamA" -o ${bamA_root}.bam
    echo "* Download ${bamB_fn} to '${bamB_root}.bam'..."
    dx download "$bamB" -o ${bamB_root}.bam

    if [ "$all" == "true" ] || [ "$allsort" == "true" ]; then 
        echo " " >>  ${log_diff_fn}
        rm -f *.sam
        if [ "$allsort" == "true" ]; then
            echo "* Sorting and comparing all reads [$bamA_fn] to [$bamB_fn]..." | tee -a ${log_diff_fn}
            sammy --all ${bamA_root}.bam --fast  # sammy sorts
            sammy --all ${bamB_root}.bam --fast
        else
            echo "* Comparing all reads [$bamA_fn] to [$bamB_fn]..." | tee -a ${log_diff_fn}
            samtools view -@ 8 ${bamA_root}.bam > ${bamA_root}.sam 
            samtools view -@ 8 ${bamB_root}.bam > ${bamB_root}.sam
        fi
        echo "- Lines:" | tee -a ${log_diff_fn} 
        wc -l *.sam | tee -a ${log_diff_fn} 
        echo "- md5sum:" | tee -a ${log_diff_fn} 
        md5sum *.sam | tee -a ${log_diff_fn}
    fi

    if [ "$uniqs" == "true" ]; then 
        echo " " >>  ${log_diff_fn}
        echo "* Comparing uniquely mapped [$bamA_fn] to [$bamB_fn]..." | tee -a ${log_diff_fn}
        sammy --uniq ${bamA_root}.bam --fast 
        sammy --uniq ${bamB_root}.bam --fast
        echo "- Lines:" | tee -a ${log_diff_fn} 
        wc -l *_uniq.sam | tee -a ${log_diff_fn} 
        echo "- md5sum:" | tee -a ${log_diff_fn} 
        md5sum *_uniq.sam | tee -a ${log_diff_fn}
        #echo "Split and diff:" | tee -a ${log_diff_fn} 
        #rm -f splitFile?_* 
        #split -l 10000000 ${bamA_root}_uniq.sam splitFile1_ 
        #split -l 10000000 ${bamB_root}_uniq.sam splitFile2_ 
        #for f in `ls splitFile1_??`; do
        #    diss $f splitFile2_${f#splitFile1_} | tee -a ${log_diff_fn}
        #done
        ##diff ${bamA_root}_uniq.sam ${bamB_root}_uniq.sam | head | wc -l | tee -a ${log_diff_fn} 
    fi

    if [ "$multis" == "true" ]; then 
        echo " " >>  ${log_diff_fn}
        echo "* Comparing multi-mapped [$bamA_fn] to [$bamB_fn]..." | tee -a ${log_diff_fn}
        sammy --multi ${bamA_root}.bam --fast 
        sammy --multi ${bamB_root}.bam --fast 
        echo "- Lines:" | tee -a ${log_diff_fn} 
        wc -l *_multi.sam | tee -a ${log_diff_fn} 
        echo "- md5sum:" | tee -a ${log_diff_fn} 
        md5sum *_multi.sam | tee -a ${log_diff_fn}
        echo "- Split and diff:" | tee -a ${log_diff_fn} 
        rm -f splitFile?_* 
        split -l 10000000 ${bamA_root}_multi.sam splitFile1_ 
        split -l 10000000 ${bamB_root}_multi.sam splitFile2_ 
        for f in `ls splitFile1_??`; do
            diss $f splitFile2_${f#splitFile1_} | tee -a ${log_diff_fn}
        done
        ##diff ${bamA_root}_multi.sam ${bamB_root}_multi.sam | head | wc -l | tee -a ${log_diff_fn} 
        #echo "- Without sam flag:" | tee -a ${log_diff_fn} 
        #cat ${bamA_root}_multi.sam | cut -f1,3- | grep -v -w NH\:i\:20 | sort > ${bamA_root}_multiClip.sam
        #cat ${bamB_root}_multi.sam | cut -f1,3- | grep -v -w NH\:i\:20 | sort > ${bamB_root}_multiClip.sam
        #diss ${bamA_root}_multiClip.sam ${bamB_root}_multiClip.sam | tee -a ${log_diff_fn} 
        #echo "- Without multi >= 20:" | tee -a ${log_diff_fn} 
        #grep -v -w NH\:i\:20 ${bamA_root}_multiClip.sam > ${bamA_root}_multiClipN20.sam
        #grep -v -w NH\:i\:20 ${bamB_root}_multiClip.sam > ${bamB_root}_multiClipN20.sam
        #diss ${bamA_root}_multiClipN20.sam ${bamB_root}_multiClipN20.sam | tee -a ${log_diff_fn} 
    fi

    if [ "$unmapped" == "true" ]; then 
        echo " " >>  ${log_diff_fn}
        echo "* Comparing unmapped [$bamA_fn] to [$bamB_fn]..." | tee -a ${log_diff_fn}
        sammy --unmapped ${bamA_root}.bam --fast 
        sammy --unmapped ${bamB_root}.bam --fast 
        echo "Lines:" | tee -a ${log_diff_fn} 
        wc -l *_unmapped.sam | tee -a ${log_diff_fn} 
        echo "md5sum:" | tee -a ${log_diff_fn} 
        md5sum *_unmapped.sam | tee -a ${log_diff_fn}
        #echo "Split and diff:" | tee -a ${log_diff_fn} 
        #rm -f splitFile?_* 
        #split -l 10000000 ${bamA_root}_unmapped.sam splitFile1_ 
        #split -l 10000000 ${bamB_root}_unmapped.sam splitFile2_ 
        #for f in `ls splitFile1_??`; do
        #    diss $f splitFile2_${f#splitFile1_} | tee -a ${log_diff_fn}
        #done
        ##diff ${bamA_root}_unmapped.sam ${bamB_root}_unmapped.sam | head | wc -l | tee -a ${log_diff_fn} 
    fi
    
    if [ "$headers" == "true" ]; then 
        echo " " >>  ${log_diff_fn}
        echo "* Comparing headers [$bamA_fn] to [$bamB.bam]..." | tee -a ${log_diff_fn}
        sammy --head ${bamA_root}.bam 
        sammy --head ${bamB_root}.bam
        grep -v \@PG ${bamA_root}_head.sam | grep -v "user command" > ${bamA_root}_clean_head.sam
        grep -v \@PG ${bamB_root}_head.sam | grep -v "user command" > ${bamB_root}_clean_head.sam
        diss ${bamA_root}_clean_head.sam ${bamB_root}_clean_head.sam | tee -a ${log_diff_fn}
        #diff ${bamA_root}_clean_head.sam ${bamB_root}_clean_head.sam | head | wc -l | tee -a ${log_diff_fn} 
        #diff ${bamA_root}_clean_head.sam ${bamB_root}_clean_head.sam | tee -a ${log_diff_fn}
    fi


    echo "* Uploading results..."
    log_diff=$(dx upload ${log_diff_fn} --brief)
    dx-jobutil-add-output log_diff "$log_diff" --class=file
    echo "* Finished."
}
