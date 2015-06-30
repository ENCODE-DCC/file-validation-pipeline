#!/bin/bash
# comp-md5sum 1.0.0

#set -x
#set +e

main() {
    echo "*****"
    echo "* Running: comp-md5sum.sh [1.0.0]"
    #echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    #echo "* also using: diff, ls, md5sum, tee, diss and sammy"
    echo "*****"

    fileA_fn=`dx describe "$fileA" --name`
    fileB_fn=`dx describe "$fileB" --name`
    fileA_root=${fileA_fn%.*}
    fileB_root=${fileB_fn%.*}
    fileA_ext=${fileA_fn#*.}
    #fileB_ext=${fileB_fn#*.}
    log_diff_fn=compMd5_${fileA_root}_${fileB_root}.txt

    if [ "${fileA_fn}" == "${fileB_fn}" ]; then
        # avoid problems with identically named files
        fileA_root="fileA"
        fileB_root="fileB"
    fi
    echo "* Download ${fileA_fn} to '${fileA_root}.${fileA_ext}'..."
    dx download "$fileA" -o ${fileA_root}.${fileA_ext}
    echo "* Download ${fileB_fn} to '${fileB_root}.${fileA_ext}'..."
    dx download "$fileB" -o ${fileB_root}.${fileA_ext}
    
     echo "* sort_first: [$sort_first]"  | tee -a ${log_diff_fn}
    if [ "$sort_first" == "true" ]; then
        echo "* Simple sorting files first..." | tee -a ${log_diff_fn}
        set -e
        sort < ${fileA_root}.${fileA_ext} > ${fileA_root}.sorted
        mv ${fileA_root}.sorted ${fileA_root}.${fileA_ext}
        sort < ${fileB_root}.${fileA_ext} > ${fileB_root}.sorted
        mv ${fileB_root}.sorted ${fileB_root}.${fileA_ext}
        set +e
    fi

    echo " " >>  ${log_diff_fn}
    echo "* Comparing by md5sum [$fileA_fn] to [$fileB_fn]..." | tee -a ${log_diff_fn}
    #echo "- Lines:" | tee -a ${log_diff_fn} 
    #wc -l *.${fileA_ext} | tee -a ${log_diff_fn} 
    echo "- md5sum:" | tee -a ${log_diff_fn} 
    md5sum *.${fileA_ext} | tee -a ${log_diff_fn}

    if [ "$sort_first" == "true" ]; then
        echo "- Lines:" | tee -a ${log_diff_fn} 
        wc -l *.${fileA_ext} | tee -a ${log_diff_fn} 
        echo "- Lines that differ:" | tee -a ${log_diff_fn} 
        diff --speed-large-files --suppress-common-lines -y ${fileA_root}.${fileA_ext} ${fileB_root}.${fileA_ext} | wc -l | tee -a ${log_diff_fn}
    fi

    echo "* Uploading results..."
    log_diff=$(dx upload ${log_diff_fn} --brief)
    dx-jobutil-add-output log_diff "$log_diff" --class=file
    echo "* Finished."
}
