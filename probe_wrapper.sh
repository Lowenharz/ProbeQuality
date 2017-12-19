#!/bin/bash
# Wrapper for probe_quality_report.pl. Intake the the path to Bevmo Variants folder, manifest.txt, and probelist.csv
# create SampleManifest.txt, and call probe_quality_report.pl
#
# If have any question, contact: qli2@illumina.com
usage(){
    echo "To use probe_wrapper: "
    echo "probe_wrapper.sh -v /PATH/TO/bevmo_version_Analysis/Variants/ -m /PATH/TO/Manifest.txt -p /PATH/TO/probelist.csv"
}
if [[ $1 == "" ]];
then
    echo "Must have arguement!"
    usage
    exit
fi

while getopts ":v:m:p:h" options; do
    case "${options}" in
        h )
            usage
            exit 0
            ;;
        v) 
            VARIANT_FOLDER=$OPTARG
            ;;
        m) 
            MANIFEST=$OPTARG
            ;;
        p)
            PROBE_LIST=$OPTARG
            ;;
    esac
done

if [ ! -f $MANIFEST ]; 
then 
    echo $MANIFEST" does not exist. Please check."
    exit
fi

if [ ! -f $PROBE_LIST ]; 
then 
    echo $PROBE_LIST" does not exist. Please check."
    exit
fi

list=()
parent_dir=
if [ -d $VARIANT_FOLDER ]; 
then 
    D2=`dirname $VARIANT_FOLDER`
    parent_dir=`dirname $D2`
    # echo $parent_dir
    samplesheet=$parent_dir"/SampleSheet.csv"
    # echo $samplesheet
    if [ -f $samplesheet ];
    then 
        # Read csv for Samplesheet.csv 
        # cat $samplesheet 
        is_sample_id=false
        rm $parent_dir"/"SampleManifest.txt
        while read line
            do
                sample_name=`echo $line | cut -d, -f2`
                if [ "$sample_name" == "Sample_ID" ]; 
                then 
                    is_sample_id=true
                fi
                if [ "$is_sample_id" = true ]; 
                then 
                    list+=($sample_name)
                fi
        done < $samplesheet 
        # echo ${list[@]}
        for (( i=1; i<${#list[@]}; i++ ));
        do 
            echo $i","${list[$i]}","$MANIFEST","$PROBE_LIST",truth">> $parent_dir"/"SampleManifest.txt
        done
    else
        $samplesheet " does not exist. Please check."
    fi
else 
    echo $VARIANT_FOLDER" does not exist. Please check."
    exit    
fi

CMD="/illumina/scratch/K2-I/software/probe_report/probe_quality_report.pl $VARIANT_FOLDER $parent_dir/SampleManifest.txt"

eval $CMD
