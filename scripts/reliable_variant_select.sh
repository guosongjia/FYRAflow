##############################################################
#  script: reliable_variant_select.sh
#  author: Guo-Song Jia
#  last edited: 2021.07.12
#  usage: bash reliable_variant_select.sh [software_name(samtools, gatk, deepvariant)] [vcfFile] [saving_PATH] [appendix] 
#  description: Script for in FYRAflow workflow. Extract reliable variant sites from different callers and then normalize them.
##############################################################

samtools_gatk_vcf_filtering () {
    whatshap unphase $1 > unphased
    cat unphased |grep "#" > header
    vcffilter -f "DP > 10" $1 |grep "1/1"  > lines
    cat header lines > vcfFile
    vcfallelicprimitives vcfFile > $2/$3.$4.final.vcf
    rm unphased header lines vcfFile
}

deepvariant_vcf_filtering () {
    whatshap unphase $1 > unphased
    cat unphased |grep "#" > header
    cat unphased |grep "PASS"|grep "1/1" > lines
    cat lines |while read line;do
        dp=$(echo $line|awk '{print $10}' |cut -d ":" -f 3)
        if [ ${dp} -ge 10 ];then
            echo $line >> filter
        fi
    done
    sudo cat header filter > $2/$3.$4.final.vcf
    rm unphased header lines filter
}


if [ $1 == "samtools" ];then
    samtools_gatk_vcf_filtering $2 $3 $4 samtools
elif [ $1 == "gatk" ];then
    samtools_gatk_vcf_filtering $2 $3 $4 gatk
elif [ $1 == "deepvariant" ];then
    deepvariant_vcf_filtering $2 $3 $4 deepvariant
else
    echo "Lack parameters!"
fi