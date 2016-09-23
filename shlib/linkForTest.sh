

## linkForTest.sh source_project target_project

source=${1}
target=${2}

mkdir -p ${target}


for samplePath in `ls -d ${source}/Sample_*/`; do
    sample=`echo ${samplePath} | cut -d"/" -f6`
    mkdir -p ${target}/${sample}/analysis
    for f in ${samplePath}/analysis/*.haplotypeCalls.er.raw.vcf*; do
	ln -s ${f} ${target}/${sample}/analysis/$(basename $f)
    done

    mkdir -p ${target}/${sample}/reports
    for f in ${samplePath}/reports/*; do
	ln -s ${f} ${target}/${sample}/reports/$(basename $f)
    done

done

if [[ -e ${source}/sample.jg.list ]]; then
    cp ${source}/sample.jg.list ${target}
fi
