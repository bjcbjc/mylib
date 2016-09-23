
bcftools view -h ${1} | grep -i gatk | gawk '{match($0, /Version=([0-9A-za-z\-\.]+)/, a); print a[1]}'
