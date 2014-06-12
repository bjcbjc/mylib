BEGIN { \
    split(col,columns, ":"); 
    split(formatname, formatfields, ":");
    #split(filterexactmatch, filtermatches, ":");
    #pat = "$1 !~ /^#/ "
    #filterpat = ""
    #for (i=1; i<=length(filtermatches); i++) {
#	if (i==1) {
#	    filterpat = (filterpat "$7 == \"" filtermatches[i] "\" ")
#	} else {
#	    filterpat = (filterpat "|| $7 == \"" filtermatches[i] "\" ")
#	}
#    }
#    if (filterpat != "") {
#        pat = ( pat "&& (" filterpat ")" )
#    }
}
$1 !~ /^#/ {					\
   oidx = 1
   delete output
   split($9, format, ":")
   for (colidx=10; colidx<=NF; colidx++) {
       split($colidx, data, ":")
       for (i=1; i<=length(format); i++) {
       	   for (j=1; j<= length(formatfields); j++) {
       	       if (format[i] == formatfields[j]) {
	       	  output[oidx++] = data[i]
	       }  
	   }
       }
   }
   for (i=1; i<=length(columns); i++) {
       printf("%s\t",$(columns[i]))
   }
   for (i=1; i<=length(output); i++) {
       if (output[i] == "") {
	   printf(".\t")
       } else {
	   printf("%s\t",output[i])
       }
   }
   printf("\n")
}

