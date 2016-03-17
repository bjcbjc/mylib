
BEGIN {
    if (length(col)==0) {
	print "need col index (1-based)"
    }
    if (length(key)==0) {
	print "need key (file with keyword)"
    }
    while ((getline line < key) > 0) {
	keytb[line] = 1
    }
}

$col in keytb { print }
