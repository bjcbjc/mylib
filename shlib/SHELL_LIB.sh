
### LOG_CMD name command
function LOG_CMD {
    echo "===== ${1} Start. =====" $(date)
    tic=`date +%s`
    echo "LOGGED_CMD=${2}"
    eval "${2}"
    status=$?
    toc=`date +%s`
    elapsed=`expr $(( ${toc} - ${tic} ))`
    if [ ${status} == 0 ]; then
	echo "===== ${1} Done in ${elapsed} sec. =====" $(date)
    else
	echo "===== ${1} Failed in ${elapsed} sec. Code: ${status} =====" $(date)
	exit 1
    fi
}

##LOG_CMD "sleep 10 test1" "sleep 10"
##LOG_CMD "sleep not working" "sleep ?"


function CHECK_PATH {
    if [ ! -d "${1}" ]; then
	mkdir -p ${1}
    fi
}
