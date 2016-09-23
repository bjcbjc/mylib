
args = commandArgs(trailingOnly= TRUE)

fn = args[1]
if (length(args) > 1) {
    outfn = args[2]
} else {
    outfn = paste0(sub('.txt$','',fn), '.pdf')
}
outtxtfn = sub('.txt.pdf$', '.pdf', outfn)
outtxtfn = sub('.pdf$', '.summary.txt', outtxtfn)


get_analysis_name <- function(jobnames) {
    n1 = nchar(jobnames[1])
    n2 = nchar(jobnames[2])
    
    for (i in 1:min(n1,n2)) {
        if (substr(jobnames[1], i, i) != substr(jobnames[2], i, i)) {
            break
        }
    }
    if (i == 1)
        return (NA)
    else
        return (sub('_$', '', substr(jobnames[1], 1, i-1)))
}

process_data <- function(fn) {
    data = read.table(fn, header= FALSE, as.is= TRUE)
    colnames(data) = c('name', 'id', 'status', 'time_hr', 'mem_GB')
    data$time_hr = as.numeric(sub('s', '', data$time_hr)) / 3600 ## sec to hr

    idx = grepl('GB', data$mem_GB)
    data$mem_GB[idx] = sub('GB', '', data$mem_GB[idx])

    idxM = grepl('MB', data$mem_GB)
    data$mem_GB[idxM] = sub('MB', '', data$mem_GB[idxM])

    idxB = grepl('B', data$mem_GB)
    data$mem_GB[idxB] = sub('B', '', data$mem_GB[idxB])

    data$mem_GB = as.numeric(data$mem_GB)
    data$mem_GB[idxB] = data$mem_GB[idxB] / 1024
    data$mem_GB[idxM | idxB] = data$mem_GB[idxM | idxB] / 1024
    data$task = 'JG'
    data$task[grepl('_reduce', data$name)] = 'reduce'
    data$task[grepl('_melt', data$name)] = 'melt'
    data$task[grepl('_extract', data$name)] = 'extract'    
    return (data)
}


library(reshape)
library(ggplot2)
library(plyr)

data = process_data(fn)
n = sum(data$task == 'melt')
analysisName = get_analysis_name(data$name)

tmp = melt(data, id.vars= c('name', 'id', 'status', 'task'))
tmp$value = round(tmp$value, digits=1)

p = ggplot(tmp, aes(task, value, color=variable)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.3) + 
    facet_wrap( variable ~ task, scales= 'free' ) +
    scale_x_discrete(breaks = NULL) + 
    ggtitle(sprintf('%s, n=%d', analysisName, n))
##print(p)


ggsave(outfn, p)


## summarize
data$project = sub('_2016-$', '', regmatches(data$name, regexpr('\\S+2016-', data$name)))
data$date = regmatches(data$name, regexpr('2016-\\d+-\\d+', data$name))
data$nSample = sum(grepl('_melt', data$name))

if (grepl('Exome', data$project[1])) {
    data$type = 'Exome'
} else {
    data$type = 'WGS'
}



tmp = ddply(data, c('project', 'nSample', 'date', 'task', 'type'), summarise, njob = length(time_hr),
    max_time_hr= max(time_hr, na.rm=TRUE), total_time_hr = sum(time_hr, na.rm=TRUE), avg_time_hr = mean(time_hr, na.rm=TRUE), 
    max_mem_GB = max(mem_GB, na.rm=TRUE), avg_mem_GB = mean(mem_GB, na.rm=TRUE))

write.table(tmp, quote = FALSE, sep = '\t', file = outtxtfn, row.names=FALSE)
