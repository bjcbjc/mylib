

from os import popen, system, makedirs
import os.path
from sys import argv
from stat import S_ISDIR
import re
import paramiko
import traceback

#filterpattern is regular expression; \S* for all non space characters; ^ for beginning of the string (full path)

def main(host = 'harlem.nygenome.org', port = 22, username='bjchen', password='nyjupit#r)@18GC', remotepath='./', localpath='./', filterpattern=''):

    if username == 'pipeline': password = 'EnWhyG33Gn0m3!'

    transport = paramiko.Transport((host,port))

    transport.connect( username=username, password=password )

    sftp = paramiko.SFTPClient.from_transport(transport)

    localpath = os.path.abspath(localpath)
    try:
        # Download
        if remotepath[-1] != '/': remotepath = remotepath + '/'
        if localpath[-1] != '/': localpath = localpath + '/'

        if filterpattern == '':
            sftp.get(remotepath, localpath)
        else:
            pathpattern = re.compile(filterpattern)
            #filelist = map(lambda(l): remotepath+l, sftp.listdir(remotepath))
            filelist = sftp.listdir(remotepath)
            while not len(filelist) == 0:
                curfile = filelist.pop()
                if S_ISDIR( sftp.stat( remotepath + curfile ).st_mode ):
                    newcontent = map(lambda(l): curfile+'/'+l, sftp.listdir( remotepath + curfile ) )
                    filelist.extend( newcontent )
                else:
                    #if match the pattern, download
                    if pathpattern.match(remotepath + curfile) != None:
                        localdir = localpath + os.path.dirname(curfile)
                        if not os.path.exists( localdir ):
                            makedirs(localdir)
                        sftp.get(remotepath + curfile, localpath + curfile)
                        print 'download %s ==> %s'%(remotepath+curfile, localpath+curfile)

    # Upload
    #filepath = '/home/foo.jpg'
    #localpath = '/home/pony.jpg'
    #sftp.put(filepath, localpath)
    except:
        sftp.close()
        transport.close()
        print 'Error! Closed connection and exit.'
        traceback.print_exc()
    # Close
    sftp.close()
    transport.close()
    
    
    
if __name__ == '__main__':
    main( **dict( map(lambda(l): l.split('='), argv[1:]) ) )
