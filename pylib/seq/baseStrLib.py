# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 13:54:53 2013

@author: bjchen
"""

import re
import numpy
from collections import Counter


"""
baseStrings are from samtools' mpileup
"""
class BaseString:
    indelOccurenceProg = re.compile('([\+-][0-9]+)')
#    cleanProg = re.compile('\^\S|\$')
    @staticmethod
    def matchIndels(baseString):
        """ return 
                indels: dict( indel, count_str)
                newBaseString: baseString without indels
        """            
        occurrences = BaseString.indelOccurenceProg.findall(baseString)
        indels = Counter()
        newBaseString = baseString
        for n in set(occurrences):
            indelprog = re.compile('\\%s[ATCGNatcgn]{%s}'%(n,n.strip('+-')))
            indelStrings = indelprog.findall(baseString)
            newBaseString = indelprog.sub('', newBaseString)
            indels.update( Counter(indelStrings) )
        return indels, newBaseString
    
    @staticmethod
    def splitBaseString(baseString):
        # indels do not have qual
        # so group indels with the previous base so we can match qual for it
#        if '^' in baseString or '$' in baseString:
#            baseString = cleanProg.sub('', baseString)
        occurrences = BaseString.indelOccurenceProg.finditer(baseString)
        startpos, endpos = [],[] # [start_pos, n]
        for m in occurrences: 
            token = m.group(0)
            n = int(token.strip('+-')) + len(token)
            s = m.start()
            startpos.append( s )
            endpos.append(s+n)
        i = 0
        tokens = []
        while i < len(baseString):
            if i in startpos:
                idx = startpos.index(i)
                indel = baseString[ startpos[idx]:endpos[idx] ]
                tokens[-1] = tokens[-1] + indel
                i = endpos[idx]
            else:
                tokens.append( baseString[i] )
                i = i + 1
        return tokens

            
""" 
baseCountStrings are strings from pileupCount.py
The following functions operate on baseCountStrings
"""  
class BaseCountString:
    indelOccurenceProg = re.compile('([\+-][0-9]+)')
    baseReProg = re.compile('([<>ATCGNatcgn])(\d*)')
    baseLabels = numpy.array([['>', 'A', 'T', 'C', 'G', 'N'],['<', 'a', 't', 'c', 'g', 'n']])
    
    @staticmethod
    def baseCountStringIndelCount(baseCountString):
        occurrences = BaseCountString.indelOccurenceProg.findall(baseCountString)
        indels = {}
        for n in set(occurrences):
            indelprog = re.compile('(\\%s[ATCGNatcgn]{%s})(\d*)'%(n,n.strip('+-')))
            indelCount = indelprog.findall(baseCountString)
            indels.update( {k:int(v) for k,v in indelCount} )
        return indels
    
    @staticmethod
    def baseCountStringRemoveIndels(baseCountString):
        occurrences = BaseCountString.indelOccurenceProg.findall(baseCountString)
        newBaseCountString = baseCountString
        for n in set(occurrences):
            indelprog = re.compile('\\%s[ATCGNatcgn]{%s}\d*'%(n,n.strip('+-')))
            newBaseCountString = indelprog.sub('', newBaseCountString)
        return newBaseCountString

    @staticmethod    
    def baseCountStringToMatrix(baseCountString):    
        matrix = numpy.zeros(BaseCountString.baseLabels.shape)
        tokens = BaseCountString.baseReProg.findall(baseCountString)
        for k, v in tokens:
            matrix[ BaseCountString.baseLabels == k ] = int(v)
        return matrix

    
    
