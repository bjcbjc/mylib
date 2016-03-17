# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 12:56:34 2013

@author: bjchen
"""

import string
import random

def randstr(size=10, chars=string.ascii_letters + string.digits):
    return ''.join(random.choice(chars) for x in xrange(size))