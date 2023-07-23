#!/usr/bin/env python3


def load_param_file(paramFilePath):
    from inspect import getsourcefile
    from os.path import abspath
    sourcePath = abspath(getsourcefile(lambda:0))
    sourceDir = sourcePath.rstrip("load_params.py")
    paramDict = dict()
    with open(paramFilePath) as F:
        for i in F:
            i = i.rstrip().split(",")
            if len(i)>1:
                paramDict[i[0]] = i[1]
    return paramDict

