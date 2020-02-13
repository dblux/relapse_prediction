#!/usr/bin/env python3
# -*- coding: utf-8 -*-

with open("coreComplexes.txt") as f:
    content = f.readlines()
cnt = 0
fw = open("entrezId1.txt", 'w')
for line in content:
    data = line.strip().split("\t")
#     print(len(data), data[0], data[6])
    fw.write("%s\t%s\t%s\n" % (data[0],data[2],data[6]))
fw.close()